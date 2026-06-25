# Speeding Up R Package Compilation — Full Change Record

All changes were made to the R-devel source tree at `/home/kapelner/R-devel`.
The build uses `MAKEFLAGS="-j <N>"` to signal the desired core count;
every parallelization reads that value via a shared helper.

---

## 0. Build Flags (`Makeconf`)

Already configured before these changes. Verify with:

```
grep -E "^CFLAGS|^FFLAGS|^LDFLAGS|openmp|OPENMP" src/main/Makeconf
```

Expected values:
```
CFLAGS = -O3 -march=native -mtune=native -fPIC -flto=auto -fprofile-generate -flto
FFLAGS = -O3 -march=native -mtune=native -fPIC -flto
R_OPENMP_CFLAGS = -fopenmp
ALL_CFLAGS = $(R_XTRA_CFLAGS) $(R_OPENMP_CFLAGS) $(MAIN_CFLAGS) $(CFLAGS)
```

---

## 1. `src/main/arithmetic.c` — OpenMP SIMD on real arithmetic

### 1a. Add OMP helper block after `#include <errno.h>`

**Before** (line ~59):
```c
#include <errno.h>
```

**After**:
```c
#include <errno.h>

#ifdef _OPENMP
# include <omp.h>
# include <ctype.h>

/* Minimum vector length before spawning parallel threads */
# define R_OMP_PARALLEL_MIN_N 100000L

/* Parse -j N / -jN / --jobs=N from MAKEFLAGS once, set omp thread count,
   and return the value for use in num_threads() clauses. */
static int R_omp_num_threads(void)
{
    static int nthreads = 0;
    if (nthreads > 0) return nthreads;
    const char *mf = getenv("MAKEFLAGS");
    if (mf) {
        const char *p = mf;
        while (*p) {
            if (p[0] == '-' && p[1] == 'j') {
                p += 2;
                while (*p == ' ') p++;
                if (isdigit((unsigned char)*p)) {
                    int n = atoi(p);
                    if (n > 0) { nthreads = n; omp_set_num_threads(n); return n; }
                }
                /* bare -j (unlimited) — fall through to default */
            } else if (strncmp(p, "--jobs=", 7) == 0) {
                p += 7;
                if (isdigit((unsigned char)*p)) {
                    int n = atoi(p);
                    if (n > 0) { nthreads = n; omp_set_num_threads(n); return n; }
                }
            } else {
                p++;
            }
        }
    }
    nthreads = omp_get_max_threads();
    return nthreads;
}
#endif /* _OPENMP */
```

### 1b. `real_unary()` MINUSOP (around line 841)

**Before**:
```c
    case MINUSOP:
        ans = NO_REFERENCES(s1) ? s1 : duplicate(s1);
        double *pa = REAL(ans);
        const double *px = REAL_RO(s1);
        n = XLENGTH(s1);
        for (i = 0; i < n; i++)
            pa[i] = -px[i];
        return ans;
```

**After**:
```c
    case MINUSOP:
        ans = NO_REFERENCES(s1) ? s1 : duplicate(s1);
        double *pa = REAL(ans);
        const double *px = REAL_RO(s1);
        n = XLENGTH(s1);
#ifdef _OPENMP
        (void) R_omp_num_threads();
#pragma omp parallel for simd schedule(static) if(n > R_OMP_PARALLEL_MIN_N)
#endif
        for (i = 0; i < n; i++)
            pa[i] = -px[i];
        return ans;
```

### 1c. `real_binary()` — equal-length REALSXP fast path (4 operators)

Applied identically for PLUSOP, MINUSOP, TIMESOP, DIVOP.
Example shown for PLUSOP (similar lines ~1031–1039):

**Before**:
```c
    else if (n1 == n2) {
        R_ITERATE_CHECK(NINTERRUPT, n, i, da[i] = dx[i] + dy[i];);
    }
```

**After**:
```c
    else if (n1 == n2) {
#ifdef _OPENMP
        (void) R_omp_num_threads();
#pragma omp parallel for simd schedule(static) if(n > R_OMP_PARALLEL_MIN_N)
        for (R_xlen_t k = 0; k < n; k++) da[k] = dx[k] + dy[k];
#else
        R_ITERATE_CHECK(NINTERRUPT, n, i, da[i] = dx[i] + dy[i];);
#endif
    }
```

Operators: `+` (PLUSOP), `-` (MINUSOP), `*` (TIMESOP), `/` (DIVOP).
POWOP is left unchanged (calls transcendental function, not safe for `parallel for`).

**Key constraint**: `R_ITERATE_CORE` uses `for (; i < n; ++i)` — no initializer —
so `#pragma omp simd` cannot be applied to it (OpenMP requires canonical form).
The fix is to bypass the macro and use a direct `for (R_xlen_t k = 0; ...)` loop.

---

## 2. `src/main/summary.c` — SIMD reduction on `rsum()`

**Before** (lines ~157–167):
```c
    ITERATE_BY_REGION(sx, x, i, nbatch, double, REAL, {
        for (R_xlen_t k = 0; k < nbatch; k++) {
            if (!narm || !ISNAN(x[k])) {
                if(!updated) updated = true;
                s += x[k];
            }
        }
    });
```

**After**:
```c
    ITERATE_BY_REGION(sx, x, i, nbatch, double, REAL, {
#if defined(_OPENMP) && HAVE_OPENMP_SIMDRED
        _Pragma("omp simd reduction(+:s)")
#endif
        for (R_xlen_t k = 0; k < nbatch; k++) {
            if (!narm || !ISNAN(x[k])) {
                if(!updated) updated = true;
                s += x[k];
            }
        }
    });
```

`HAVE_OPENMP_SIMDRED` is a configure-time flag that guards SIMD reduction
support (not all compilers handle `reduction` on `LDOUBLE`).

---

## 3. `src/main/relop.c` — SIMD on numeric comparison loops

### Add `ISNA_INT` macro and rewrite `NR_HELPER` (around line 424)

**Before**:
```c
#define ISNA_INT(x) x == NA_INTEGER

#define NR_HELPER(OP, type1, ACCESSOR1, ISNA1, type2, ACCESSOR2, ISNA2) do { \
    type1 x1, *px1 = ACCESSOR1(s1);                                           \
    type2 x2, *px2 = ACCESSOR2(s2);                                           \
    int *pa = LOGICAL(ans);                                                    \
    MOD_ITERATE2(n, n1, n2, i, i1, i2, {                                      \
        x1 = px1[i1]; x2 = px2[i2];                                           \
        if (ISNA1(x1) || ISNA2(x2)) pa[i] = NA_LOGICAL;                      \
        else pa[i] = (x1 OP x2);                                              \
    });                                                                        \
} while (0)
```

**After**:
```c
#define ISNA_INT(x) x == NA_INTEGER

#define NR_HELPER(OP, type1, ACCESSOR1, ISNA1, type2, ACCESSOR2, ISNA2) do { \
    type1 x1, *px1 = ACCESSOR1(s1);                                           \
    type2 x2, *px2 = ACCESSOR2(s2);                                           \
    int *pa = LOGICAL(ans);                                                    \
    if (n1 == n2) {                                                            \
        R_xlen_t _k;                                                           \
        _Pragma("omp simd")                                                    \
        for (_k = 0; _k < n; _k++) {                                          \
            type1 _v1 = px1[_k];                                               \
            type2 _v2 = px2[_k];                                               \
            pa[_k] = (ISNA1(_v1) || ISNA2(_v2)) ?                             \
                NA_LOGICAL : (_v1 OP _v2);                                     \
        }                                                                      \
    } else {                                                                   \
        MOD_ITERATE2(n, n1, n2, i, i1, i2, {                                 \
            x1 = px1[i1]; x2 = px2[i2];                                       \
            if (ISNA1(x1) || ISNA2(x2)) pa[i] = NA_LOGICAL;                  \
            else pa[i] = (x1 OP x2);                                          \
        });                                                                    \
    }                                                                          \
} while (0)
```

Covers all 6 operators (==, !=, <, >, <=, >=) via the `NUMERIC_RELOP` macro.

---

## 4. `src/main/logic.c` — SIMD on logical operations

### 4a. `lunary()` — NOT operator

In the `LGLSXP`, `INTSXP`, `REALSXP`, and `RAWSXP` branches, add
`_Pragma("omp simd")` before each `for` loop. Example:

**Before**:
```c
    for (i = 0; i < n; i++) {
        int x = LOGICAL_ELT(s1, i);
        LOGICAL(ans)[i] = (x == NA_LOGICAL) ? NA_LOGICAL : !x;
    }
```

**After**:
```c
    _Pragma("omp simd")
    for (i = 0; i < n; i++) {
        int x = LOGICAL_ELT(s1, i);
        LOGICAL(ans)[i] = (x == NA_LOGICAL) ? NA_LOGICAL : !x;
    }
```

### 4b. `binaryLogic()` — AND / OR

**Before**:
```c
    switch (code) {
    case 1:  /* & : AND */
        MOD_ITERATE2(n, n1, n2, i, i1, i2, {
            x1 = px1[i1]; x2 = px2[i2];
            if (x1 == 0 || x2 == 0)                    pa[i] = 0;
            else if (x1 == NA_LOGICAL || x2 == NA_LOGICAL) pa[i] = NA_LOGICAL;
            else                                        pa[i] = 1;
        });
        break;
    case 2:  /* | : OR */
        MOD_ITERATE2(n, n1, n2, i, i1, i2, { ... });
        break;
    }
```

**After**:
```c
    switch (code) {
    case 1:  /* & : AND */
        if (n1 == n2) {
            R_xlen_t _k;
            _Pragma("omp simd")
            for (_k = 0; _k < n; _k++) {
                int _a = px1[_k], _b = px2[_k];
                if (_a == 0 || _b == 0)                        pa[_k] = 0;
                else if (_a == NA_LOGICAL || _b == NA_LOGICAL) pa[_k] = NA_LOGICAL;
                else                                           pa[_k] = 1;
            }
        } else {
            MOD_ITERATE2(n, n1, n2, i, i1, i2, { /* original */ });
        }
        break;
    case 2:  /* | : OR */
        if (n1 == n2) {
            R_xlen_t _k;
            _Pragma("omp simd")
            for (_k = 0; _k < n; _k++) {
                int _a = px1[_k], _b = px2[_k];
                if ((_a != NA_LOGICAL && _a) || (_b != NA_LOGICAL && _b)) pa[_k] = 1;
                else if (_a == 0 && _b == 0)                              pa[_k] = 0;
                else                                                       pa[_k] = NA_LOGICAL;
            }
        } else {
            MOD_ITERATE2(n, n1, n2, i, i1, i2, { /* original */ });
        }
        break;
    }
```

### 4c. `binaryLogic2()` — RAW AND / OR

**Before**:
```c
    switch (code) {
    case 1:
        MOD_ITERATE2(n, n1, n2, i, i1, i2, {
            x1 = pr1[i1]; x2 = pr2[i2]; pa[i] = x1 & x2;
        });
        break;
    case 2:
        MOD_ITERATE2(n, n1, n2, i, i1, i2, {
            x1 = pr1[i1]; x2 = pr2[i2]; pa[i] = x1 | x2;
        });
        break;
    }
```

**After** (pointer hoisting + n1==n2 fast path added before each MOD_ITERATE2):
```c
    const Rbyte *pr1 = RAW_RO(s1);
    const Rbyte *pr2 = RAW_RO(s2);
    Rbyte *pa = RAW(ans);
    switch (code) {
    case 1:
        if (n1 == n2) {
            R_xlen_t _k;
            _Pragma("omp simd")
            for (_k = 0; _k < n; _k++) pa[_k] = pr1[_k] & pr2[_k];
        } else {
            MOD_ITERATE2(n, n1, n2, i, i1, i2, {
                x1 = pr1[i1]; x2 = pr2[i2]; pa[i] = x1 & x2;
            });
        }
        break;
    case 2:
        if (n1 == n2) {
            R_xlen_t _k;
            _Pragma("omp simd")
            for (_k = 0; _k < n; _k++) pa[_k] = pr1[_k] | pr2[_k];
        } else {
            MOD_ITERATE2(n, n1, n2, i, i1, i2, {
                x1 = pr1[i1]; x2 = pr2[i2]; pa[i] = x1 | x2;
            });
        }
        break;
    }
```

---

## 5. `src/main/eval.c` — NOT modified

`bcEval_loop()` fast-path scalar ops are per-element stack operations with no
arrays to vectorize. No SIMD applies. Confirmed zero `omp`/`simd`/`pragma`
references added.

---

## 6. `.makeflags_ncores()` — shared R helper

Parses `-j N`, `-jN`, `--jobs=N` from `MAKEFLAGS`. Returns an integer.
This function is defined in **three** places (each file that needs it
independently — two because of how the tools build pipes files directly to R):

- `src/library/compiler/R/cmp.R` (before `cmpframe`)
- `src/library/tools/R/install.R` (before `.convertRdfiles`)
- `src/library/tools/R/makeLazyLoad.R` (MUST be here — the tools byte-compile
  step pipes `makeLazyLoad.R` directly to R via `cat makeLazyLoad.R | R`,
  so functions from other files in the package are not available)

```r
.makeflags_ncores <- function() {
    mf <- Sys.getenv("MAKEFLAGS", unset = "")
    if (!nzchar(mf)) return(1L)
    m <- regmatches(mf, regexpr("-j *([0-9]+)", mf, perl = TRUE))
    if (length(m) && nzchar(m)) {
        n <- suppressWarnings(as.integer(sub("-j *", "", m)))
        if (!is.na(n) && n > 0L) return(n)
    }
    m <- regmatches(mf, regexpr("--jobs=([0-9]+)", mf, perl = TRUE))
    if (length(m) && nzchar(m)) {
        n <- suppressWarnings(as.integer(sub("--jobs=", "", m)))
        if (!is.na(n) && n > 0L) return(n)
    }
    1L
}
```

---

## 7. `src/library/compiler/R/cmp.R` — parallel `cmpframe()`

`cmpframe()` is called by `cmpfile()` / `cmplib()` (not the main `R CMD INSTALL`
byte-compile path — see §8 for that).

**Before**:
```r
cmpframe <- function(inpos, file) {
    ...
    all_names <- ls(pos = inpos, all.names = TRUE)
    closures  <- Filter(function(f) typeof(get(f, pos=inpos))=="closure", all_names)
    for (f in closures) {
        def <- get(f, pos = inpos)
        cat(gettextf("compiling '%s'", f), "\n", sep="")
        assign(f, cmpfun(def), pos = outpos)
    }
    ...
}
```

**After**:
```r
cmpframe <- function(inpos, file) {
    ...
    all_names <- ls(pos = inpos, all.names = TRUE)
    closures  <- Filter(function(f) typeof(get(f, pos = inpos)) == "closure",
                        all_names)
    ncores <- .makeflags_ncores()

    compile_one <- function(f) {
        def <- get(f, pos = inpos)
        cat(gettextf("compiling '%s'", f), "\n", sep = "")
        cmpfun(def)
    }

    compiled <- if (ncores > 1L &&
                    requireNamespace("parallel", quietly = TRUE)) {
        parallel::mclapply(closures, compile_one, mc.cores = ncores)
    } else {
        lapply(closures, compile_one)
    }

    for (j in seq_along(closures))
        assign(closures[[j]], compiled[[j]], pos = outpos)
    ...
}
```

---

## 8. `src/library/tools/R/makeLazyLoad.R` — parallel `code2LazyLoadDB()`

This is the **main `R CMD INSTALL` byte-compile path**:
`install.R` → sets `compilePKGS(1L)` → `makeLazyLoading()` → `code2LazyLoadDB()`.

**Critical note**: The tools build system compiles `tools` itself by running:
```sh
(cat src/library/tools/R/makeLazyLoad.R; echo 'makeLazyLoading("tools")') | R
```
This means `makeLazyLoad.R` is **sourced in isolation** — functions from other
files in the `tools` package (including `install.R`) are NOT available.
Therefore `.makeflags_ncores()` MUST be defined at the top of `makeLazyLoad.R`.

**Before**:
```r
code2LazyLoadDB <-
    function(package, lib.loc = NULL, ...)
{
    ...
    ns <- suppressPackageStartupMessages(loadNamespace(
              package = package, ..., partial = TRUE))
    makeLazyLoadDB(ns, dbbase, compress = compress,
                   set.install.dir = set.install.dir)
    ...
}
```

**After** (full replacement — `.makeflags_ncores` defined above it):
```r
.makeflags_ncores <- function() { ... }   # same as §6

code2LazyLoadDB <-
    function(package, lib.loc = NULL, ...)
{
    ...
    ## Disable auto-compilation during loading so we can compile in parallel.
    old_pkgs <- compiler::compilePKGS(0L)
    on.exit(compiler::compilePKGS(old_pkgs), add = TRUE)

    ns <- suppressPackageStartupMessages(loadNamespace(
              package = package, ..., partial = TRUE))

    if (old_pkgs != 0L) {
        all_names <- ls(ns, all.names = TRUE)
        closures  <- all_names[vapply(all_names, function(f)
            typeof(get(f, envir = ns, inherits = FALSE)) == "closure",
            logical(1L))]

        ## "Simple" closures have ns as their environment — safe to compile in
        ## forks. "Complex" ones (defined in local() etc.) compiled sequentially.
        simple  <- closures[vapply(closures, function(f)
            identical(environment(get(f, envir = ns, inherits = FALSE)), ns),
            logical(1L))]
        complex <- closures[!closures %in% simple]

        ns_assign <- function(nm, fn) {
            locked <- bindingIsLocked(nm, ns)
            if (locked) unlockBinding(nm, ns)
            assign(nm, fn, envir = ns)
            if (locked) lockBinding(nm, ns)
        }

        ncores <- .makeflags_ncores()
        if (ncores > 1L && length(simple) > 0L &&
                requireNamespace("parallel", quietly = TRUE)) {
            cfns <- parallel::mclapply(simple, function(f)
                compiler::cmpfun(get(f, envir = ns, inherits = FALSE)),
                mc.cores = ncores)
            ## Each fork compiled against a copy of ns. Fix the environment
            ## to point to the parent's ns before saving.
            for (i in seq_along(simple)) {
                fn <- cfns[[i]]
                environment(fn) <- ns
                ns_assign(simple[[i]], fn)
            }
        } else {
            for (f in simple)
                ns_assign(f, compiler::cmpfun(get(f, envir = ns, inherits = FALSE)))
        }

        for (f in complex)
            ns_assign(f, compiler::cmpfun(get(f, envir = ns, inherits = FALSE)))
    }

    makeLazyLoadDB(ns, dbbase, compress = compress,
                   set.install.dir = set.install.dir)
    ...
}
```

**Why `environment(fn) <- ns`**: `mclapply` forks the process. Each fork's
`ns` is a copy. After the fork returns the compiled function, its `environment`
still points to the fork's copy of `ns` (now gone). Reassigning it to the
parent's `ns` fixes variable lookup. This is safe for top-level closures
(whose original environment IS `ns`). Closures with sub-environments (from
`local()`) are compiled sequentially in the `complex` path.

**Why `ns_assign` with `unlockBinding`**: When byte-compiling the `compiler`
package itself, its namespace is already loaded and sealed before
`code2LazyLoadDB` is called (R bootstraps `compiler` early). In that case
`loadNamespace(..., partial = TRUE)` silently returns the already-sealed
namespace, making plain `assign(..., envir = ns)` fail with
`"cannot change value of locked binding"`. `ns_assign` unlocks each binding
before writing and re-locks it after, handling both sealed and unsealed
namespaces uniformly. The same helper is used for both the parallel and
sequential paths.

---

## 9. `src/library/tools/R/install.R` — parallel `.convertRdfiles()` and `** inst`

### 9a. `.makeflags_ncores()` — add before `.convertRdfiles` (~line 3177)

Same implementation as §6.

### 9b. `.convertRdfiles()` — parallel Rd-to-HTML/latex/example conversion

The "installing help indices" sub-step of `** help`.

**Before**: a `for (nf in files)` loop with local mutable state (`shown`,
`.messages`) converted each Rd file sequentially.

**After**: refactored to `convert_one_rd(nf)` closure returning
`list(bf, did, msgs)`, dispatched via `mclapply`, output printed sequentially
after all workers complete. Dead-code closures (`showtype`, `.whandler`,
`.ehandler`, `.convert`) removed.

Key structure:
```r
convert_one_rd <- function(nf) {
    msgs <- character()
    wh <- function(e) { msgs <<- c(msgs, ...); tryInvokeRestart("muffleWarning") }
    conv <- function(expr) withCallingHandlers(tryCatch(expr, error=...), warning=wh)
    Rd <- db[[nf]]; bf <- sub("\\.[Rr]d$", "", basename(nf))
    did <- character()
    if ("html" %in% types) { ff <- ...; if (needs update) { conv(Rd2HTML(...)); did <- c(did,"html") } }
    if ("latex" %in% types) { ... }
    if ("example" %in% types) { ... }
    list(bf = bf, did = did, msgs = unique(msgs))
}

ncores <- .makeflags_ncores()
results <- if (ncores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(files, convert_one_rd, mc.cores = ncores)
} else {
    lapply(files, convert_one_rd)
}

for (res in results) {   # print in collation order
    if (length(res$did)) { cat("   ", res$bf, ...); cat("\n") }
    if (length(res$msgs)) writeLines(res$msgs)
}
```

### 9c. `** inst` — parallel file copy

**Before**:
```r
i2_files <- gsub("^inst", quote_replacement(instdir), i_files)
file.copy(i_files, i2_files)
```

**After**:
```r
i2_files <- gsub("^inst", quote_replacement(instdir), i_files)
ncores_inst <- .makeflags_ncores()
if (ncores_inst > 1L && length(i_files) > ncores_inst &&
        requireNamespace("parallel", quietly = TRUE)) {
    chunks <- parallel::splitIndices(length(i_files), ncores_inst)
    parallel::mclapply(chunks, function(idx)
        file.copy(i_files[idx], i2_files[idx]),
        mc.cores = ncores_inst)
} else {
    file.copy(i_files, i2_files)
}
```

Directories are still created sequentially before this (they must exist first).
`Sys.chmod` stays sequential (fast vectorized call).

---

## 10. `src/library/tools/R/admin.R` — parallel indices and R-file parsing

### 10a. `.vinstall_package_indices()` — parallel across base packages

Called during `make install` for the base package set (not per user package).

**Before**:
```r
for (p in pkgs)
    .install_package_indices(file.path(src_dir, p), file.path(out_dir, p))
```

**After**:
```r
pkgs <- unlist(strsplit(packages, "[[:space:]]+"))
ncores <- .makeflags_ncores()
if (ncores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(pkgs, function(p)
        .install_package_indices(file.path(src_dir, p),
                                 file.path(out_dir,  p)),
        mc.cores = ncores)
} else {
    for (p in pkgs)
        .install_package_indices(file.path(src_dir, p),
                                 file.path(out_dir,  p))
}
```

### 10b. `.install_package_code_files()` — parallel R-file parsing (`** R` step)

**Branch 1 (package has `Encoding:` in DESCRIPTION)**:

**Before**: sequential `for(f in codeFiles)` loop doing read+iconv+parse+write.

**After**: per-file processing extracted to `process_one_enc(f)` returning
`list(f, line1, tmp, bad_utf8, bad2, err)`, dispatched via `mclapply`,
then results written sequentially to `con` in collation order.

```r
process_one_enc <- function(f) {
    lines <- readLines(f, warn = FALSE)
    bad_utf8 <- enc == "UTF-8" && any(!validUTF8(lines))
    tmp <- iconv(lines, from = enc, to = "")
    bad <- which(is.na(tmp))
    if (length(bad)) tmp <- iconv(lines, from = enc, to = "", sub = "byte")
    comm <- grep("^[^#'\"]*#", lines[bad], invert = TRUE, useBytes = TRUE)
    bad2 <- bad[comm]
    line1 <- paste0("#line 1 \"", f, "\"")
    err <- tryCatch(testParse(text = c(line1, tmp)), error = function(e) e)
    list(f = f, line1 = line1, tmp = tmp, bad_utf8 = bad_utf8, bad2 = bad2, err = err)
}
results <- if (ncores > 1L && length(codeFiles) > 1L && ...) {
    parallel::mclapply(codeFiles, process_one_enc, mc.cores = ncores)
} else lapply(codeFiles, process_one_enc)

con <- file(outFile, "a"); on.exit(close(con))
for (res in results) {
    if (res$bad_utf8) warning(...)
    if (length(res$bad2)) warning(...)
    if (inherits(res$err, "error")) stop(conditionMessage(res$err), ...)
    writeLines(res$line1, con); writeLines(res$tmp, con)
}
```

**Branch 2 (no explicit encoding — common case)**:

**Before**:
```r
lapply(codeFiles, testParse)
if (!all(.file_append_ensuring_LFs(outFile, codeFiles)))
    stop("unable to write code files")
```

**After**:
```r
if (ncores > 1L && length(codeFiles) > 1L &&
        requireNamespace("parallel", quietly = TRUE)) {
    results <- parallel::mclapply(codeFiles, function(f)
        tryCatch(testParse(f), error = function(e) e),
        mc.cores = ncores)
    errs <- Filter(function(x) inherits(x, "error"), results)
    if (length(errs))
        stop(conditionMessage(errs[[1L]]), domain = NA, call. = FALSE)
} else {
    lapply(codeFiles, testParse)
}
if (!all(.file_append_ensuring_LFs(outFile, codeFiles)))
    stop("unable to write code files")
```

---

## 11. `src/library/tools/R/Rd.R` — parallel Rd parsing (`** help` step)

### `.build_Rd_db()` — parallel `prepare_Rd()` per file

**Before** (~line 408):
```r
if (length(files)) {
    db1 <- lapply(files, .fetch_Rd_object, stages)
    names(db1) <- names(files)
    db <- c(db, db1)
}
```

**After**:
```r
if (length(files)) {
    ncores <- .makeflags_ncores()
    db1 <- if (ncores > 1L &&
               requireNamespace("parallel", quietly = TRUE)) {
        parallel::mclapply(files, .fetch_Rd_object, stages,
                           mc.cores = ncores)
    } else {
        lapply(files, .fetch_Rd_object, stages)
    }
    names(db1) <- names(files)
    db <- c(db, db1)
}
```

Each `prepare_Rd()` call parses one `.Rd` file independently (CPU-bound).
This parallelizes the `** help` step directly.

---

## 12. Build and install

**Critical**: `R CMD INSTALL` uses whichever `R` is first on `PATH`. If the
source tree was configured with a prefix other than `/usr/local` (check with
`grep ^S\[.prefix.\] config.status`), the `sudo make install` will NOT update
`/usr/local/bin/R`. Reconfigure with the correct prefix before installing:

```bash
# If the configured prefix is not /usr/local, reconfigure first:
./configure --prefix=/usr/local --enable-R-shlib --with-blas --with-lapack
# (add any other flags that were in the original configure invocation)
```

```bash
# Full build (all changes):
MAKEFLAGS="-j48" make -j48

# Install (goes to the prefix set by configure — default /usr/local):
sudo make install -j48

# If only tools R files changed (faster):
MAKEFLAGS="-j48" make -j48 -C src/library/tools
sudo cp library/tools/R/tools.rdb \
        library/tools/R/tools.rdx \
        library/tools/R/tools \
        /usr/local/lib/R/library/tools/R/
```

---

## 13. Steps NOT parallelized (and why)

| Step | Reason |
|---|---|
| `** byte-compile` (bcEval_loop) | Per-scalar stack ops, no arrays |
| `testing if installed package keeps a record of temporary installation path` | Single subprocess check on one package; no loop |
| `** building package indices` (single package) | Sequential pipeline: load db → `Rd_contents()` → 4 saveRDS; each step depends on previous; all < 0.5s total |

---

## 14. Verification

```r
/usr/local/bin/R --vanilla -q -e '
f <- getFromNamespace("code2LazyLoadDB",         "tools"); b <- deparse(body(f))
f2 <- getFromNamespace(".convertRdfiles",        "tools"); b2 <- deparse(body(f2))
f3 <- getFromNamespace(".vinstall_package_indices","tools"); b3 <- deparse(body(f3))
f4 <- getFromNamespace(".build_Rd_db",           "tools"); b4 <- deparse(body(f4))
f5 <- getFromNamespace(".install_package_code_files","tools"); b5 <- deparse(body(f5))
f6 <- getFromNamespace("cmpframe",               "compiler"); b6 <- deparse(body(f6))
cat("code2LazyLoadDB parallel:         ", any(grepl("mclapply", b)),  "\n")
cat("convertRdfiles parallel:          ", any(grepl("mclapply", b2)), "\n")
cat("vinstall_package_indices parallel:", any(grepl("mclapply", b3)), "\n")
cat("build_Rd_db parallel:             ", any(grepl("mclapply", b4)), "\n")
cat("install_package_code_files parallel:", any(grepl("mclapply", b5)), "\n")
cat("cmpframe parallel:                ", any(grepl("mclapply", b6)), "\n")
'
```

All six should print `TRUE`.
