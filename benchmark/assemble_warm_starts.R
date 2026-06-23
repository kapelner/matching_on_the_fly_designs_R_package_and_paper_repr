
pkgload::load_all("EDI", quiet = TRUE)
library(data.table)

N_VALS = c(100L, 200L, 500L, 1000L)

OP_FOR_COL = list(
    Rand     = "rand",
    Boot     = "non_param_boot",
    Bayesian = "bayesian_boot",
    JK       = "jackknife",
    PB       = "param_boot"
)

cell_style = function(val) {
    if (val == "(D)") return('background:#f0f0f0; color:#999; font-style:italic;')
    base = sub(" \\(D\\)$", "", val)
    if (base == "N/S")   return('background:#f5f5f5; color:#777;')
    if (grepl("^same",   base)) return('background:#f5f5f5; color:#777;')
    if (grepl("^< 2%",   base)) return('background:#fffdf0; color:#8a6d00;')
    if (grepl("^-",      base)) return('background:#fdecec; color:#9f1d1d;')
    'background:#eaf6ec; color:#1b6b32;'
}

annotate_D = function(val, cls, op, n) {
    if (is.na(val) || val == "N/S" || val == "(D)") return(val)
    enabled = tryCatch(
        edi_warm_start_dispatch_policy(cls, op, n = n),
        error = function(e) TRUE
    )
    if (!isTRUE(enabled) && !grepl("\\(D\\)", val)) paste0(val, " (D)") else val
}

make_table_html = function(merged_n, n_val) {
    lines = c(
        '<table border="1" style="border-collapse: collapse; width: 100%;">',
        '<thead>',
        '<tr>',
        '<th style="text-align:left; padding:4px;">Path</th>',
        '<th style="text-align:center; padding:4px;">Random-<br>ization</th>',
        '<th style="text-align:center; padding:4px;">Nonparam<br> Bootstrap</th>',
        '<th style="text-align:center; padding:4px;">Bayesian Bootstrap</th>',
        '<th style="text-align:center; padding:4px;">Jack-<br>knife</th>',
        '<th style="text-align:center; padding:4px;">Param. Bootstrap</th>',
        '</tr>',
        '</thead>',
        '<tbody>'
    )
    for (i in seq_len(nrow(merged_n))) {
        r = merged_n[i]
        vals = list(
            Rand     = annotate_D(r$Rand,     r$Path, "rand",          n_val),
            Boot     = annotate_D(r$Boot,     r$Path, "non_param_boot", n_val),
            Bayesian = annotate_D(r$Bayesian, r$Path, "bayesian_boot",  n_val),
            JK       = annotate_D(r$JK,       r$Path, "jackknife",      n_val),
            PB       = annotate_D(r$PB,       r$Path, "param_boot",     n_val)
        )
        lines = c(lines, '<tr>')
        lines = c(lines, sprintf('<td style="padding:4px; font-family:monospace;">%s</td>', r$Path))
        for (col in c("Rand", "Boot", "Bayesian", "JK", "PB")) {
            v = if (is.na(vals[[col]])) "N/S" else vals[[col]]
            s = cell_style(v)
            lines = c(lines, sprintf('<td style="text-align:center; padding:4px; %s">%s</td>', s, v))
        }
        lines = c(lines, '</tr>')
    }
    c(lines, '</tbody>', '</table>')
}

all_tables = list()
for (n_val in N_VALS) {
    f_csv  = sprintf("warm_starts_final_results_n%d.csv", n_val)
    bj_csv = sprintf("warm_starts_bayes_jackknife_results_n%d.csv", n_val)

    if (!file.exists(f_csv) || !file.exists(bj_csv)) {
        cat(sprintf("Missing CSVs for n=%d, skipping.\n", n_val))
        next
    }

    f_dt  = fread(f_csv,  col.names = c("Path", "Rand", "Boot", "JK_f", "PB"))
    bj_dt = fread(bj_csv, col.names = c("Path", "Bayesian", "JK"))

    merged = merge(f_dt, bj_dt, by = "Path", all = TRUE)
    merged[, JK_f := NULL]
    merged[is.na(Rand),     Rand     := "N/S"]
    merged[is.na(Boot),     Boot     := "N/S"]
    merged[is.na(Bayesian), Bayesian := "N/S"]
    merged[is.na(JK),       JK       := "N/S"]
    merged[is.na(PB),       PB       := "N/S"]
    setorder(merged, Path)

    all_tables[[as.character(n_val)]] = make_table_html(merged, n_val)
    cat(sprintf("Built table for n=%d (%d paths).\n", n_val, nrow(merged)))
}

md_file = "package_metadata/warm_starts.md"
md = readLines(md_file, warn = FALSE)

for (n_val in N_VALS) {
    tbl = all_tables[[as.character(n_val)]]
    if (is.null(tbl)) next

    section_header = sprintf("### n = %d", n_val)
    start_idx = which(md == section_header)
    if (length(start_idx) == 0L) { cat(sprintf("Section '%s' not found.\n", section_header)); next }
    start_idx = start_idx[1L]

    table_start = start_idx + which(grepl("^<table", md[(start_idx + 1L):length(md)]))[1L]
    table_end   = table_start + which(grepl("^</table>$", md[(table_start + 1L):length(md)]))[1L]

    md = c(md[seq_len(table_start - 1L)], tbl, md[(table_end + 1L):length(md)])
    cat(sprintf("Replaced n=%d table in warm_starts.md.\n", n_val))
}

# Update "12 timing repetitions" to "100 timing repetitions"
md = gsub("12 timing repetitions", "100 timing repetitions", md, fixed = TRUE)

writeLines(md, md_file)
cat("Done. warm_starts.md updated.\n")
