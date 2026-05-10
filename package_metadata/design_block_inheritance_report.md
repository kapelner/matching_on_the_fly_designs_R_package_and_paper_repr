# DesignBlock Inheritance Report

## Question

Would it be possible to create a `DesignBlock` class that inherits `Design` and centralizes blocking data and methods such as strata information, `m`, and blocking helpers?

## Short Answer

Yes, but not as a trivial extraction.

A `DesignBlock` class is conceptually valid in this codebase because blocking is already treated as a capability that spans both fixed and sequential designs. However, the current implementation splits blocking logic across:

- `Design`, which already stores `m` and exposes several blocking-facing methods
- `DesignFixed`, which owns the main strata-to-block construction helper
- individual sequential blocking subclasses, which each implement their own strata-state logic

So this refactor is possible, but it is a medium-sized architectural change rather than a simple cleanup.

## Why The Idea Makes Sense

The current code already suggests that "blocking" is its own abstraction:

- `Design` already owns `private$m`
- `Design` already exposes `assert_blocking_design()`, `assert_equal_block_sizes()`, `set_m()`, `get_block_ids()`, and `summarize_blocks()`
- blocking is used by both fixed and sequential designs, not just one branch

That means the current hierarchy is already behaving as if there were a shared blocking layer, but that layer is only partially formalized.

In other words, `DesignBlock` would not introduce a new concept. It would make an existing concept explicit.

## What Prevents A Clean Extraction Today

### 1. Blocking support is detected by a hard-coded whitelist

`Design` currently decides whether a design "supports blocking" via `private$design_is_supported_blocking()`, which checks concrete classes by name.

That is a sign that the inheritance tree does not currently encode the abstraction correctly. A real `DesignBlock` base would allow this logic to become something much simpler, such as checking `is(self, "DesignBlock")`.

### 2. `get_strata_keys()` currently lives under `DesignFixed`

The main helper for turning covariates into block keys is currently in `DesignFixed`, not `Design`. That is awkward if blocking is meant to be shared across branches.

This helper also contains fixed-design-specific assumptions:

- it works off the full `Xraw`
- it uses greedy column accumulation toward `B_target`
- it bins continuous covariates using quantiles
- it enforces `exact_num_blocks` and `equal_block_sizes`

Some of that is generic blocking behavior, but some of it is really "fixed covariate-defined blocking".

### 3. Sequential blocking designs do not use the same representation

The sequential blocking classes do not currently rely on `get_strata_keys()`. Instead they use per-subject strata key builders and stateful block queues, for example:

- `DesignSeqOneByOneSPBR`
- `DesignSeqOneByOneRandomBlockSize`

So `DesignBlock` cannot just be "move `get_strata_keys()` upward". It would need a more abstract API that both fixed and sequential subclasses can satisfy.

### 4. Not all block-like designs share the same semantics

There are at least three different blocking patterns here:

- fixed blocking based on covariate strata
- fixed optimal blocks computed by an algorithm
- sequential blocking with within-stratum state

They all produce block identifiers, but they do not all construct them the same way.

That means a useful `DesignBlock` should own the common interface and state, not force one block-construction algorithm on every subclass.

## What Could Reasonably Move Into `DesignBlock`

These pieces are good candidates:

- `private$m`
- `public$assert_blocking_design()`
- `public$assert_equal_block_sizes()`
- `public$set_m()`
- `public$get_block_ids()`
- `public$summarize_blocks()`
- blocking-related shared fields such as `strata_cols`

Potentially also:

- shared bootstrap helpers for "within blocks" vs "resample whole blocks" where the only requirement is a valid block-id vector

## What Should Probably Stay Out Of A Generic `DesignBlock`

These are not uniformly shared enough to belong in the base:

- the greedy `B_target` logic from `DesignFixed`
- quantile binning for continuous covariates
- sequential per-stratum refill state such as permuted block queues
- optimal-block solver machinery

Those belong in subclasses or in narrower intermediate subclasses.

## A Better Shape For The Hierarchy

The most defensible version is:

```text
Design
├─ DesignBlock
│  ├─ DesignFixed
│  │  ├─ DesignFixedBlocking
│  │  ├─ DesignFixedOptimalBlocks
│  │  └─ DesignFixedBlockedCluster
│  └─ DesignSeqOneByOne
│     ├─ DesignSeqOneByOneSPBR
│     └─ DesignSeqOneByOneRandomBlockSize
└─ non-blocking designs
```

But there is an important catch: today `DesignFixed` and `DesignSeqOneByOne` both directly inherit `Design`.

Because R6 only gives you single inheritance, adopting `DesignBlock` would require deciding whether:

1. `DesignFixed` itself should inherit `DesignBlock`, or
2. only the blocking fixed subclasses should inherit `DesignBlock`, or
3. `DesignSeqOneByOne` itself should inherit `DesignBlock`, or
4. only the blocking sequential subclasses should inherit `DesignBlock`

The cleanest choice is not obvious.

## Recommended Refactor Strategy

I would not start by inserting `DesignBlock` immediately.

I would do this in two phases.

### Phase 1: Extract the interface without changing the tree much

First, make the blocking abstraction explicit at the method level:

- replace the class-name whitelist with a capability test
- introduce a small internal generic such as `private$compute_block_ids()` or `private$get_block_ids_impl()`
- let each blocking subclass provide its own implementation
- keep `get_block_ids()` in one shared place

This gets most of the architectural value without forcing a large inheritance rewrite.

### Phase 2: Add `DesignBlock` only if the shared surface is stable

Once the common API is real and the concrete subclasses fit it naturally, then move the shared blocking fields and methods into `DesignBlock`.

At that point, the new class will reflect an abstraction that has already been tested in the code, rather than guessed in advance.

## My Recommendation

`DesignBlock` is possible and directionally correct, but I would not implement it as the first move.

The current codebase does support the idea, but only partially:

- the shared blocking state already exists
- the shared blocking public API already mostly exists
- the block-construction logic is still fragmented by design family

So the right conclusion is:

`DesignBlock` is feasible, but only after first separating:

- common blocking state and user-facing methods
- family-specific block construction logic

## Bottom Line

Yes, this would be possible.

But the right version of `DesignBlock` is a thin abstraction over block identity and block-aware utilities, not a place to dump every blocking algorithm. If you try to move the fixed-design greedy strata builder, sequential block-state logic, and optimal-block solvers into one base class immediately, the new class will become a grab bag and the refactor will likely make the hierarchy less clear, not more.
