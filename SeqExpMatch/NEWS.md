# SeqExpMatch 0.1.1

Resubmission to CRAN after archival on 2025-09-23 ("check problems were not corrected despite reminders"). This release fixes the outstanding `R CMD check` NOTEs:

* Fixed malformed Rd markup in `SeqDesignInference$compute_confidence_interval()` documentation: literal braces in inline formulas (e.g. `t_{alpha/2, ...}`) were being parsed as Rd markup delimiters instead of text, producing a NOTE. Braces are now escaped.
* Updated the arXiv citation in `DESCRIPTION` from the deprecated `<arXiv:2010.05980>` form to the published paper's DOI form `<doi:10.1111/biom.13561>` (Kapelner and Krieger, *Biometrics*, 2021), per current CRAN citation requirements.
* Added an `Authors@R` field to `DESCRIPTION` (previously used the older `Author`/`Maintainer` fields only).
* Added `inst/CITATION` so `citation("SeqExpMatch")` points readers directly to the published paper (https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13561) instead of an auto-generated package citation.

No user-facing functionality changed in this release.

# SeqExpMatch 0.1.0

Initial CRAN release. Provides sequential two-arm experimental designs (completely randomized, balanced completely randomized, Efron's biased coin, Atkinson's covariate-adjusted biased coin, and Kapelner and Krieger's covariate-adjusted matching-on-the-fly designs including the CARA variants) along with estimation, testing, and confidence interval procedures for the resulting data.
