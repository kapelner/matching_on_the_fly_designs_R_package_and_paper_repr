import re
import os

def fix_duplicate_private(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    # Find all private = list( ... ) blocks
    # This regex is a bit fragile but we'll try to find the start and ends
    # Actually, let's just do it manually for the reported files since they are few.
    print(f"Manually check: {file_path}")

files = [
    "EDI/R/inference_incidence_KK_marginal.R",
    "EDI/R/inference_incidence_KK_clogit_plus_glmm.R",
    "EDI/R/inference_ordinal_cauchit.R",
    "EDI/R/inference_count_zero_inflated.R",
    "EDI/R/inference_incidence_gcomp.R",
    "EDI/R/inference_incidence_KK_clogit.R",
    "EDI/R/inference_count_hurdle_poisson_KK_ivwc_abstract.R",
    "EDI/R/inference_count_hurdle.R"
]

for f in files:
    fix_duplicate_private(f)
