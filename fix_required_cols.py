import re
import os
import glob

def process_file(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    # Case 1: required_cols = 1L in calls, where we likely want to match "w" or "treatment" or "beta_T"
    # We need to be careful about the variable name.
    
    # Generic replacement for common patterns in this codebase
    new_content = content
    
    # Files using "w" as treatment column in design matrix
    if "KK" in file_path or "gee" in file_path or "glmm" in file_path or "clogit" in file_path:
        new_content = re.sub(r'required_cols\s*=\s*1L', 'required_cols = match("w", colnames(X_full))', new_content)
    elif "logit.R" in file_path or "log_binomial" in file_path or "modified_poisson" in file_path or "risk_diff" in file_path or "binomial_identity" in file_path or "robust_regr" in file_path:
         new_content = re.sub(r'required_cols\s*=\s*1L', 'required_cols = match("treatment", colnames(X_full))', new_content)
    elif "one_lik" in file_path:
         new_content = re.sub(r'required_cols\s*=\s*j_beta_T', 'required_cols = j_beta_T', new_content) # already fine
    
    # Specific fix for abstract method signature
    if "inference_all_abstract.R" in file_path:
        new_content = re.sub(r'required_cols\s*=\s*1L', 'required_cols = 2L', new_content)

    if new_content != content:
        with open(file_path, 'w') as f:
            f.write(new_content)
        print(f"Updated {file_path}")

files = glob.glob('EDI/R/inference_*.R')
for f in files:
    process_file(f)
