
dir_path = "SeqExpMatch/R/"
files = list.files(dir_path, pattern = "^inference_.*\\.R$", full.names = TRUE)

abstract_classes_to_asymp = c(
  "InferenceAbstractKKGEE",
  "InferenceAbstractKKGLMM",
  "InferenceAbstractKKClogitIVWC",
  "InferenceAbstractKKClogitCombinedLikelihood",
  "InferenceKKPassThrough",
  "InferenceMLEorKMforGLMs",
  "InferenceIncidGCompAbstract",
  "InferencePropGCompAbstract",
  "InferenceOrdinalGCompAbstract",
  "InferenceSurvivalStratCoxPHAbstract",
  "InferenceAbstractKKWilcoxBaseIVWC",
  "InferenceAbstractKKWilcoxRegrIVWC",
  "InferenceAbstractKKLWACoxIVWC",
  "InferenceAbstractKKLWACoxCombinedLikelihood",
  "InferenceAbstractKKClaytonCopulaIVWC",
  "InferenceAbstractKKClaytonCopulaCombinedLikelihood",
  "InferenceAbstractKKWeibullFrailtyIVWC",
  "InferenceAbstractKKWeibullFrailtyCombinedLikelihood",
  "InferenceAbstractKKPoissonCPoissonIVWC",
  "InferenceAbstractKKPoissonCPoissonCombinedLikelihood",
  "InferenceAbstractKKHurdlePoissonIVWC",
  "InferenceAbstractKKHurdlePoissonCombinedLikelihood",
  "InferenceCountHurdleNegBinAbstract",
  "InferenceCountZeroAugmentedPoissonAbstract",
  "InferencePropZeroOneInflatedBetaAbstract",
  "InferenceOrdinalPartialProportionalOddsAbstract",
  "InferenceIncidLogBinomialAbstract",
  "InferenceIncidConstrainedBinomialAbstract",
  "InferenceIncidBinomialIdentityAbstract"
)

inherit_tags = c(
  "#' @inherit InferenceRand methods",
  "#' @inherit InferenceBoot methods",
  "#' @inherit InferenceAsymp methods"
)

for (file in files) {
  lines = readLines(file, warn = FALSE)
  changed = FALSE
  
  # 1. Update inheritance for specific abstract classes
  for (cls in abstract_classes_to_asymp) {
    pattern = paste0("^", cls, "\\s*=\\s*R6::R6Class\\(\"", cls, "\",")
    idx = grep(pattern, lines)
    if (length(idx) > 0) {
      for (i in seq(idx + 1, min(idx + 15, length(lines)))) {
        if (grepl("inherit\\s*=", lines[i])) {
          if (!grepl("InferenceAsymp", lines[i])) {
            lines[i] = gsub("inherit\\s*=\\s*[A-Za-z0-9_]+\\s*([,)])$", "inherit = InferenceAsymp\\1", lines[i])
            changed = TRUE
          }
          break
        }
      }
    }
  }

  # 2. Update classes inheriting from Inference to InferenceAsymp
  r6_indices = grep("^[A-Za-z0-9_]+\\s*=\\s*R6::R6Class\\(", lines)
  for (idx in r6_indices) {
    for (i in seq(idx + 1, min(idx + 15, length(lines)))) {
       if (grepl("inherit\\s*=\\s*Inference\\s*([,)])$", lines[i])) {
         lines[i] = gsub("inherit\\s*=\\s*Inference\\s*([,)])$", "inherit = InferenceAsymp\\1", lines[i])
         changed = TRUE
         break
       }
    }
  }

  # 3. Add inherit tags for exported classes
  keep_going = TRUE
  while(keep_going) {
    keep_going = FALSE
    r6_indices = grep("^[A-Za-z0-9_]+\\s*=\\s*R6::R6Class\\(", lines)
    for (idx in r6_indices) {
      cls_name = gsub("^([A-Za-z0-9_]+)\\s*=\\s*R6::R6Class\\(.*", "\\1", lines[idx])
      if (startsWith(cls_name, "Inference")) {
        found_export = FALSE
        export_idx = -1
        # Search backwards from class definition
        for (i in seq(idx - 1, max(1, idx - 40))) {
          if (grepl("^#'\\s*@export", lines[i])) {
            found_export = TRUE
            export_idx = i
            break
          }
        }
        
        if (found_export) {
          # Check if tags are already there in the roxygen block before this class
          found_asymp = FALSE
          for (i in seq(idx - 1, max(1, idx - 40))) {
             if (grepl("@inherit InferenceAsymp methods", lines[i])) {
               found_asymp = TRUE
               break
             }
             if (grepl("^#'\\s*@description", lines[i])) {
                # We reached description, if we haven't found asymp yet, it's probably missing
                # But let's continue searching until start of block
             }
             if (!startsWith(lines[i], "#'")) {
                # End of roxygen block
                break
             }
          }
          
          if (!found_asymp) {
            lines = c(lines[1:(export_idx-1)], inherit_tags, lines[export_idx:length(lines)])
            changed = TRUE
            keep_going = TRUE
            break
          }
        }
      }
    }
  }
  
  if (changed) {
    writeLines(lines, file)
    cat("Updated", file, "\n")
  }
}
