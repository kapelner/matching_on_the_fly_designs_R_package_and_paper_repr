pacman::p_load(data.table, Rcpp, stringr)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/res"))

all_data_files = str_subset(dir(), regex("result_"))
load(all_data_files[1])
res = data.frame(matrix(NA_real_, nrow = max(as.numeric(str_match(all_data_files, regex("_(\\d*)"))[,2])), ncol = length(res_row)))

res_row$weights = NA_character_
list_to_turn_to_df = lapply(res_row, function(x){rep(x, length(all_data_files))})
list_to_turn_to_df$filename = rep(NA_character_, length(all_data_files))
res = data.table(data.frame(list_to_turn_to_df))
names_res = names(res)

Rcpp::cppFunction('
	void load_all_data_cpp(SEXP res, CharacterVector fields, CharacterVector all_data_files, Rcpp::Function func) {
	  Function R_load("load");
	  Environment env = Environment::global_env();
		for (int i = 0; i < all_data_files.length(); i++){
  		if (i % 1000 == 0){
  		  Rcout << "i: " << i << " filename: " << all_data_files[i] << std::endl;
  		}
		  String s(all_data_files[i]);
		  R_load(s);
		  List res_row = env["res_row"];
		  res_row("filename") = s;
		  func(res, i, fields, res_row);
		}
	}
', depends = "data.table")

load_all_data_cpp(res, names_res, all_data_files, data.table::set)
save(res, file = "all_sim_data.RData")


load("all_sim_data.RData")
res[, inference_method := sub("result_\\d*_(.*$)", "\\1", filename, perl = TRUE)]
res[design == "KK14" & test_type == "randomization-exact" & inference_method == "survival_multivariate_weibull_regression" & nsim == 1][, .(inference_type = sub("result_\\d*_(.*$)", "\\1", filename,perl=TRUE))]
res[, design := factor(design, levels = c("CRD", "iBCRD", "Efron", "Atkinson", "KK14", "KK21", "KK21stepwise"))]
res = res[!is.na(res$p_val), ]
res_summ_pval = res[, .(avg_p_val = mean(p_val < 0.05), num_sim = .N), by = c("n", "design", "test_type", "inference_method", "prob_of_adding_response", "betaT")]
res_summ_pval[, pval_test := 2 * pnorm(-abs(avg_p_val - 0.05) / sqrt(0.05 * .95 / num_sim))]
# res_summ_pval[, pval_bad := ifelse(pval_test < 0.05 / .N, "*", "")]
res_summ_pval[pval_test < 0.05 / .N & !grepl("simple", inference_method) & !grepl("matching_dummies", inference_method) & betaT == 0]

table(res$estimand)
res_summ_ci = res[test_type == "MLE-or-KM-based", .(avg_coverage = mean(ci_a <= estimand & estimand <= ci_b), num_sim = .N), by = c("n", "design", "inference_method", "prob_of_adding_response", "betaT")]
res_summ_ci[, pval_test := 2 * pnorm(-abs(avg_coverage - 0.95) / sqrt(0.95 * .05 / num_sim))]
# res_summ_ci[, pval_bad := ifelse(pval_test < 0.05 / .N, "*", "")]
res_summ_ci[pval_test < 0.05 / .N]

res_summ_est = res[, .(med_sq_err = round(10000 * median((estimate - estimand)^2)), num_sim = .N), by = c("n", "response_type", "design", "test_type", "inference_method", "prob_of_adding_response", "betaT")]
res_summ_est = res_summ_est[!(!grepl("KK21", design) & prob_of_adding_response != 0.5), ]
data.frame(res_summ_est[grepl("simple", inference_method) | grepl("KK_compound_mean", inference_method)][betaT == 1][order(response_type, test_type, inference_method, design)])
data.frame(res_summ_est[grepl("covar", inference_method)][order(response_type, test_type, inference_method, design)])
data.frame(res_summ_est[grepl("univ", inference_method)][order(response_type, test_type, inference_method, design)])
data.frame(res_summ_est[grepl("multi", inference_method)][order(response_type, test_type, inference_method, design)])


