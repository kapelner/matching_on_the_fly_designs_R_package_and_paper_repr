pacman::p_load(stringr, data.table, Rcpp)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/res"))

all_data_files = str_subset(dir(), regex("result_"))
load(all_data_files[1])
res = data.frame(matrix(NA_real_, nrow = max(as.numeric(str_match(all_data_files, regex("_(\\d*)"))[,2])), ncol = length(res_row)))
names_res_row = names(res_row)
res_row$weights = NA_character_
list_to_turn_to_df = lapply(res_row, function(x){rep(x, length(all_data_files))})
res = data.table(data.frame(list_to_turn_to_df))

Rcpp::cppFunction('
	void load_all_data_cpp(SEXP res, CharacterVector fields, CharacterVector all_data_files, Rcpp::Function func) {
	  Function R_load("load"); 
	  //Function DT_set("data.table::set");
	  Environment env = Environment::global_env();
		for (int i = 0; i < all_data_files.length(); i++){
  		if (i % 1000 == 0){
  		  Rcout << "i: " << i << " filename: " << all_data_files[i] << std::endl;
  		}
		  String s(all_data_files[i]);
		  R_load(s);
		  List res_row = env["res_row"];
		  func(res, i, fields, res_row);
		}
	}
', depends="data.table")

load_all_data_cpp(res, names_res_row, all_data_files, data.table::set)


res_df[, .(avg_pval = mean(p_val)), by = n]
