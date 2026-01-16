#include <Rcpp.h>
#include <unordered_set>

using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector columns_have_missingness_cpp(List df) {
  int n_cols = df.size();
  LogicalVector has_missing(n_cols);

  for (int j = 0; j < n_cols; ++j) {
    SEXP col = df[j];
    int col_type = TYPEOF(col);
    bool found_na = false;

    switch(col_type) {
      case REALSXP: {
        NumericVector vec(col);
        for (int i = 0; i < vec.size(); ++i) {
          if (NumericVector::is_na(vec[i])) {
            found_na = true;
            break;
          }
        }
        break;
      }
      case INTSXP: {
        IntegerVector vec(col);
        for (int i = 0; i < vec.size(); ++i) {
          if (IntegerVector::is_na(vec[i])) {
            found_na = true;
            break;
          }
        }
        break;
      }
      case STRSXP: {
        CharacterVector vec(col);
        for (int i = 0; i < vec.size(); ++i) {
          if (CharacterVector::is_na(vec[i])) {
            found_na = true;
            break;
          }
        }
        break;
      }
      case LGLSXP: {
        LogicalVector vec(col);
        for (int i = 0; i < vec.size(); ++i) {
          if (LogicalVector::is_na(vec[i])) {
            found_na = true;
            break;
          }
        }
        break;
      }
      default:
        found_na = false;
    }

    has_missing[j] = found_na;
  }

  // Set names if the input has names
  SEXP names_attr = Rf_getAttrib(df, R_NamesSymbol);
  if (names_attr != R_NilValue) {
    has_missing.names() = CharacterVector(names_attr);
  }

  return has_missing;
}

// [[Rcpp::export]]
List create_missingness_indicators_cpp(List df, IntegerVector col_indices) {
  int n_cols = col_indices.size();
  int n_rows = 0;

  if (df.size() > 0) {
    SEXP first_col = df[0];
    n_rows = Rf_length(first_col);
  }

  List result(n_cols);
  CharacterVector result_names(n_cols);

  for (int j = 0; j < n_cols; ++j) {
    int col_idx = col_indices[j] - 1; // R to C++ indexing
    SEXP col = df[col_idx];
    int col_type = TYPEOF(col);

    IntegerVector indicator(n_rows);

    switch(col_type) {
      case REALSXP: {
        NumericVector vec(col);
        for (int i = 0; i < n_rows; ++i) {
          indicator[i] = NumericVector::is_na(vec[i]) ? 1 : 0;
        }
        break;
      }
      case INTSXP: {
        IntegerVector vec(col);
        for (int i = 0; i < n_rows; ++i) {
          indicator[i] = IntegerVector::is_na(vec[i]) ? 1 : 0;
        }
        break;
      }
      case STRSXP: {
        CharacterVector vec(col);
        for (int i = 0; i < n_rows; ++i) {
          indicator[i] = CharacterVector::is_na(vec[i]) ? 1 : 0;
        }
        break;
      }
      case LGLSXP: {
        LogicalVector vec(col);
        for (int i = 0; i < n_rows; ++i) {
          indicator[i] = LogicalVector::is_na(vec[i]) ? 1 : 0;
        }
        break;
      }
      default:
        // Unknown type, mark all as 0
        std::fill(indicator.begin(), indicator.end(), 0);
    }

    result[j] = indicator;

    // Generate name for this indicator column
    CharacterVector col_names = df.names();
    if (col_names.size() > col_idx) {
      std::string base_name = Rcpp::as<std::string>(col_names[col_idx]);
      result_names[j] = base_name + "_is_missing";
    } else {
      result_names[j] = "V" + std::to_string(col_idx + 1) + "_is_missing";
    }
  }

  result.names() = result_names;
  return result;
}

// [[Rcpp::export]]
IntegerVector count_unique_values_cpp(List df) {
  int n_cols = df.size();
  IntegerVector unique_counts(n_cols);

  for (int j = 0; j < n_cols; ++j) {
    SEXP col = df[j];
    int col_type = TYPEOF(col);
    int unique_count = 0;

    switch(col_type) {
      case REALSXP: {
        NumericVector vec(col);
        std::unordered_set<double> unique_vals;
        for (int i = 0; i < vec.size(); ++i) {
          if (!NumericVector::is_na(vec[i])) {
            unique_vals.insert(vec[i]);
          }
        }
        unique_count = unique_vals.size();
        // Add 1 if there are any NAs (NA counts as a unique value)
        for (int i = 0; i < vec.size(); ++i) {
          if (NumericVector::is_na(vec[i])) {
            unique_count++;
            break;
          }
        }
        break;
      }
      case INTSXP: {
        IntegerVector vec(col);
        std::unordered_set<int> unique_vals;
        bool has_na = false;
        for (int i = 0; i < vec.size(); ++i) {
          if (IntegerVector::is_na(vec[i])) {
            has_na = true;
          } else {
            unique_vals.insert(vec[i]);
          }
        }
        unique_count = unique_vals.size() + (has_na ? 1 : 0);
        break;
      }
      case STRSXP: {
        CharacterVector vec(col);
        std::unordered_set<std::string> unique_vals;
        bool has_na = false;
        for (int i = 0; i < vec.size(); ++i) {
          if (CharacterVector::is_na(vec[i])) {
            has_na = true;
          } else {
            unique_vals.insert(Rcpp::as<std::string>(vec[i]));
          }
        }
        unique_count = unique_vals.size() + (has_na ? 1 : 0);
        break;
      }
      case LGLSXP: {
        LogicalVector vec(col);
        std::unordered_set<int> unique_vals;
        bool has_na = false;
        for (int i = 0; i < vec.size(); ++i) {
          if (LogicalVector::is_na(vec[i])) {
            has_na = true;
          } else {
            unique_vals.insert(vec[i]);
          }
        }
        unique_count = unique_vals.size() + (has_na ? 1 : 0);
        break;
      }
      default:
        unique_count = 0;
    }

    unique_counts[j] = unique_count;
  }

  // Set names if the input has names
  SEXP names_attr = Rf_getAttrib(df, R_NamesSymbol);
  if (names_attr != R_NilValue) {
    unique_counts.names() = CharacterVector(names_attr);
  }

  return unique_counts;
}
