#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector get_column_types_cpp(List df) {
  int n_cols = df.size();
  CharacterVector types(n_cols);

  for (int j = 0; j < n_cols; ++j) {
    SEXP col = df[j];

    // Get the class attribute
    SEXP class_attr = Rf_getAttrib(col, R_ClassSymbol);

    if (class_attr != R_NilValue && Rf_length(class_attr) > 0) {
      // Return the first class as a string
      CharacterVector col_class(class_attr);
      types[j] = col_class[0];
    } else {
      // If no class attribute, determine from SEXP type
      int sexp_type = TYPEOF(col);
      switch(sexp_type) {
        case INTSXP:
          types[j] = "integer";
          break;
        case REALSXP:
          types[j] = "numeric";
          break;
        case STRSXP:
          types[j] = "character";
          break;
        case LGLSXP:
          types[j] = "logical";
          break;
        default:
          types[j] = "unknown";
      }
    }
  }

  // Set names on the result vector
  SEXP names_attr = Rf_getAttrib(df, R_NamesSymbol);
  if (names_attr != R_NilValue) {
    types.names() = CharacterVector(names_attr);
  }

  return types;
}
