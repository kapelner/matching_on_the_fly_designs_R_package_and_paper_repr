#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector efron_redraw_cpp(int t, double prob_T, double weighted_coin_prob) {
	std::vector<double> results_vec(t);
    double* w_ptr = results_vec.data();

	// RNG scope for proper random number generation
	RNGScope scope;

	// Assign treatments sequentially using Efron's biased coin
	int n_T = 0;
	int n_C = 0;

	for (int i = 0; i < t; ++i) {
	double assignment_prob;

	// Calculate imbalance scores
	double score_T = n_T * prob_T;
	double score_C = n_C * (1.0 - prob_T);

	if (score_T > score_C) {
		// More treatments than expected, bias towards control
		assignment_prob = 1.0 - weighted_coin_prob;
	} else if (score_T < score_C) {
		// More controls than expected, bias towards treatment
		assignment_prob = weighted_coin_prob;
	} else {
		// Balanced, use prob_T
		assignment_prob = prob_T;
	}

	// Flip biased coin
	double u = R::unif_rand();
	if (u < assignment_prob) {
		w_ptr[i] = 1.0;
		n_T++;
	} else {
		w_ptr[i] = 0.0;
		n_C++;
	}
	}

	return wrap(results_vec);
}
