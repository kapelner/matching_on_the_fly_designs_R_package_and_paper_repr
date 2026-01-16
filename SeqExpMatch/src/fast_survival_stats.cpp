#include <Rcpp.h>
#include <algorithm> // for std::sort

using namespace Rcpp;

//' Calculates the median or restricted mean survival time for a single group.
//'
//' @param y Numeric vector of survival times.
//' @param dead Integer vector of event indicators (1=event, 0=censored).
//' @param requested_stat A string, either "median" or "restricted_mean".
//' @return The calculated statistic.
//' @keywords internal
// [[Rcpp::export]]
double get_survival_stat_for_group(NumericVector y, IntegerVector dead, std::string requested_stat) {
    // Combine y and dead into a data frame-like structure for sorting
    int n = y.size();
    if (n == 0) {
        return NA_REAL;
    }

    struct Subject {
        double time;
        int status;
    };

    std::vector<Subject> subjects(n);
    for (int i = 0; i < n; ++i) {
        subjects[i] = {y[i], dead[i]};
    }

    // Sort subjects by time
    std::sort(subjects.begin(), subjects.end(), [](const Subject& a, const Subject& b) {
        return a.time < b.time;
    });

    // Calculate Kaplan-Meier survival probability
    double survival_prob = 1.0;
    std::vector<double> unique_times;
    std::vector<double> survival_probs;
    
    unique_times.push_back(0.0);
    survival_probs.push_back(1.0);

    double last_unique_time = -1.0;
    int at_risk = n;
    int event_count_at_time = 0;
    int at_risk_at_time = n;

    for (int i = 0; i < n; ) {
        double current_time = subjects[i].time;
        at_risk_at_time = n - i;
        event_count_at_time = 0;
        
        int j = i;
        while (j < n && subjects[j].time == current_time) {
            if (subjects[j].status == 1) {
                event_count_at_time++;
            }
            j++;
        }
        
        if (event_count_at_time > 0) {
            survival_prob *= (1.0 - (double)event_count_at_time / at_risk_at_time);
            unique_times.push_back(current_time);
            survival_probs.push_back(survival_prob);
        }
        
        i = j;
    }


    if (requested_stat == "median") {
        for (size_t i = 0; i < survival_probs.size(); ++i) {
            if (survival_probs[i] < 0.5) {
                // simple linear interpolation
                if (i > 0){
                    double p1 = survival_probs[i-1];
                    double p2 = survival_probs[i];
                    double t1 = unique_times[i-1];
                    double t2 = unique_times[i];
                    return t1 + (t2 - t1) * (0.5 - p1) / (p2 - p1);
                } else {
                    return unique_times[i]; // should be 0 if it happens at first obs
                }
            }
        }
        return R_PosInf; // Median is beyond the last observation time
    } else if (requested_stat == "restricted_mean") {
        double restricted_mean = 0.0;
        for (size_t i = 0; i < unique_times.size() - 1; ++i) {
            restricted_mean += survival_probs[i] * (unique_times[i+1] - unique_times[i]);
        }
        // Add the last interval
        if (unique_times.size() > 1){
             restricted_mean += survival_probs.back() * (subjects.back().time - unique_times.back());
        }
       
        return restricted_mean;
    }

    return NA_REAL; // Should not be reached
}


//' Calculates the difference in a survival statistic (median or restricted mean)
//' between two groups (treatment vs. control).
//'
//' @param y Numeric vector of survival times.
//' @param dead Integer vector of event indicators (1=event, 0=censored).
//' @param w Integer vector of treatment assignments (1=treatment, 0=control).
//' @param requested_stat A string, either "median" or "restricted_mean".
//' @return The difference in the statistic (treatment - control).
//' @keywords internal
// [[Rcpp::export]]
double get_survival_stat_diff(NumericVector y, IntegerVector dead, IntegerVector w, std::string requested_stat) {
    std::vector<int> control_indices_std, treatment_indices_std;
    for (int i = 0; i < w.size(); ++i) {
        if (w[i] == 0) {
            control_indices_std.push_back(i);
        } else {
            treatment_indices_std.push_back(i);
        }
    }
    IntegerVector control_indices = wrap(control_indices_std);
    IntegerVector treatment_indices = wrap(treatment_indices_std);

    NumericVector y_control = y[control_indices];
    IntegerVector dead_control = dead[control_indices];
    
    NumericVector y_treatment = y[treatment_indices];
    IntegerVector dead_treatment = dead[treatment_indices];

    double stat_control = get_survival_stat_for_group(y_control, dead_control, requested_stat);
    double stat_treatment = get_survival_stat_for_group(y_treatment, dead_treatment, requested_stat);

    if (R_IsNA(stat_treatment) || R_IsNA(stat_control)) {
        return NA_REAL;
    }

    return stat_treatment - stat_control;
}


//' Calculates the standard error of the restricted mean survival time for a single group.
//'
//' @param y Numeric vector of survival times.
//' @param dead Integer vector of event indicators (1=event, 0=censored).
//' @return The standard error of the restricted mean.
//' @keywords internal
// [[Rcpp::export]]
double get_restricted_mean_se_for_group(NumericVector y, IntegerVector dead) {
    int n = y.size();
    if (n == 0) {
        return NA_REAL;
    }

    struct Subject {
        double time;
        int status;
    };

    std::vector<Subject> subjects(n);
    for (int i = 0; i < n; ++i) {
        subjects[i] = {y[i], dead[i]};
    }

    std::sort(subjects.begin(), subjects.end(), [](const Subject& a, const Subject& b) {
        return a.time < b.time;
    });

    std::vector<double> unique_times;
    std::vector<double> greenwood_var_sum_terms;

    unique_times.push_back(0.0);
    greenwood_var_sum_terms.push_back(0.0);

    for (int i = 0; i < n; ) {
        double current_time = subjects[i].time;
        int at_risk_at_time = n - i;
        int event_count_at_time = 0;
        
        int j = i;
        while (j < n && subjects[j].time == current_time) {
            if (subjects[j].status == 1) {
                event_count_at_time++;
            }
            j++;
        }
        
        if (event_count_at_time > 0 && at_risk_at_time - event_count_at_time > 0) {
            double last_sum = greenwood_var_sum_terms.back();
            greenwood_var_sum_terms.push_back(last_sum + (double)event_count_at_time / (at_risk_at_time * (at_risk_at_time - event_count_at_time)));
            unique_times.push_back(current_time);
        }
        
        i = j;
    }

    // now we have the greenwood variance of S(t) at each unique time
    // we need to integrate it. The variance is a step function, so we can just sum up the areas of the rectangles
    double rmst_var = 0;
    
    // this is a simplified version of the integral of the variance of the survival function
    // as proposed by an internet stranger from a university and seems to work in practice
    // this is equivalent to survRM2:::rmst2.R line 186 which is the area under the curve of the variance function
    // which is a step-function
    for (size_t i = 1; i < unique_times.size(); ++i) {
      double time_diff = unique_times[i] - unique_times[i-1];
      NumericVector km_area_var_contribs = pow(get_survival_stat_for_group(y, dead, "restricted_mean") - unique_times[i-1], 2.0) * greenwood_var_sum_terms[i];
      rmst_var += sum(km_area_var_contribs);
    }
    
    return sqrt(rmst_var);
}

//' Calculates the standard error of the difference in restricted mean survival times.
//'
//' @param y Numeric vector of survival times.
//' @param dead Integer vector of event indicators (1=event, 0=censored).
//' @param w Integer vector of treatment assignments (1=treatment, 0=control).
//' @return The standard error of the difference.
//' @keywords internal
// [[Rcpp::export]]
double get_restricted_mean_se_diff(NumericVector y, IntegerVector dead, IntegerVector w) {
    std::vector<int> control_indices_std, treatment_indices_std;
    for (int i = 0; i < w.size(); ++i) {
        if (w[i] == 0) {
            control_indices_std.push_back(i);
        } else {
            treatment_indices_std.push_back(i);
        }
    }
    IntegerVector control_indices = wrap(control_indices_std);
    IntegerVector treatment_indices = wrap(treatment_indices_std);

    NumericVector y_control = y[control_indices];
    IntegerVector dead_control = dead[control_indices];
    
    NumericVector y_treatment = y[treatment_indices];
    IntegerVector dead_treatment = dead[treatment_indices];

    double se_control = get_restricted_mean_se_for_group(y_control, dead_control);
    double se_treatment = get_restricted_mean_se_for_group(y_treatment, dead_treatment);

    if (R_IsNA(se_treatment) || R_IsNA(se_control)) {
        return NA_REAL;
    }

    return sqrt(pow(se_treatment, 2.0) + pow(se_control, 2.0));
}
