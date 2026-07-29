#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace Rcpp;

namespace {

inline uint64_t bounded_rand(std::mt19937_64& rng, uint64_t s) {
	uint64_t x = rng();
	__uint128_t m = static_cast<__uint128_t>(x) * s;
	uint64_t l = static_cast<uint64_t>(m);
	if (l < s) {
		uint64_t t = (-s) % s;
		while (l < t) {
			x = rng();
			m = static_cast<__uint128_t>(x) * s;
			l = static_cast<uint64_t>(m);
		}
	}
	return static_cast<uint64_t>(m >> 64);
}

inline std::mt19937_64 make_local_rng() {
	return std::mt19937_64(static_cast<uint64_t>(
		R::unif_rand() * static_cast<double>(std::numeric_limits<uint64_t>::max())));
}

std::string char_at_or_na(const CharacterVector& x, int i) {
	if (i < 0 || i >= x.size() || STRING_ELT(x, i) == NA_STRING) return "NA";
	return Rcpp::as<std::string>(x[i]);
}

std::vector<int> proportional_allocation(int total_size, const std::vector<std::vector<int> >& groups, bool replace) {
	const int G = static_cast<int>(groups.size());
	std::vector<int> alloc(G, 0);
	std::vector<double> raw(G, 0.0);
	std::vector<double> rem(G, 0.0);
	int n_total = 0;
	for (int g = 0; g < G; ++g) n_total += static_cast<int>(groups[g].size());
	if (n_total <= 0) stop("No exchangeable units are available.");

	int assigned = 0;
	for (int g = 0; g < G; ++g) {
		raw[g] = static_cast<double>(total_size) * static_cast<double>(groups[g].size()) / static_cast<double>(n_total);
		alloc[g] = static_cast<int>(std::floor(raw[g]));
		rem[g] = raw[g] - static_cast<double>(alloc[g]);
		assigned += alloc[g];
	}

	int remaining = total_size - assigned;
	std::vector<int> ord(G);
	std::iota(ord.begin(), ord.end(), 0);
	std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return rem[a] > rem[b]; });
	for (int i = 0; i < remaining; ++i) ++alloc[ord[i]];

	if (!replace) {
		for (int g = 0; g < G; ++g) {
			if (alloc[g] > static_cast<int>(groups[g].size())) {
				stop("Subsampling stratum is too small for the requested subsample size.");
			}
		}
	}
	return alloc;
}

void sample_group_units(
	const std::vector<int>& ids,
	int k,
	bool replace,
	std::mt19937_64& rng,
	std::vector<int>& selected
) {
	if (k <= 0) return;
	const int n = static_cast<int>(ids.size());
	if (n <= 0) stop("No exchangeable units are available.");
	if (!replace && k > n) stop("Subsampling group is too small for the requested subsample size.");

	if (replace) {
		const uint64_t un = static_cast<uint64_t>(n);
		for (int i = 0; i < k; ++i) {
			selected.push_back(ids[static_cast<int>(bounded_rand(rng, un))]);
		}
		return;
	}

	std::unordered_set<int> seen;
	seen.reserve(static_cast<size_t>(k) * 2U + 1U);
	const uint64_t un = static_cast<uint64_t>(n);
	const int target_size = static_cast<int>(selected.size()) + k;
	while (static_cast<int>(selected.size()) < target_size) {
		int pos = static_cast<int>(bounded_rand(rng, un));
		if (seen.insert(pos).second) selected.push_back(ids[pos]);
	}
}

IntegerVector renumber_positive_match_ids(const std::vector<int>& raw_m) {
	std::vector<int> positives;
	positives.reserve(raw_m.size());
	for (int v : raw_m) {
		if (v > 0) positives.push_back(v);
	}
	std::sort(positives.begin(), positives.end());
	positives.erase(std::unique(positives.begin(), positives.end()), positives.end());

	std::unordered_map<int, int> map;
	map.reserve(positives.size() * 2U + 1U);
	for (int i = 0; i < static_cast<int>(positives.size()); ++i) {
		map.emplace(positives[i], i + 1);
	}

	IntegerVector out(raw_m.size());
	for (int i = 0; i < static_cast<int>(raw_m.size()); ++i) {
		int v = raw_m[i];
		if (v > 0) out[i] = map[v];
		else out[i] = 0;
	}
	return out;
}

List build_draw(
	const std::vector<std::vector<int> >& units,
	const std::vector<int>& unit_lengths,
	const std::vector<int>& selected_unit_ids,
	const std::vector<std::string>& unit_kind,
	const IntegerVector& m_vec_full,
	bool has_unit_kind,
	bool has_m_vec_full,
	bool preserve_order,
	bool identity_singletons,
	const std::string& unit_type,
	const std::string& size_label,
	int n_units
) {
	if (identity_singletons && !has_unit_kind && !has_m_vec_full) {
		const int n_selected = static_cast<int>(selected_unit_ids.size());
		IntegerVector i_b(n_selected);
		IntegerVector unit_ids(n_selected);
		for (int i = 0; i < n_selected; ++i) {
			const int row = selected_unit_ids[i] + 1;
			i_b[i] = row;
			unit_ids[i] = row;
		}
		if (preserve_order && n_selected > 1) {
			std::sort(i_b.begin(), i_b.end());
		}
		List draw = List::create(
			Named("i_b") = i_b,
			Named("m_vec_b") = R_NilValue,
			Named("unit_ids") = unit_ids,
			Named("unit_type") = unit_type,
			Named("n_units") = n_units
		);
		draw.push_back(n_selected, size_label);
		return draw;
	}

	int total_rows = 0;
	for (int unit_id : selected_unit_ids) total_rows += unit_lengths[unit_id];

	std::vector<int> rows;
	rows.reserve(total_rows);
	std::vector<int> raw_m;
	bool has_m = has_unit_kind || has_m_vec_full;
	if (has_m) raw_m.reserve(total_rows);

	int pair_id = 0;
	for (int unit_id : selected_unit_ids) {
		const std::vector<int>& unit_rows = units[unit_id];
		if (has_unit_kind) {
			const bool is_pair = unit_kind[unit_id] == "pair";
			if (is_pair) ++pair_id;
			const int out_m = is_pair ? pair_id : 0;
			for (int row : unit_rows) {
				rows.push_back(row);
				raw_m.push_back(out_m);
			}
		} else if (has_m_vec_full) {
			for (int row : unit_rows) {
				rows.push_back(row);
				int m_val = 0;
				const int pos = row - 1;
				if (pos >= 0 && pos < m_vec_full.size() && m_vec_full[pos] != NA_INTEGER) {
					m_val = m_vec_full[pos];
				}
				raw_m.push_back(m_val > 0 ? m_val : 0);
			}
		} else {
			rows.insert(rows.end(), unit_rows.begin(), unit_rows.end());
		}
	}

	if (preserve_order && rows.size() > 1U) {
		std::vector<int> ord(rows.size());
		std::iota(ord.begin(), ord.end(), 0);
		std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return rows[a] < rows[b]; });
		std::vector<int> sorted_rows(rows.size());
		std::vector<int> sorted_m;
		if (has_m) sorted_m.resize(raw_m.size());
		for (int i = 0; i < static_cast<int>(ord.size()); ++i) {
			sorted_rows[i] = rows[ord[i]];
			if (has_m) sorted_m[i] = raw_m[ord[i]];
		}
		rows.swap(sorted_rows);
		if (has_m) raw_m.swap(sorted_m);
	}

	IntegerVector i_b(rows.size());
	std::copy(rows.begin(), rows.end(), i_b.begin());

	IntegerVector unit_ids(selected_unit_ids.size());
	for (int i = 0; i < static_cast<int>(selected_unit_ids.size()); ++i) {
		unit_ids[i] = selected_unit_ids[i] + 1;
	}

	RObject m_vec_b = R_NilValue;
	if (has_m) {
		m_vec_b = renumber_positive_match_ids(raw_m);
	}

	List draw = List::create(
		Named("i_b") = i_b,
		Named("m_vec_b") = m_vec_b,
		Named("unit_ids") = unit_ids,
		Named("unit_type") = unit_type,
		Named("n_units") = n_units
	);
	draw.push_back(static_cast<int>(selected_unit_ids.size()), size_label);
	return draw;
}

} // namespace

// [[Rcpp::export]]
List exchangeable_resampling_draws_cpp(
	List units,
	SEXP strata_ids_sexp,
	SEXP unit_kind_sexp,
	SEXP m_vec_full_sexp,
	int B,
	int size,
	bool replace,
	bool stratified,
	bool preserve_order,
	std::string unit_type,
	std::string size_label
) {
	if (B < 0) stop("B must be non-negative.");
	if (size < 0) stop("Resampling size must be non-negative.");
	const int n_units = units.size();
	if (n_units <= 0) stop("No exchangeable units are available.");

	std::vector<std::vector<int> > unit_rows(n_units);
	std::vector<int> unit_lengths(n_units);
	bool identity_singletons = true;
	for (int i = 0; i < n_units; ++i) {
		IntegerVector u = units[i];
		unit_rows[i].assign(u.begin(), u.end());
		unit_lengths[i] = static_cast<int>(unit_rows[i].size());
		if (unit_lengths[i] != 1 || unit_rows[i][0] != i + 1) {
			identity_singletons = false;
		}
	}

	const bool has_unit_kind = !Rf_isNull(unit_kind_sexp);
	std::vector<std::string> unit_kind;
	if (has_unit_kind) {
		CharacterVector uk(unit_kind_sexp);
		if (uk.size() != n_units) stop("unit_kind must have one entry per exchangeable unit.");
		unit_kind.reserve(n_units);
		for (int i = 0; i < n_units; ++i) unit_kind.push_back(char_at_or_na(uk, i));
	}

	const bool has_m_vec_full = !Rf_isNull(m_vec_full_sexp);
	IntegerVector m_vec_full;
	if (has_m_vec_full) m_vec_full = IntegerVector(m_vec_full_sexp);

	std::vector<std::vector<int> > groups;
	if (stratified && !Rf_isNull(strata_ids_sexp)) {
		CharacterVector strata(strata_ids_sexp);
		if (strata.size() != n_units) stop("strata_ids must have one entry per exchangeable unit.");
		std::unordered_map<std::string, int> group_index;
		group_index.reserve(static_cast<size_t>(n_units) * 2U + 1U);
		for (int i = 0; i < n_units; ++i) {
			std::string key = char_at_or_na(strata, i);
			auto it = group_index.find(key);
			if (it == group_index.end()) {
				const int g = static_cast<int>(groups.size());
				group_index.emplace(key, g);
				groups.push_back(std::vector<int>());
				groups[g].push_back(i);
			} else {
				groups[it->second].push_back(i);
			}
		}
	} else {
		groups.resize(1);
		groups[0].resize(n_units);
		std::iota(groups[0].begin(), groups[0].end(), 0);
	}

	std::vector<int> alloc = proportional_allocation(size, groups, replace);
	List out(B);
	auto rng = make_local_rng();
	for (int b = 0; b < B; ++b) {
		std::vector<int> selected;
		selected.reserve(size);
		for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
			sample_group_units(groups[g], alloc[g], replace, rng, selected);
		}
		out[b] = build_draw(
			unit_rows,
			unit_lengths,
			selected,
			unit_kind,
			m_vec_full,
			has_unit_kind,
			has_m_vec_full && !has_unit_kind,
			preserve_order,
			identity_singletons,
			unit_type,
			size_label,
			n_units
		);
	}
	return out;
}
