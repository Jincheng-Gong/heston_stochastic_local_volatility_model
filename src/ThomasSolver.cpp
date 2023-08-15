#include "ThomasSolver.h"

ThomasSolver::ThomasSolver(const vector<double>& lower_diag, const vector<double>& center_diag, const vector<double>& upper_diag, const vector<double>& rhs)
	: _lower_diagonal(lower_diag), _center_diagonal(center_diag), _upper_diagonal(upper_diag), _right_hand_side(rhs)
{
}

vector<double> ThomasSolver::solve() const
{
	size_t rhs_size = _right_hand_side.size();

	// 1) Compute new coeffs [sweep step]
	vector<double> upper_diag_prime; // c'_j;
	upper_diag_prime.push_back(_upper_diagonal[0] / _center_diagonal[0]);
	for (size_t j = 1; j < rhs_size - 1; ++j)
		upper_diag_prime.push_back(_upper_diagonal[j] / (_center_diagonal[j] - _lower_diagonal[j - 1] * upper_diag_prime[j - 1]));

	vector<double> rhs_prime;
	rhs_prime.push_back(_right_hand_side[0] / _center_diagonal[0]);
	for (size_t j = 1; j < rhs_size; ++j)
		rhs_prime.push_back((_right_hand_side[j] - _lower_diagonal[j - 1] * rhs_prime[j - 1]) / (_center_diagonal[j] - _lower_diagonal[j - 1] * upper_diag_prime[j - 1]));

	// 2) Back substitution step:
	vector<double> solution(rhs_size);
	solution[rhs_size - 1] = rhs_prime[rhs_size - 1];
	// for (size_t j = rhs_size - 2; j >= 0; --j) Does not work -> infinite loop -> -1 = + infty [max number of size_t]
	for (size_t j = rhs_size - 2; j < rhs_size - 1; --j)
		solution[j] = rhs_prime[j] - upper_diag_prime[j] * solution[j + 1];

	return solution;

}