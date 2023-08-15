#pragma once

#include<vector>
using namespace std;

class ThomasSolver
{
public:
	ThomasSolver(const vector<double>& lower_diag, const vector<double>& center_diag, const vector<double>& upper_diag, const vector<double>& rhs);

	// returns [X_1, X_2, ...]
	vector<double> solve() const;

private:
	vector<double> _lower_diagonal;   // a_j
	vector<double> _center_diagonal;  // b_j
	vector<double> _upper_diagonal;   // c_j
	vector<double> _right_hand_side;  // R_j
};
