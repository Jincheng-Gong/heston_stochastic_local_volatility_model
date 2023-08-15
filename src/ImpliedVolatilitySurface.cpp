#include "ImpliedVolatilitySurface.h"
#include<cmath>
#include <algorithm>
#include "ThomasSolver.h"
using namespace std;


// Initialization List !!!
ImpliedVolatilitySurface::ImpliedVolatilitySurface(const vector<double> & maturities, const vector<double> & strikes, const vector< vector<double> > & implied_vols, const double & risk_free_rate)
	: _maturities(maturities), _strikes(strikes), _market_implied_vols(implied_vols), _risk_free_rate(risk_free_rate)
{
	if (!check_ordered_vectors())
		throw "The maturities and strikes must be ordered and positive [strictly positive for maturities]";
	if (!check_dimensions())
		throw "The dimensions of the implied vol matrix must be (maturity dim x strike dim)";
	if (!check_vols_positivity())
		throw "The implied volatilities must be greater than 0.";

	// _delta_K[j] = K_{j+1} - K_j
	for (size_t j = 0; j < _strikes.size() - 1; ++j)
		_delta_K.push_back(_strikes[j + 1] - _strikes[j]);

	evaluate_cubic_spline_coefficients();

}

bool ImpliedVolatilitySurface::check_ordered_vectors() const
{
	return is_sorted(_maturities.begin(), _maturities.end()) && is_sorted(_strikes.begin(), _strikes.end()) && _maturities[0] > 0 && _strikes[0] >= 0;
}

bool ImpliedVolatilitySurface::check_dimensions() const
{
	return _market_implied_vols.size() == _maturities.size() && _market_implied_vols[0].size() == _strikes.size();
}

bool ImpliedVolatilitySurface::check_vols_positivity() const
{
	size_t maturity_dim = _market_implied_vols.size();
	size_t strike_dim = _market_implied_vols[0].size();

	for (size_t mat_idx = 0; mat_idx < maturity_dim; ++mat_idx)
		for (size_t strike_idx = 0; strike_idx < strike_dim; ++strike_idx)
			if (_market_implied_vols[mat_idx][strike_idx] <= 0.)
				return false;

	// If the loop has finished, it means that no negative value has been found
	return true;
}

double ImpliedVolatilitySurface::implied_volatility(const double & maturity, const double & strike) const
{
	// Locate the maturity in the correct maturity interval from _maturities : T in [T_i, T_{i+1}[ or T <= T_1 or T >= T_M
	size_t left_maturity_index;
	size_t right_maturity_index;

	// left extrapolation case : T <= T_1
	if (maturity <= _maturities[0])
	{
		left_maturity_index = 0;
		double left_strike = strike * std::exp(_risk_free_rate * (_maturities[left_maturity_index] - maturity));
		double left_extrapolated_implied_vol = compute_smile_implied_vol(left_maturity_index, left_strike);
		return left_extrapolated_implied_vol;
	}
	else
	{
		// right extrapolation case : T >= T_M
		if (maturity >= _maturities[maturity_size() - 1])
		{
			right_maturity_index = maturity_size() - 1;
			double right_strike = strike * std::exp(_risk_free_rate * (_maturities[right_maturity_index] - maturity));
			double right_extrapolated_implied_vol = compute_smile_implied_vol(right_maturity_index, right_strike);
			return right_extrapolated_implied_vol;
		}
		else
			// interpolation case : T \in [T_i, T_{i+1}[	
		{
			right_maturity_index = 0;
			// the loop stops when T < T_{i+1} so that T_i <= T < T_{i+1} 
			while (maturity >= _maturities[right_maturity_index])
				right_maturity_index++;
			left_maturity_index = right_maturity_index - 1;

			// Evaluate K^(i) and K^(i+1) - No need for S_0, just rate is needed, hence passed as a private data of the object.
			double left_strike = strike * exp(_risk_free_rate * (_maturities[left_maturity_index] - maturity));
			double right_strike = strike * exp(_risk_free_rate * (_maturities[right_maturity_index] - maturity));

			// Compute v(T_i, k_F) and v(T_{i+1}, k_F) : use of cubic spline interpolation/extrapolation
			double left_implied_vol = compute_smile_implied_vol(left_maturity_index, left_strike);
			double left_total_variance = left_implied_vol * left_implied_vol * _maturities[left_maturity_index];

			double right_implied_vol = compute_smile_implied_vol(right_maturity_index, right_strike);
			double right_total_variance = right_implied_vol * right_implied_vol * _maturities[right_maturity_index];

			// Linear interpolation in variance along maturities
			double implied_total_variance = left_total_variance + ((right_total_variance - left_total_variance) / (_maturities[right_maturity_index] - _maturities[left_maturity_index])) * (maturity - _maturities[left_maturity_index]);

			double impliedVol = sqrt(implied_total_variance / maturity);
			return impliedVol;
		}
	}
}

void ImpliedVolatilitySurface::evaluate_cubic_spline_coefficients()
{
	for (size_t mat_idx = 0; mat_idx < _maturities.size(); ++mat_idx)
	{
		vector<double> lower_diag;
		vector<double> central_diag;
		vector<double> upper_diag;
		vector<double> rhs;

		for (size_t j = 0; j < strike_size() - 3; ++j)
		{
			lower_diag.push_back(_delta_K[j + 1]);
			upper_diag.push_back(_delta_K[j + 1]);
		}

		vector<double> Y;
		for (size_t j = 0; j < strike_size(); ++j)
			Y.push_back(_market_implied_vols[mat_idx][j]);

		for (size_t j = 0; j < strike_size() - 2; ++j)
		{
			central_diag.push_back(2. * (_delta_K[j] + _delta_K[j + 1]));
			rhs.push_back(3. * ((Y[j + 2] - Y[j + 1]) / _delta_K[j + 1] - (Y[j + 1] - Y[j]) / _delta_K[j]));
		}

		ThomasSolver thomas_solver(lower_diag, central_diag, upper_diag, rhs);
		vector<double> solution = thomas_solver.solve();

		vector<double> beta; // maturity T[mat_idx]
		beta.push_back(0.);
		for (size_t j = 1; j < strike_size() - 1; ++j)
			beta.push_back(solution[j - 1]);
		_beta_coefficients.push_back(beta);

		vector<double> alpha; // maturity T[mat_idx]
		for (size_t j = 0; j < _strikes.size() - 2; ++j)
			alpha.push_back((beta[j + 1] - beta[j]) / (3. * _delta_K[j]));
		alpha.push_back(-beta[_strikes.size() - 2] / (3. * _delta_K[_strikes.size() - 2]));
		_alpha_coefficients.push_back(alpha);

		vector<double> gamma; // maturity T[mat_idx]
		for (size_t j = 0; j < _strikes.size() - 1; ++j)
			gamma.push_back((Y[j + 1] - Y[j]) / _delta_K[j] - alpha[j] * _delta_K[j] * _delta_K[j] - beta[j] * _delta_K[j]);
		_gamma_coefficients.push_back(gamma);

	}
}

double ImpliedVolatilitySurface::compute_smile_implied_vol(const size_t & maturity_index, const double & strike) const
{
	double implied_vol;

	// Locate the strike in the right strike interval from _strikes : : K in [K_j, K_{j+1}] or K < K_1 or K > K_N

	// left extrapolation case
	if (strike <= _strikes[0])
	{
		// sigma(T_i, K) = sigma(T_i, K_1) + gamma_1 * (K - K_1)
		double deriv_left = _gamma_coefficients[maturity_index][0];
		implied_vol = _market_implied_vols[maturity_index][0] + deriv_left * (strike - _strikes[0]);
	}
	else
	{
		// right extrapolation case
		if (strike >= _strikes[strike_size() - 1])
		{
			// sigma(T_i, K) = sigma(T_i, K_N) + (3 alpha_{N-1} Delta^2 x_{N-1} + 2 beta_{N-1} Delta x_{N-1} + gamma_{N-1}) * (K - K_N)
			double deriv_right = 3. * _alpha_coefficients[maturity_index][strike_size() - 2] * _delta_K[strike_size() - 2] * _delta_K[strike_size() - 2] + 2. * _beta_coefficients[maturity_index][strike_size() - 2] * _delta_K[strike_size() - 2] + _gamma_coefficients[maturity_index][strike_size() - 2];
			implied_vol = _market_implied_vols[maturity_index][strike_size() - 1] + deriv_right * (strike - _strikes[strike_size() - 1]);
		}
		else
			// interpolation case	
		{
			size_t right_strike_index = 0;
			// the loop stops when K < K_{j+1} so that K_j <= K < K_{j+1} 
			while (strike >= _strikes[right_strike_index])
				right_strike_index++;

			size_t left_strike_index = right_strike_index - 1;
			// Cubic spline index is "left_strike_index"

			// Evaluate that cubic spline function at strike K
			double alpha = _alpha_coefficients[maturity_index][left_strike_index];
			double beta = _beta_coefficients[maturity_index][left_strike_index];
			double gamma = _gamma_coefficients[maturity_index][left_strike_index];

			double dx = strike - _strikes[left_strike_index];
			implied_vol = alpha * dx*dx*dx + beta * dx*dx + gamma * dx + _market_implied_vols[maturity_index][left_strike_index];
		}
	}

	return implied_vol;
}
