#include "DupireLocalVolatilitySurface.h"
#include<cmath>
#include <iostream>
using namespace std;
DupireLocalVolatilitySurface::DupireLocalVolatilitySurface(const ImpliedVolatilitySurface & implied_vol_surface, const double& eps_mat, const double& eps_strike, const double& init_spot)
	: _implied_volatility_surface(implied_vol_surface), _epsilon_maturity(eps_mat), _epsilon_strike(eps_strike), _initial_spot(init_spot)
{
}

// return sigma_D(T,K)
double DupireLocalVolatilitySurface::local_volatility(const double & maturity, const double & strike) const
{ // TODO Local Volatility (PASS with BS Vol)
	// sigma*(T,K)
	double sigma = _implied_volatility_surface.implied_volatility(maturity, strike);

	double d_sigma_dT = first_order_derivative_impliedvol_maturity(maturity, strike);
	double d_sigma_dK = first_order_derivative_impliedvol_strike(maturity, strike);
	double d2_sigma_dK2 = second_order_derivative_impliedvol_strike(maturity, strike);
	double risk_free_rate = _implied_volatility_surface.risk_free_rate();

	double d1 = (log(_initial_spot / strike) + (risk_free_rate + 0.5 * sigma * sigma)*maturity) / (sigma* sqrt(maturity));
	double d2 = d1 - sigma * sqrt(maturity);
	
	double numerator = 1. + (2.*maturity / sigma) * (d_sigma_dT + risk_free_rate * d_sigma_dK);

	double k_s_sqrtT = strike * sigma * sqrt(maturity);
	double k_ds_dK_sqrtT = strike * d_sigma_dK * sqrt(maturity);
	double k_d2s_dK2_sqrtT = strike * d2_sigma_dK2 * sqrt(maturity);

	double denominator = 1. + 2. * d1 * k_ds_dK_sqrtT + d1*d2*k_ds_dK_sqrtT*k_ds_dK_sqrtT + k_d2s_dK2_sqrtT * k_s_sqrtT;

	// if (denominator < 0){ // abs
	// 	denominator = 1;
	// 	// cout << denominator << " ";
	// }

	// if (numerator < 0) { // abs
	// 	numerator = 1;
	// 	// cout << numerator << " ";
	// }

	// denominator < 0 or numerator < 0 ---> dupire_volatility < 0
	// When lv < 0, use bs_vol.
	if (denominator < 0 || numerator < 0) {
		denominator = 1;
		numerator = 1;
	}

	double dupire_variance = sigma * sigma * numerator / denominator;
	double dupire_volatility = sqrt(dupire_variance);

	return dupire_volatility;
}

double DupireLocalVolatilitySurface::first_order_derivative_impliedvol_maturity(const double & maturity, const double & strike) const
{
	double sigma_up = _implied_volatility_surface.implied_volatility(maturity + _epsilon_maturity, strike);
	double sigma_down = _implied_volatility_surface.implied_volatility(maturity - _epsilon_maturity, strike);
	double deriv = 0.5 * (sigma_up - sigma_down) / (_epsilon_maturity);
	return deriv;
}

double DupireLocalVolatilitySurface::first_order_derivative_impliedvol_strike(const double & maturity, const double & strike) const
{
	double sigma_up = _implied_volatility_surface.implied_volatility(maturity, strike + _epsilon_strike);
	double sigma_down = _implied_volatility_surface.implied_volatility(maturity, strike - _epsilon_strike);
	double deriv = 0.5 * (sigma_up - sigma_down) / (_epsilon_strike);
	return deriv;
}

double DupireLocalVolatilitySurface::second_order_derivative_impliedvol_strike(const double & maturity, const double & strike) const
{
	double sigma_up = _implied_volatility_surface.implied_volatility(maturity, strike + _epsilon_strike);
	double sigma_base = _implied_volatility_surface.implied_volatility(maturity, strike);
	double sigma_down = _implied_volatility_surface.implied_volatility(maturity, strike - _epsilon_strike);
	double deriv = (sigma_up - 2.*sigma_base + sigma_down) / (_epsilon_strike*_epsilon_strike);
	return deriv;
}
