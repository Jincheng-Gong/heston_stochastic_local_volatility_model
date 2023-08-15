#include "Blackscholesformulas.h"
#include "Gaussienne.h"
#include <cmath>

double BlackScholesCallvanilla(double Spot,double Strike,double r,double d,double Vol,double Expiry)
{
	double standardDeviation = Vol * sqrt(Expiry);
	double moneyness = log(Spot / Strike);
	double d1 = (moneyness + (r - d) * Expiry +
		0.5 * standardDeviation * standardDeviation) / standardDeviation;
	double d2 = d1 - standardDeviation;
	return Spot * exp(-d * Expiry) * NormalCDF(d1) -Strike * exp(-r * Expiry) * NormalCDF(d2);
}

double BlackScholesCallvanillaVega(double Spot,double Strike,double r,double d,double Vol,double Expiry)
{
	double standardDeviation = Vol * sqrt(Expiry);
	double moneyness = log(Spot / Strike);
	double d1 = (moneyness + (r - d) * Expiry +
		0.5 * standardDeviation * standardDeviation) / standardDeviation;
	return Spot * exp(-d * Expiry) * sqrt(Expiry) * NormalCDF(d1);
}

double BlackScholesCallforwardVega(double Spot, double Strike, double r, double d, double Vol, double Expiry, double start_time, double proportion)
{
	double time_difference = Expiry - start_time;
	double b = r - d;
	double squarred_vol = Vol * Vol;
	double b_plus_squarred_vol = b + squarred_vol / 2;
	return pow(time_difference, 1.0 / 2) + (pow(time_difference, -1.0 / 2) / squarred_vol) * (log(1 / proportion) + b_plus_squarred_vol * time_difference);
}

double BlackScholesCallforward(double Spot, double Strike, double r, double d, double Vol, double Expiry,double start_time,double proportion)
{
	double b = r - d;
	double s_t = Spot * exp(b - r) * start_time;
	double numerator = log(1 / proportion) + (b + Vol * Vol / 2) * (Expiry - start_time);
	double denominator = Vol * sqrt(Expiry - start_time);
	double d1 = numerator / denominator;
	double d2 = d1 - denominator;
	return s_t * (exp((b-r) * (Expiry-start_time)) * NormalCDF(d1) - proportion*exp(-r*(Expiry-start_time)) * NormalCDF(d2));
}
