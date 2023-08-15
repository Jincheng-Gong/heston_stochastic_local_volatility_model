#include"Payoff.h"
// #include <minmax.h>

PayOff::PayOff(const double& strike) : _strike(strike)
{
}

double PayOff::strike()
{
	return _strike;
}

PayOffCall::PayOffCall(const double& strike) : PayOff(strike)
{
}

double PayOffCall::operator()(const vector<double> Spot) const
{
	// _strike++;
	return max(Spot[0] - _strike, 0.0);
}

PayOffCall* PayOffCall::clone() const
{
	return new PayOffCall(*this);
}

PayOffPut::PayOffPut(const double& strike): PayOff(strike)
{
}

double PayOffPut::operator()(const vector<double>Spot) const
{
	return max(_strike - Spot[0], 0.0);
}

PayOffPut* PayOffPut::clone() const
{
	return new PayOffPut(*this);
}

PayoffForward::PayoffForward(const double& proportion): PayOff(proportion)
{
}

double PayoffForward::operator()(const vector<double> spot) const
{
	return max(spot[0]-_strike*spot[1],0.0);
}

PayoffForward* PayoffForward::clone() const
{
	return new PayoffForward(*this);
}
