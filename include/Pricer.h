#pragma once
#include "Payoff.h"
#include "PathSimulator.h"


// Monte Carlo way of computing the Price of options in a Model
class MonteCarloPricer 
{
public:
	MonteCarloPricer(const PayOff& payoff,const double& maturity, const PathSimulator& path1,const size_t& number_of_simulations);
	MonteCarloPricer* clone() const;
	MonteCarloPricer(const MonteCarloPricer& pricer);
	MonteCarloPricer& operator=(const MonteCarloPricer& pricer);
	double price() const;
	~MonteCarloPricer();
private:
	PathSimulator* _path;// The simulation scheme choosed
	int _number_of_simulations;
	PayOff* _payoff;// the option Payoff
	double _maturity;// the maturity of the option
};
