#pragma once
#include <vector>
using namespace std;

// Designed class for options's payoffs construction
class PayOff
{
public:
	PayOff(const double & strike);
	virtual double operator()(const vector<double> Spot) const = 0;// Gives the payoff
	virtual PayOff* clone() const = 0;
	virtual ~PayOff() {}
	//Getter method
	double strike();
protected:
	// strike of the option
	double _strike;
};

// Vanilla call's Payoff
class PayOffCall : public PayOff
{
public:
	PayOffCall(const double& Strike_);
	double operator()(const vector<double> Spot) const override;
	virtual ~PayOffCall() {}
	virtual PayOffCall* clone() const override;
};

// Vanilla put's payoff
class PayOffPut : public PayOff
{
public:
	PayOffPut(const double& Strike_);
	double operator()(const vector<double> Spot) const override;
	virtual ~PayOffPut() {}
	PayOffPut* clone() const override;
};

class PayoffForward:public PayOff
{
public:
	PayoffForward(const double& proportion);
	double operator()(const vector<double> spot) const override;
	virtual ~PayoffForward() {}
	PayoffForward* clone() const override;
};
