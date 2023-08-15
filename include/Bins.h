#pragma once
// Designed class for the computation of conditionnal expectation with the bins
#include<vector>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


class bins
{
protected:
	// Number of simulations for the computation of conditional expectation
	size_t _number_of_simulations;
	// Number of bins for the computation of conditional expectation
	size_t _number_of_bins;
	// Create the bins from a sample according to the choosen method
	virtual vector<double> bins_creation_manner(const vector<double>& spot) const = 0;
public:
	// constructor with parameters
	bins(const size_t & number_of_simulations, const size_t& number_of_bins);
	// Virtual destructor
	virtual ~bins() {};
	// Enables the object to makes copies of it's own
	virtual bins* clone() const = 0;
	//Computes the conditional expectation according to the bin to which the simulated spot value belongs
	double expectation(vector< pair<double, double> >& _spot_function_variance, const double &spot) const;
	// Getter 
	size_t number_of_simulations();
};

// Derived class from bins for the equidistant method of creating bins
class equidistant : public bins
{
protected:
	vector<double> bins_creation_manner(const vector<double>& spot) const override;
public:
	equidistant(const int& number_of_simulations, const int& number_of_bins);
	equidistant* clone() const override;
};

// Derived class from bins for the equal number method of creating bins
class equal_number : public bins
{
protected:
	vector<double> bins_creation_manner(const vector<double>& spot) const override;
public:
	equal_number(const int& number_of_simulations, const int& number_of_bins);
	equal_number* clone() const override;
};
