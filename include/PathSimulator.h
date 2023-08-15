#pragma once
#include "Modeles.h"
#include "Gaussienne.h"
#include <random>
#include <chrono>

// Designed class for simulating the models
class PathSimulator
{
public:
	virtual PathSimulator* clone() const = 0;
	PathSimulator(const Model_with_vol& model, const vector<double> time_points);
	PathSimulator(const PathSimulator& path_simulator);
	PathSimulator& operator=(const PathSimulator& path_simulator);
	virtual ~PathSimulator();
	// draw the sample for the computation of conditionnal expectation
	virtual Pair next_step(const size_t& time_idx, const Pair& spot_variance, vector<Pair>& spot_variance_sample) const = 0;
	// gives maturity path
	Pair path_maturity();
	// Getter Method
	Model_with_vol *get_model();
	size_t index_maturity() const;
protected:
	// Model to simulate
	Model_with_vol* model;
	// Time's points discretization
	vector<double> _time_points; // [t_0 = 0, t_1, ..., t_M]
};

// Euler's scheme for discretization
class PathSimulatorEuler : public PathSimulator
{
public:
	PathSimulator* clone() const override;
	PathSimulatorEuler(const Heston_local_sto_vol_Model& model, const vector<double> time_points);
	Pair next_step(const size_t& time_idx, const Pair& spot_variance, vector<Pair>& spot_variance_sample) const override;
};

// Scheme QE for discretization
class PathSimulatorQE : public PathSimulator
{
public:
	PathSimulatorQE* clone() const override;
	PathSimulatorQE(const Model_with_vol& model, const vector<double> time_points);
	Pair next_step(const size_t& time_idx, const Pair& spot_variance, vector<Pair>& spot_variance_sample) const override;
};
