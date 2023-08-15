#include "PathSimulator.h"
#include <algorithm>
#include <cstddef>
// #include <ppl.h> // 微软并行库

PathSimulator::PathSimulator(const Model_with_vol& model, const vector<double> time_points) :model(model.clone()), _time_points(time_points)
{
}

PathSimulator::PathSimulator(const PathSimulator& path_simulator)
	: model((path_simulator.model)->clone()), _time_points(path_simulator._time_points)
{
}

PathSimulator& PathSimulator::operator=(const PathSimulator& path_simulator)
{
	if (this != &path_simulator)
	{
		_time_points = path_simulator._time_points;
		delete model;
		model = path_simulator.model->clone();
	}
	return *this;
}

PathSimulator::~PathSimulator()
{
	delete model;
}

// Pair PathSimulator::path_maturity()
// {
// 	size_t number_of_simulations = model->get_bin()->number_of_simulations();
// 	Pair spot_variance = model->init_spot_variance();
// 	vector<Pair> samples(number_of_simulations, spot_variance);
// 	vector<Pair> v;
// 	size_t k{};
// 	for (size_t i = 1; i < _time_points.size(); i++)
// 	{
// 		v = samples;
// 		for (size_t j = 0; j < number_of_simulations && i != _time_points.size() - 1; j++)
// 		{
// 			samples[j] = next_step(i, spot_variance, v);
// 		}
// 		for (k = 0; k < number_of_simulations && samples[k].second == 0 && i != _time_points.size() - 1; k++);
// 		spot_variance = samples[k]; 
// 		if (i == _time_points.size() - 1)
// 			spot_variance = next_step(i, spot_variance, v);
// 	}
// 	return spot_variance;
// }

Pair PathSimulator::path_maturity()
{ // TODO Path Maturity (DONE)
	size_t number_of_simulations = model->get_bin()->number_of_simulations();
	Pair spot_variance = model->init_spot_variance();
	vector<Pair> samples(number_of_simulations, spot_variance);
	vector<Pair> v;
	size_t k{};
	for (size_t i = 1; i < _time_points.size(); i++)
	{
		v = samples;
		for (size_t j = 0; j < number_of_simulations && i != _time_points.size() - 1; j++)
		{
			spot_variance = samples[j];
			samples[j] = next_step(i, spot_variance, v);
		}
		if (i == _time_points.size() - 1){
			for (size_t j = 0; j < number_of_simulations; j++){
				spot_variance.first += samples[j].first;
				spot_variance.second += samples[j].second;
			}
			spot_variance.first /= number_of_simulations;
			spot_variance.second /= number_of_simulations;
			spot_variance = next_step(i, spot_variance, v);
		}
	}
	return spot_variance;
}

Model_with_vol* PathSimulator::get_model()
{
	return model;
}

size_t PathSimulator::index_maturity() const
{
	return _time_points.size();
}

PathSimulator* PathSimulatorEuler::clone() const
{
	return new PathSimulatorEuler(*this);
}

PathSimulatorEuler::PathSimulatorEuler(const Heston_local_sto_vol_Model& model, const vector<double> time_points) :PathSimulator(model, time_points)
{
}

Pair PathSimulatorEuler::next_step(const size_t& time_idx, const Pair& spot_variance, vector<Pair>& spot_variance_sample) const
{ // TODO Next Step (PASS with BS Vol)
	if (time_idx == 0)
		return model->init_spot_variance();

	double next1, next2;
	//Gives the discretization step
	double delta_t = _time_points[time_idx] - _time_points[time_idx - 1];
	// Asset's diffusion brownian
	double z1 = Gausienne();
	// Volatility's diffusion Brownian
	double z2 = Gausienne();
	double z3 = model->correlation() * z1 + sqrt(1 - pow(model->correlation(), 2)) * z2;
	//Gives the volatility discretization 
	next2 = spot_variance.second + model->variance_drift(_time_points[time_idx - 1], max(spot_variance.second, 0.0)) * delta_t + model->variance_diffusion(_time_points[time_idx - 1], max(spot_variance.second, 0.0)) * sqrt(delta_t) * z3;
	// To avoid negative values of the variance (Full Truncation)
	next2 = max(next2, 0.0);
	// Gives the asset discretization
	if (time_idx == 1)
	{
		next1 = spot_variance.first + model->risk_free_rate() * delta_t * spot_variance.first + spot_variance.second * spot_variance.first * model->psi_function(max(spot_variance.second, 0.0)) * sqrt(delta_t) * z1;
		return Pair(max(next1,10e-100), next2);
	}

	// Applying psi_function to variance's sample
	for (size_t i = 0; i < spot_variance_sample.size(); i++)
		spot_variance_sample[i].second = pow(model->psi_function(spot_variance_sample[i].second), 2);

	next1 = spot_variance.first + model->risk_free_rate() * delta_t * spot_variance.first + sqrt(model->local_volatility(_time_points[time_idx - 1], spot_variance_sample, spot_variance.first)) * spot_variance.first * model->psi_function(max(spot_variance.second, 0.0)) * sqrt(delta_t) * z1;

	// if (isnan(next1)) { // 如果出现local volatility是nan的情况，用上一步的volatility（spot_variance.second）计算下一步spot
	// 	next1 = spot_variance.first + model->risk_free_rate() * delta_t * spot_variance.first + spot_variance.second * spot_variance.first * model->psi_function(max(spot_variance.second, 0.0)) * sqrt(delta_t) * z1;
	// 	double test = 0;
	// 	test = sqrt(model->local_volatility(_time_points[time_idx - 1], spot_variance_sample, spot_variance.first));
	// }

	return Pair(max(next1, 10e-100), next2);
}

PathSimulatorQE* PathSimulatorQE::clone() const
{
	return new PathSimulatorQE(*this);
}

PathSimulatorQE::PathSimulatorQE(const Model_with_vol& model, const vector<double> time_points) :PathSimulator(model, time_points)
{
}

Pair PathSimulatorQE::next_step(const size_t& time_idx, const Pair& spot_variance, vector<Pair>& spot_variance_sample) const
{
	if (time_idx == 0)
		return model->init_spot_variance();
	//Gives the discretization step
	double delta_t = _time_points[time_idx] - _time_points[time_idx - 1];
	// Brownian of the volatility
	double z1 = Gausienne();
	// QE algorithm to simulate de volatility of the Model
	double m = model->theta() + (spot_variance.second - model->theta()) * exp(-model->kappa() * delta_t);
	double s1 = spot_variance.second * pow(model->vol_of_vol(), 2) * exp(-model->kappa() * delta_t) / model->kappa();
	double s2 = model->theta() * pow(model->vol_of_vol(), 2) / (2 * model->kappa());
	double s3 = 1 - exp(-model->kappa() * delta_t);
	double s_quarred = s1 * s3 + s2 * pow(s3, 2);
	double psi = s_quarred / pow(m, 2);
	double u_v = (rand() % RAND_MAX + 1.) / (RAND_MAX + 1.);
	double next2{};
	if (psi <= 1.5)
	{
		double b = 2 * 1 / psi - 1 + sqrt(2 / psi) * sqrt(2 / psi - 1);
		double a = m / (1 + pow(b, 2));
		double z_v = NormalCDFInverse(u_v);
		next2 = a * pow(b + z_v, 2);
	}
	else
	{
		double p = (psi - 1) / (psi + 1);
		double beta = (1 - p) / m;
		if (u_v >= 0 && u_v <= p)
			next2 = 0;
		else if (u_v > p && u_v <= 1)
			next2 = (1 / beta) * log((1 - p) / (1 - u_v));
	}
	// Bownian of the asset's diffusion
	double z2 = Gausienne();
	double z3 = model->correlation() * z1 + sqrt(1 - pow(model->correlation(), 2)) * z2;
	// Transforming the stock in log in other to simulate
	double current_log_spot = log(spot_variance.first);
	double c1 = model->kappa() * delta_t - 1;
	double rho1 = sqrt(1 - pow(model->correlation(), 2));
	// QE algorithm to simulate the stock
	if (time_idx == 1)
	{
		double r1 = next2 - model->kappa() * model->theta() * delta_t + spot_variance.second * c1;
		double r2 = sqrt(spot_variance.second * delta_t);
		double next_log_spot = current_log_spot + model->risk_free_rate() * delta_t - 0.5 * spot_variance.second * delta_t + (model->correlation() / model->vol_of_vol()) * r1 + rho1 * r2 * z3;
		double next1 = exp(next_log_spot);
		return Pair(next1, next2);
	}
	// Applying psi function to the variance samples
	for (size_t i = 0; i < spot_variance_sample.size(); i++)
		spot_variance_sample[i].second = pow(model->psi_function(spot_variance_sample[i].second), 2);

	double r1 = sqrt(model->local_volatility(_time_points[time_idx - 1], spot_variance_sample, spot_variance.first)) * (next2 - model->kappa() * model->theta() * delta_t + spot_variance.second * c1);
	double r2 = sqrt(model->local_volatility(_time_points[time_idx - 1], spot_variance_sample, spot_variance.first) * spot_variance.second * delta_t);
	double next_log_spot = current_log_spot + model->risk_free_rate() * delta_t - 0.5 * model->local_volatility(_time_points[time_idx - 1], spot_variance_sample, spot_variance.first) * spot_variance.second * delta_t + (model->correlation() / model->vol_of_vol()) * r1 + rho1 * r2 * z3;
	double next1 = exp(next_log_spot);
	return Pair(next1, next2);
}
