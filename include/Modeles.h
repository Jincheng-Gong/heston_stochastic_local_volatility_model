#pragma once
#include "DupireLocalVolatilitySurface.h"
#include "Bins.h"
#include <utility>
using Pair = std::pair<double, double>;

// Designed class for models with local volatility like Heston and SABR Local stochastic volatility Models 
class Model_with_vol 
{
protected:
	// Parameters
	double _kappa;
	double _risk_free_rate;
	double _correlation;
	double _vol_of_vol;
	double _theta;
	Pair _init_spot_variance;
	// Dupire Local volatility component
	DupireLocalVolatilitySurface _dupire_local_volatility;
	// Bins object in other to compute the conditional expectation
	bins* bin;
public:
	Model_with_vol(const double& risk_free_rate, const double& correlation, const Pair& init_spot_variance, const double& kappa, const double& vol_of_vol, const double& theta, const DupireLocalVolatilitySurface& dupirelocalvolatility, const bins& b);
	virtual~Model_with_vol();
	Model_with_vol(const Model_with_vol& model);
	Model_with_vol& operator=(const Model_with_vol& model);
	virtual Model_with_vol* clone() const = 0;
	// Gives the drift pair
	Pair drift_pair(const double& time, const Pair& spot_variance);
	// Gives the local volatility component of the model
	double local_volatility(const double& time, vector<Pair>& spot_function_variance, const double& spot) const;
	// Gives the diffusion pair
	virtual Pair diffusion_pair(const double& time, const Pair& spot_variance, vector<Pair>& spot_function_variance) const = 0;
	// virtuals functions
	virtual double psi_function(const double& variance) const = 0;
	virtual double variance_drift(const double& time, const double& variance) const = 0;
	virtual double variance_diffusion(const double& time, const double& variance) const = 0;
	// Getter methods
	double kappa() const;
	double risk_free_rate() const;
	double correlation() const;
	double vol_of_vol() const;
	double theta() const;
	Pair init_spot_variance() const;
	bins* get_bin() const;
	DupireLocalVolatilitySurface get_dupire_vol() const;
};

// Heston local stochastic volatility Model
class Heston_local_sto_vol_Model :public Model_with_vol
{
public:
	Heston_local_sto_vol_Model(const double& risk_free_rate, const double& correlation, const Pair& init_spot_variance, const double& kappa, const double& vol_of_vol, const double& theta, const DupireLocalVolatilitySurface& dupirelocalvolatility, const bins& b);
	Pair diffusion_pair(const double& time, const Pair& spot_variance, vector<Pair>& spot_function_variance) const override;
	double psi_function(const double& variance) const override;
	double variance_drift(const double& time, const double& variance) const override;
	double variance_diffusion(const double& time, const double& variance) const override;
	Heston_local_sto_vol_Model* clone() const override;
};

// Designes class for an european call option
class BSCall
{
public:
	BSCall(double r_, double d_, double T, double Spot_, double Strike_);
	virtual double Price(double Vol) const = 0;
	virtual double Vega(double Vol) const = 0;
protected:
	double r;
	double d;
	double T;
	double Spot;
	double Strike;
};

class BsCallvanilla :public BSCall
{
public:
	BsCallvanilla(double r_, double d_, double T, double Spot_, double Strike_);
	double Price(double vol) const override;
	double Vega(double vol) const override;
};

class BsCallforwardoption :BSCall
{
public:
	BsCallforwardoption(double r_, double d_, double T, double Spot_, double Strike_,double start_time,double proportion);
	double Price(double vol) const override;
	double Vega(double vol) const override;
private:
	double _start_time;
	double _proportion;
};