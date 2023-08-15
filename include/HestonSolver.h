#include<complex>
#include <utility>
#pragma once
using Pair = std::pair<double,double>;
double HestonSolver(const double& tau, const double& kappa, const double& theta, const double& vol_of_vol, const double& rho, const double& strike, const double& init_spot, const double& init_variance, const double& risk_free_rate);
