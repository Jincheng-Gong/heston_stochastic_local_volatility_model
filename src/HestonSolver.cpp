#include "HestonSolver.h"
#include <iostream>
#include<cmath>
using namespace std;

double HestonSolver(const double& tau, const double& kappa, const double& theta, const double& vol_of_vol, const double& rho, const double& strike, const double & init_spot ,const double& init_variance, const double& risk_free_rate)
{
    int N = 10000, N_max = 1000, i;
    double aa, bb, u2, du = static_cast<double>(N_max) / N;
    complex<double> I(0, 1), P(0, 0), u1(0, 0), a1(0, 0), a2(0, 0), kap(kappa, 0), lamb2(vol_of_vol * vol_of_vol, 0), d1(0, 0),
        d2(0, 0), g1(0, 0), g2(0, 0), un(1, 0), t(-tau, 0), con(log(init_spot / strike) + risk_free_rate * tau, 0),
        v0(init_variance, 0), b1(0, 0), b2(0, 0), phi1(0, 0), phi2(0, 0), a(theta * kappa * tau / pow(vol_of_vol, 2), 0),
        d(du, 0), st(strike, 0), deux(2, 0), con1(init_spot / strike - exp(-risk_free_rate * tau)), pi(3.14, 0), rho_vol(rho * vol_of_vol, 0);
    for (i = 1; i < N; i++)
    {
        aa = theta * kappa * tau / pow(vol_of_vol, 2);
        bb = -2 * theta * kappa / pow(vol_of_vol, 2);
        u2 = i * du;
        u1 = complex<float>(u2, -1);
        a1 = rho_vol * u1 * I;
        a2 = rho_vol * u2 * I;
        d1 = sqrt((a1 - kap) * (a1 - kap) + lamb2 * (u1 * I + u1 * u1));
        d2 = sqrt((a2 - kap) * (a2 - kap) + lamb2 * (u2 * I + u2 * u2));
        g1 = (kap - a1 - d1) / (kap - a1 + d1);
        g2 = (kap - a2 - d2) / (kap - a2 + d2);
        b1 = exp(u1 * I * (con)) * pow(((un - g1 * exp(d1 * t)) / (un - g1)), bb);
        b2 = exp(u2 * I * (con)) * pow(((un - g2 * exp(d2 * t)) / (un - g2)), bb);
        phi1 = b1 * exp(a * (kap - a1 - d1) + v0 * (kap - a1 - d1) * (un - exp(d1 * t)) / (un - g1 * exp(d1 * t)) / lamb2);
        phi2 = b2 * exp(a * (kap - a2 - d2) + v0 * (kap - a2 - d2) * (un - exp(d2 * t)) / (un - g2 * exp(d2 * t)) / lamb2);
        P += ((phi1 - phi2) / (u2 * I)) * d;
    }
    return strike * real((con1) / deux + P / pi);
}
