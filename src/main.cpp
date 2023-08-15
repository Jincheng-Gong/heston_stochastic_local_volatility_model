#include "HestonSolver.h"
#include "NewtonRaphson.h"
#include "Pricer.h"
#include "Blackscholesformulas.h"
#include "Payoff.h"

#include <iostream>

using namespace std;

int main() {

    // * Payoff Test
    PayOffCall TestVanillaCall(50);
    vector<double> Spot = {100.0};
    cout << "The Payoff of VanillaCall: " << TestVanillaCall(Spot) << endl;

    // * BS Pricing Test
    BsCallvanilla theCall(0, 0, 5, 100, 100);
    double vol = NewtonRaphson<BsCallvanilla, &BsCallvanilla::Price, &BsCallvanilla::Vega>(21.71805807260404, 0.5, 0.0001, theCall);
    double PriceTwo = BlackScholesCallvanilla(100, 100, 0, 0, vol, 5);
    cout << "The Newton-Raphson Method's Vol: " << vol << endl;
    cout << "Repricing by Black-Scholes Model: " << PriceTwo << endl;

    // * Data Input
    srand(clock());
    Pair init_spot_variance(100, 0.0945);
    equal_number b(100, 10); // equidistant b(100, 5);
    vector<double> s(9);
    vector<double> t(8);
    vector<vector<double>> sigma(8, vector<double>(9));
    s[0] = 20; s[1] = 40; s[2] = 60; s[3] = 80; s[4] = 100; s[5] = 120; s[6] = 140; s[7] = 160; s[8] = 180;
    t[0] = 0.25; t[1] = 0.5; t[2] = 0.75; t[3] = 1; t[4] = 2; t[5] = 3; t[6] = 4; t[7] = 5;
    sigma[0][0] = 0.39; sigma[0][1] = 0.31; sigma[0][2] = 0.24; sigma[0][3] = 0.22; sigma[0][4] = 0.16; sigma[0][5] = 0.19; sigma[0][6] = 0.23; sigma[0][7] = 0.29; sigma[0][8] = 0.38;
    sigma[1][0] = 0.44; sigma[1][1] = 0.36; sigma[1][2] = 0.27; sigma[1][3] = 0.21; sigma[1][4] = 0.17; sigma[1][5] = 0.21; sigma[1][6] = 0.27; sigma[1][7] = 0.35; sigma[1][8] = 0.40;
    sigma[2][0] = 0.45; sigma[2][1] = 0.30; sigma[2][2] = 0.25; sigma[2][3] = 0.21; sigma[2][4] = 0.18; sigma[2][5] = 0.22; sigma[2][6] = 0.29; sigma[2][7] = 0.37; sigma[2][8] = 0.45;
    sigma[3][0] = 0.48; sigma[3][1] = 0.42; sigma[3][2] = 0.34; sigma[3][3] = 0.28; sigma[3][4] = 0.20; sigma[3][5] = 0.26; sigma[3][6] = 0.31; sigma[3][7] = 0.42; sigma[3][8] = 0.50;
    sigma[4][0] = 0.52; sigma[4][1] = 0.43; sigma[4][2] = 0.34; sigma[4][3] = 0.26; sigma[4][4] = 0.21; sigma[4][5] = 0.27; sigma[4][6] = 0.38; sigma[4][7] = 0.45; sigma[4][8] = 0.55;
    sigma[5][0] = 0.54; sigma[5][1] = 0.46; sigma[5][2] = 0.34; sigma[5][3] = 0.27; sigma[5][4] = 0.23; sigma[5][5] = 0.28; sigma[5][6] = 0.36; sigma[5][7] = 0.49; sigma[5][8] = 0.58;
    sigma[6][0] = 0.57; sigma[6][1] = 0.50; sigma[6][2] = 0.46; sigma[6][3] = 0.35; sigma[6][4] = 0.25; sigma[6][5] = 0.32; sigma[6][6] = 0.45; sigma[6][7] = 0.54; sigma[6][8] = 0.60;
    sigma[7][0] = 0.60; sigma[7][1] = 0.52; sigma[7][2] = 0.41; sigma[7][3] = 0.31; sigma[7][4] = 0.26; sigma[7][5] = 0.34; sigma[7][6] = 0.40; sigma[7][7] = 0.55; sigma[7][8] = 0.62;

    // * Implied Volatility Surface
    ImpliedVolatilitySurface surface(t, s, sigma, 0.01);
    cout << "Risk Free Rate: " << surface.risk_free_rate() << endl;
    cout << "Implied Volatility (0.79, 157.310665): " << surface.implied_volatility(0.79, 156.310665) << endl;

    // * Dupire Local Volatility Surface
    DupireLocalVolatilitySurface dlvs(surface, 0.0001, 0.0001, 100);
    cout << "Risk Free Rate: " << dlvs.risk_free_rate() << endl;
    cout << "Local Volatility (0.8, 158.310665): " << dlvs.local_volatility(0.8, 158.310665) << endl;

    // * Heston Local Stochastic Volatility Model
    Heston_local_sto_vol_Model model(0, -0.315, init_spot_variance, 1.05, 0.95, 0.0855, dlvs, b);
    vector<double> v{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};

    // * HestonSLV Monte Carlo Method
    // PathSimulatorEuler p(model, v); 
    PathSimulatorQE p(model, v);
    MonteCarloPricer MC(TestVanillaCall, 0.5, p, 10000); // MonteCarloPricer MC(*TestVanillaCall.clone(), 1, p, 1000);
    double mc_price = MC.price();
    cout << "Monte Carlo Price for VanillaCall: " << mc_price << endl;

    // * BS Pricing Method with Implied Volatility
    double vol_iv = surface.implied_volatility(1, 50);
    double PriceBS_iv = BlackScholesCallvanilla(100, 50, 0.01, 0, vol_iv, 1);
    cout << "BS Vol: " << vol_iv << endl;
    cout << "Repricing by Black-Scholes Model with iv: " << PriceBS_iv << endl;

    // * BS Pricing Method with Dupire Local Volatility
    double vol_lv = dlvs.local_volatility(1, 50);
    double PriceBS_lv = BlackScholesCallvanilla(100, 50, 0.01, 0, vol_lv, 1);
    cout << "Local Vol: " << vol_lv << endl;
    cout << "Repricing by Black-Scholes Model with lv: " << PriceBS_lv << endl;

    // * Program Over
    return 0;
}
