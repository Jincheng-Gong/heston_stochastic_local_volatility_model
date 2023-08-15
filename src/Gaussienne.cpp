// simulation de la loi gaussienne par Box Muller
#include"Gaussienne.h"
#include<cmath>
#include <iostream>
#include<cmath>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include<cstdlib>
using namespace std;
double Gausienne()
{
	double resultat;
	double x;
	double y;
	double module = 0.0;
	do
	{
		x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		module = x * x + y * y;
	} while (module >= 1.0);
	resultat = x * sqrt(-2 * log(module) / module);
	return resultat;
}


double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = { 2.515517, 0.802853, 0.010328 };
    double d[] = { 1.432788, 0.189269, 0.001308 };
    return t - ((c[2] * t + c[1]) * t + c[0]) /
        (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        /*std::stringstream os;
        os << "Invalid input argument (" << p<< "); must be larger than 0 but less than 1.";
        throw std::invalid_argument(os.str());*/
        cout << "Invalid input argument (" << p << "); must be larger than 0 but less than 1.";
    }

    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation(sqrt(-2.0 * log(p)));
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation(sqrt(-2.0 * log(1 - p)));
    }
}

double NormalCDF(double k)
{
    return erfc(-k / std::sqrt(2)) / 2;
}
