#pragma once
// This header file contains all the functions relative to the gaussian distribution needed in other to performs the computations
double Gausienne();// draw a random number by the Box-Muller method
double RationalApproximation(double t);// useful method to compute the inverse of normal cumulative distribution function
double NormalCDFInverse(double p);// Compute the inverse of normal cumulative distribution function
double NormalCDF(double k);//Compute the normal cumulative distribution function
