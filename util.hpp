

#include <iostream>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

// stat/math functions
double log_factorial(int n);
double log_bin_coeff(int n, int k);
double log_binomial_likelihood(int x, int n, double mu);
double log_beta(double a, double b);
void dirichlet_sample(int size, double alpha[], double *x,gsl_rng *r);
double logsumexp(double x[], int nx);