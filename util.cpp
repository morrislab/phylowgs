#include <iostream>
#include <fstream>
#include "util.hpp"

using namespace std;

// stat/math functions
double log_factorial(int n){
	return lgamma(n + 1);
}

double log_bin_coeff(int n, int k){
	return log_factorial(n) - log_factorial(k) - log_factorial(n - k);
}

double log_binomial_likelihood(int x, int n, double mu){
	return  x * log(mu) + (n - x) * log(1 - mu);
}

double log_beta(double a, double b){
	return lgamma(a) + lgamma(b) - lgamma(a + b);
}

void dirichlet_sample(int size, double alpha[], double *x,gsl_rng *r){
	gsl_ran_dirichlet(r,size,alpha,x);
	for(int i=0;i<size;i++)
		x[i]=x[i]+0.0001;
	double sum=0.0;	
	for(int i=0;i<size;i++)
		sum+=x[i];
	for(int i=0;i<size;i++)
		x[i]/=sum;
}

double logsumexp(double x[], int nx){
	double maxes=x[0], sum=0.0;
	for (int i=1;i<nx;i++)
		if (x[i]>maxes)
			maxes=x[i];
	for (int i=0;i<nx;i++)
		sum+=exp(x[i]-maxes);       
	return log(sum) + maxes;
}

