#include <stdio.h>
#include <malloc.h>
#include <numeric>
#include <csignal>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <boost/bind.hpp>
#include "kde.h"

using namespace std;


double KernelDensityEstimate::GaussianKernel(double x) { 
	// N(0,1) variable
	double a = 1;
	return a * exp( - 0.5 * x * x);
}

KernelDensityEstimate::KernelDensityEstimate() {

}

KernelDensityEstimate::~KernelDensityEstimate() {

}

double KernelDensityEstimate::KernelSmoother (const boost::function<double (double)>& kernel, const vector<float>& array, double x, double h) {
	// The sum of the terms of the kernel density function
	h = 1.0 / h;
	long n = array.size();
	double result = kernel(h * (x - array[0]));
	for (long j = 1; j < n; ++j) {
		result += kernel(h * (x - array[j]));
	}
	return result;
}
