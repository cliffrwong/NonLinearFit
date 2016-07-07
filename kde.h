#ifndef __kdegauss__kde__
#define __kdegauss__kde__

#include <boost/function.hpp>
#include <cmath>
#include <vector>

using namespace std;


class KernelDensityEstimate
{
	public:
		//constructor
		KernelDensityEstimate();

		~KernelDensityEstimate();
		
		double KernelSmoother (const boost::function<double (double)>& kernel, const std::vector<float>& array, double x, double h);

		static double GaussianKernel(double x);	

	private:

};

#endif