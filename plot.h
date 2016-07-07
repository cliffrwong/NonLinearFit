//  filesyncalgo.hpp
//  subfilesync
//
//  Created by Cliff Wong on 10/5/14.
//  Copyright (c) 2014 aosinnov. All rights reserved.
//

#ifndef __wavscore__plot__
#define __wavscore__plot__

#include <boost/function.hpp>
#include <cmath>
#include <vector>

using namespace std;

class Plot
{
	public:
		Plot();

		~Plot();

		void plotwork(vector<double> &x, vector<double> &y, bool printLines);

		void fittedcurve (vector<double> &x, vector<double> &yfit, double &A, double &mu, double &sig, double &y0);

	private:
		void minMaxVal(vector<double> &v, double &max, double &min, double &ytics);


};

#endif 