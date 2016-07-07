/////////////////////////////////////////////////////////////////////////////// 
// 
// Copyright (c) 2016 Cliff Wong 
// 
// This code is licensed under the MIT License (MIT). 
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
// THE SOFTWARE. 
// 
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <malloc.h>
#include <numeric>
#include <csignal>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <fstream>
#include "gaussfit.h"
#include "plot.h"
#include "kde.h"

using namespace std;


//-------------------------------------------------------------------
// Constructor.
//
// PURPOSE                                                    
// -------   
// Initialize the class. 
// Read the parameters.
//
// PARAMETERS                                                      
// ----------
// NSources		      --> number of sources, N.
// MTargets		      --> number of targets, M.
// pSources		      --> pointer to sources, px(N).
// pTargets             --> pointer to the targets, py(M).
// Bandwidth		--> the source bandwidth, h.
// Order                --> order of the derivative, r.  
// pDensityDerivative   --> pointer the the evaluated Density 
//-------------------------------------------------------------------



// void generate(vector<int> &xcenter2, vector<double> &bdiff2, double &totalsc, double &sig, double &mu, Gnuplot &gp, double &max_error, bool &saveImg){
// 	uint N = 100;
// 	double A, y0;
// 	vector<double> yfit(N);
// 	gaussfit(bdiff2, A, mu, sig, y0);
// 	sig = abs(sig);
// 	fittedcurve(xcenter2, yfit, A, mu, sig, y0);
// 	//plot yfit
// 	if(saveImg)
// 		gp.send1d(boost::make_tuple(xcenter2, yfit));
// 	// double height = A - y0;
// 	// cout << "height " << height << " mu " << mu << " sig " << sig << " adjerror " << adjerror << endl;
// }

void kgpFunction(int argc, char *argv[]) {
	// ofstream myfile;
	// myfile.open ("data.txt");
	// myfile.close();
	fstream myfile("data.txt", ios_base::in);
	float a;
	vector<float> array;
    while (myfile >> a)
    {
    	array.push_back(a);
        // printf("%f ", a);
    }
    int size = 30;
    vector<double> x, y, yfit;
    double A, mu, sig, y0;
    double min = 0;
    double max = 30;
    double incr = (max - min)/size;
    KernelDensityEstimate kde;
    Plot plot;
    for (int i = 1; i < size; i++) {
    	x.push_back(i*incr+min);
    	y.push_back(kde.KernelSmoother(KernelDensityEstimate::GaussianKernel, array, i*incr+min, incr/5));
    }

    gaussfit(x, y, A, mu, sig, y0);

    plot.fittedcurve(x, yfit, A, mu, sig, y0);

	plot.plotwork(x, y, false);
	plot.plotwork(x, yfit, true);
}


//The gateway function
int main(int argc, char *argv[]) {
	kgpFunction(argc, argv);
	return 0;
}