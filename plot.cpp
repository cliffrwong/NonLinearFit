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
#include "plot.h"
#include "gnuplot-iostream.h"

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



Gnuplot gp;


Plot::Plot()
{
    // gp = {};
    char const* imagefp = "figure.png";
    char buff[256];
    sprintf(buff, "set term pngcairo size 800, %d\n",15*2*300);
    gp << buff;
    int num = sprintf(buff, "set output \"%s\"\n", imagefp);
    if (num >= 255 || num < 0)
        printf("Buffer overflow\n");
    gp << buff;
    sprintf(buff, "set multiplot layout %d,1\nset lmargin at screen 0.15\nset "
                  "rmargin at screen 0.95\n",2*15);
    gp << buff;
}

Plot::~Plot()
{

}

void Plot::plotwork(vector<double> &x, vector<double> &y, bool printLines)
{
	char buff[256];
    double max;
    double min;
    double ytics;
    minMaxVal(y, max, min, ytics);
    string printLineStr;
    if(printLines)
        printLineStr = "lines";
    else
        printLineStr = "points";
    sprintf(buff, "set xrange [%f:%f]\nset yrange [%f:%f]\nset grid ytics lt 0"
    			  " lw 1 lc rgb \"#bbbbbb\"\nset ytics %f\nplot '-' using 1:2 "
                  "with %s\n", x.front(), x.back(), min, max, ytics, \
                  printLineStr.c_str());
	gp << buff;
    gp.send1d(boost::make_tuple(x, y));
}

void Plot::minMaxVal(vector<double> &v, double &max, double &min, double &ytics)
{
	double margin = .05;
    max = v[0];       // start with max = first element
    min = v[0];       // start with min = first element

    for(uint i = 1; i<v.size(); i++)
    {
    	if(v[i] > max)
        	max = v[i];
     	else if(v[i] < min)
     		min = v[i];
     }           // return highest value in array
     int range = max - min;
     ytics = range/5;
     max = max+range*margin;
     min = min-range*margin;
}



void Plot::fittedcurve (vector<double> &x, vector<double> &yfit, double &A, double &mu, double &sig, double &y0)
{
	for (uint i = 0; i < x.size(); ++i){
		yfit.push_back(A * exp (-(x[i]-mu)*(x[i]-mu)/(2*(sig*sig))) + y0);
    }
}

