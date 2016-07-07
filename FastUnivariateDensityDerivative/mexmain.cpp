//-------------------------------------------------------------------
// The code was written by Vikas C. Raykar 
// and is copyrighted under the Lessr GPL: 
//
// Copyright (C) 2006 Vikas C. Raykar
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 or later.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details. 
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
// MA 02111-1307, USA.  
//
// The author may be contacted via email at: vikas(at)cs(.)umd(.)edu 
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// File    : mexmain.cpp
// Purpose : Interface between MATLAB and C++
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : September 17, 2005
//-------------------------------------------------------------------

#include "mex.h"
#include "UnivariateDensityDerivative.h"

//The gateway function

void mexFunction(int nlhs,				// Number of left hand side (output) arguments
				 mxArray *plhs[],		// Array of left hand side arguments
				 int nrhs,              // Number of right hand side (input) arguments
				 const mxArray *prhs[])  // Array of right hand side arguments
{

   //check for proper number of arguments 
 
    if(nrhs != 7) mexErrMsgTxt("7 input  arguments required.");
	if(nlhs != 1) mexErrMsgTxt("1 output argument  required.");

   //////////////////////////////////////////////////////////////
  // Input arguments
  //////////////////////////////////////////////////////////////
  
 
  //------ the first input argument: NSources ---------------//

  int argu = 0;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'N' must be a scalar.");
  }

  /*  get the scalar input NSources */
  int NSources = (int) mxGetScalar(prhs[argu]);
  if (NSources <= 0) mexErrMsgTxt("Input 'N' must be a positive number.");

   //------ the second input argument: MTargets ---------------//

  argu = 1;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'M' must be a scalar.");
  }

  /*  get the scalar input NSources */
  int MTargets = (int) mxGetScalar(prhs[argu]);
  if (MTargets <= 0) mexErrMsgTxt("Input 'M' must be a positive number.");

  //----- the third input argument: pSources--------------//
 
  argu = 2;

  /*  create a pointer to the input vector pSources */
  double *pSources = mxGetPr(prhs[argu]);
  
  int mrows = mxGetM(prhs[argu]); //mrows
  int ncols = mxGetN(prhs[argu]); //ncols
  if ( mrows != 1 && ncols != NSources)  mexErrMsgTxt("Input 'X' must be a 1 x N matrix");

   //----- the fourth input argument: pTargets--------------//
  argu = 3;

  /*  create a pointer to the input vector pTargets */
  double *pTargets = mxGetPr(prhs[argu]);
 
  mrows = mxGetM(prhs[argu]); //mrows
  ncols = mxGetN(prhs[argu]); //ncols
  if ( mrows != 1 && ncols != MTargets)  mexErrMsgTxt("Input 'Y' must be a 1 x M matrix");

  //----- the fifth input argument: Bandwidth--------------//

  argu = 4;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'h' must be a scalar.");
  }

  double  Bandwidth = (double) mxGetScalar(prhs[argu]);
  if ( Bandwidth <= 0.0) mexErrMsgTxt("Input 'h' must be a positive number.");

  //------ the sixth input argument: Order ---------------//

  argu = 5;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'r' must be a scalar.");
  }

  /*  get the scalar input NSources */
  int Order = (int) mxGetScalar(prhs[argu]);
  if (Order < 0) mexErrMsgTxt("Input 'r' must be a positive integer.");

  //----- the seventh input argument: epsilon--------------//

  argu = 6;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'epsilon' must be a scalar.");
  }

  double  epsilon = (double) mxGetScalar(prhs[argu]);
  if ( epsilon <= 0.0) mexErrMsgTxt("Input 'epsilon' must be a positive number.");

  

  //////////////////////////////////////////////////////////////
  // Output arguments
  //////////////////////////////////////////////////////////////

   /*  set the output pointer to the output result(vector) */
  plhs[0] = mxCreateDoubleMatrix(1,MTargets,mxREAL);
  
  /*  create a C pointer to a copy of the output result(vector)*/
  double *pDensityDerivative = mxGetPr(plhs[0]);

  //////////////////////////////////////////////////////////////
  // function calls;
  //////////////////////////////////////////////////////////////

  UnivariateDensityDerivative* pUDD=new UnivariateDensityDerivative(
	        NSources,
			MTargets,
			pSources,
			pTargets,
			Bandwidth,
		    Order,
			epsilon,
			pDensityDerivative);

  pUDD->Evaluate();

  delete pUDD;

  return;
  
}
