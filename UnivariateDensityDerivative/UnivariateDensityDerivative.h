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



//-------------------------------------------------------------
// File    : UnivariateDensityDerivative.h
// Purpose : Header file for UnivariateDensityDerivative.cpp
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : September 17, 2005
//-------------------------------------------------------------
// Direct implementation of the r^{th} kernel density derivative
// estimate based on the Gaussian kernel.
// [HANDLES ONLY UNIVARIATE CASE]
//
// Implementation based on:
//
// V. C. Raykar and R. Duraiswami 'Very fast optimal bandwidth 
// selection for univariate kernel density estimation'
// Technical Report CS-TR-4774, Dept. of Computer 
// Science, University of Maryland, College Park.
// ------------------------------------------------------------
//
// INPUTS [6] 
// ----------------
// NSources		  --> number of sources, N.
// MTargets		  --> number of targets, M.
// pSources		  --> pointer to sources, px(N).
// pTargets       --> pointer to the targets, py(M).
// Bandwidth	  --> the source bandwidth, h.
// Order          --> order of the derivative, r.  
//
// OUTPUTS [1]
// ----------------
// pDensityDerivative --> pointer the the evaluated Density 
//						 Derivative, pD(M).
//-------------------------------------------------------------------

#ifndef UNIVARIATE_DENSITY_DERIVATIVE_H
#define UNIVARIATE_DENSITY_DERIVATIVE_H


class UnivariateDensityDerivative{
	public:
		//constructor 
		UnivariateDensityDerivative(int NSources,
			int MTargets,
			double *pSources,
			double *pTargets,
			double Bandwidth,
		    int Order,
			double *pDensityDerivative);

		//destructor
		~UnivariateDensityDerivative();

		//function to evaluate the Density Derivative
		void Evaluate();

		//function to evaluate the Hermite polynomial.
		double hermite(double x, int r);

	private:
		int N;				//number of sources.
		int M;				//number of targets.
		double *px;			//pointer to sources, (N).
		double *py;		    //pointer to the targets, (M).
		double  h;			//the source bandwidth.
		int r;              //the rth density derivative.
		double *pD;         //pointer the the evaluated Density Derivative, (M).
  
};


#endif