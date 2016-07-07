function [D]=UnivariateDensityDerivative(N,M,X,Y,h,r)

% UnivariateDensityDerivative.
% Direct Evaluation.
% C++ Implementation.
% Loads UnivariateDensityDerivative.dll
%--------------------------------------------------------------------------
% Author     :     Vikas.C.Raykar 
% Date        :     17 September 2005
% Contact    :     vikas@cs.umd.edu
%--------------------------------------------------------------------------
% INPUTS
%--------------------------------------------------------------------------
% N                 --> number of source points.
% M                 --> number of target points.
% X                 --> 1 x N matrix of N source points.
% Y                 --> 1 x M matrix of M target points.
% h                 --> source bandwidth or scale.
% r                 --> derivative order.
%--------------------------------------------------------------------------
% OUTPUTS
%--------------------------------------------------------------------------
% D                --> 1 x M vector of the density derivative
%                          evaluated at  each target point.
%--------------------------------------------------------------------------
% Direct implementation of the r^{th} kernel density derivative
% estimate based on the Gaussian kernel.
% [HANDLES ONLY UNIVARIATE CASE]
%
% Implementation based on:
%
% V. C. Raykar and R. Duraiswami 'Very fast optimal bandwidth 
% selection for univariate kernel density estimation'
% Technical Report CS-TR-4774, Dept. of Computer 
% Science, University of Maryland, College Park.
%--------------------------------------------------------------------------

