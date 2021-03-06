#include <stdlib.h>
#include <iostream>
#include <stdio.h>
// #include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include "plot.h"

using namespace std;


/* expfit.c -- model functions for exponential + background */

struct data {
  size_t n;
  double * x;
  double * y;
};

// Original Function
int expb_f (const gsl_vector * x, void *data, gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *xPos = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  
  double A = gsl_vector_get (x, 0);
  double mu = gsl_vector_get (x, 1);
  double sig = gsl_vector_get (x, 2);
  double y0 = gsl_vector_get (x, 3);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      double t = xPos[i];
      double Yi = A * exp (-pow(t-mu,2)/(2*sig*sig)) + y0;
      //Add penalization factor here to keep peak positive
      if(A>0)
        gsl_vector_set (f, i, (Yi - y[i]));
      else
        gsl_vector_set (f, i, (Yi - y[i])-10*A);  
    }
  return GSL_SUCCESS;
}

// Derivative of function
int expb_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;
  double *xPos = ((struct data *)data)->x;
  double A = gsl_vector_get (x, 0);
  double mu = gsl_vector_get (x, 1);
  double sig = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      double t = xPos[i];
      double e = exp(-pow(t-mu,2)/(2*sig*sig));
      gsl_matrix_set (J, i, 0, e); 
      gsl_matrix_set (J, i, 1, e*A*(t-mu)/(sig*sig));
      gsl_matrix_set (J, i, 2, e*A*pow(t-mu,2)/(sig*sig*sig));
      gsl_matrix_set (J, i, 3, 1);
    }
  return GSL_SUCCESS;
}

// Function and derivative function
int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
  expb_f (x, data, f);
  expb_df (x, data, J);
  return GSL_SUCCESS;
}

bool distfn(std::pair<double, int> i, std::pair<double, int> j) { return i.first<j.first; }

// Take top 5 highest value, remove the one that's farthest from the mean, and then return the average x,y of the remaining 4
void getinit(vector<double> &x, vector<double> &y, double &maxVal, double &x_maxVal){
  std::priority_queue<std::pair<double, int> > q;
  std::vector<std::pair<double, int> > v;
  for (int i = 0; i < int(y.size()); ++i) {
    q.push(std::pair<double, int>(y[i], x[i]));
  }
  int k = 5; // number of indices we need
  double mean = 0;
  for (int i = 0; i < k; ++i) {
    mean = q.top().second + mean;
    v.push_back(make_pair(q.top().first,q.top().second));
    q.pop();
  }
  mean = mean/k;
  for (int i = 0; i < k; ++i) {
    v[i].first = abs(v[i].second-mean);
  }
  v.erase(std::max_element(v.begin(),v.end(),distfn));
  maxVal = 0;
  x_maxVal = 0;
  for (int i = 0; i < k-1; ++i) {
    x_maxVal += v[i].second;
    maxVal += y[v[i].second];
  }
  x_maxVal = x_maxVal/(k-1);
  maxVal = maxVal/(k-1);
}

void gaussfit (vector<double> &x, vector<double> &y, double &A, double &mu, double &sig, double &y0){
  uint N = 250;
  double FIT_RATIO = 1;
  double FIT_OFFSET = 0;

  double maxVal;
  double x_maxVal;
  double sig_init = 5;
  getinit(x, y, maxVal, x_maxVal);
  double x_init[4] = {maxVal, x_maxVal, sig_init, 0.0};
  
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int iter = 0;
  const size_t n = N;

  // Number of variables
  const size_t p = 4;

  struct data d = {n,  &x[0], &y[0]};
  gsl_multifit_function_fdf f;
  gsl_vector_view xOrig = gsl_vector_view_array (x_init, p);
  // const gsl_rng_type * type;
  // gsl_rng * r;

  gsl_rng_env_setup();

  // Random number generator
  // type = gsl_rng_default;
  // r = gsl_rng_alloc (type);

  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  // Set T to the Levenberg–Marquardt derivative solver
  T = gsl_multifit_fdfsolver_lmsder;

  // Allocate new instance of LM derivative solver
  s = gsl_multifit_fdfsolver_alloc (T, n, p);

  // Initialize solver with function f and initial vector x
  gsl_multifit_fdfsolver_set (s, &f, &xOrig.vector);
  
  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
      if (status)
        break;
      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);
  #define FIT(i) gsl_vector_get(s->x, i)
  A=FIT(0);
  mu=FIT(1)*FIT_RATIO-FIT_OFFSET;
  sig=FIT(2)*FIT_RATIO;
  y0=FIT(3);

  // free all the memory associated with the solver s.
  gsl_multifit_fdfsolver_free (s);
  // gsl_rng_free (r);
}
