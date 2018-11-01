#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>
#include "target.h"

// target  function definitions
target* set_up_target(int n_dims, int n_peaks, double spacing, double width, double height_ratio, double shape_param){
  target* the_target = (target*)malloc(sizeof(target));
  the_target->n_peaks = n_peaks;
  the_target->peaks = (peak**)malloc(n_peaks*sizeof(peak*)); // array of pointers to peaks.
  double peak_height = 1.0;
  for(int i=0; i<n_peaks; i++){
    peak* a_peak = (peak*)malloc(sizeof(peak));
    double* peak_position = (double*)malloc(n_dims*sizeof(double));
    for(int j=0; j<n_dims; j++){
      peak_position[j] = (j == 0)? spacing*(-1 + 2*i) : 0.0; // peaks at -A, A, 3A, ..
    }
    a_peak->n_dims = n_dims;
    a_peak->position = peak_position;
    a_peak->width = width;
    a_peak->height = peak_height;
    a_peak->shape = shape_param;
    the_target->peaks[i] = a_peak; 
    peak_height *= height_ratio;
  }
  the_target->pi_evaluation_count = 0;
  return the_target;
}

void print_target_info(target* the_target){
}

// the density to be sampled:
double pi(target* the_target, double* x){
  double result = 0.0;
  for(int k = 0; k < the_target->n_peaks; k++){   
    double peak_result;
    peak* a_peak = the_target->peaks[k];
    double rsq = 0;
    for(int i = 0; i < a_peak->n_dims; i++){
      double dx = (x[i] - a_peak->position[i])/a_peak->width;
      rsq += dx*dx;
    }
    if(1){
      peak_result = exp(-0.5*rsq);
    }else{     
      peak_result = pow((1.0 + rsq), -1.0*a_peak->shape); // -0.5*(n_dims + 1));
    }
    result += a_peak->height * peak_result;
  }
  the_target->pi_evaluation_count++;
  return result;
}
