#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "globals.h"
#include "util.h"
#include "target.h"
// #include "mcmc2body.h"

long pi_evaluation_count;

// target  function definitions

target* set_up_target_from_argv(char** argv, int* i){
  int i_arg = *i;
  int n_dims = atoi(argv[i_arg++]);
  int n_peaks = atoi(argv[i_arg++]);
  target* the_target = (target*)malloc(sizeof(target));
  the_target->n_dims = n_dims;
  the_target->n_peaks = n_peaks;
  the_target->peaks = (peak**)malloc(n_peaks*sizeof(peak*));
  the_target->mean_x = (double*)calloc(n_dims, sizeof(double));
  double sum_weights = 0;
  double* x_wsum = (double*)calloc(n_dims, sizeof(double));
  for(int i_peak = 0; i_peak < n_peaks; i_peak++){
      
    peak* a_peak = (peak*)malloc(sizeof(peak));
    a_peak->n_dims = n_dims;
    a_peak->type = atoi(argv[i_arg++]);
    //    printf("i_peak, %d; type: %d \n", i_peak, a_peak->type);
    double* position = (double*)malloc(n_dims*sizeof(double));
    for(int i = 0; i < n_dims; i++){
      position[i] = atof(argv[i_arg++]);
      //     printf("i_peak: %d ; i, position[i]: %d %10.6g \n", i_peak, i, position[i]);
    }
    a_peak->position = position;
    a_peak->height = atof(argv[i_arg++]);
    a_peak->width = atof(argv[i_arg++]);
    a_peak->shape = atof(argv[i_arg++]);
    //    printf("   h, w, s: %10.6g %10.6g %10.6g \n", a_peak->height, a_peak->width, a_peak->shape);
    the_target->peaks[i_peak] = a_peak;
    double weight = a_peak->height * pow(a_peak->width, n_dims);
    for(int i = 0; i < n_dims; i++){
      x_wsum[i] += weight * position[i];
    }
    sum_weights += weight;
  }
  for(int i = 0; i < n_dims; i++){
    the_target->mean_x[i] = x_wsum[i] / sum_weights;
  }

  //  the_target->pi_evaluation_count = 0;
  *i = i_arg;
  return the_target;
}

target* set_up_target_from_peaks(int n_peaks, peak** the_peaks){
  target* the_target = (target*)malloc(sizeof(target));
  the_target->n_peaks = n_peaks;
  the_target->peaks = (peak**)malloc(n_peaks*sizeof(peak*)); // array of pointers to peaks.
  for(int i=0; i<n_peaks; i++){
    the_target->peaks[i] = the_peaks[i];
  }
  //  the_target->pi_evaluation_count = 0;
  return the_target;
}

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
  //  the_target->pi_evaluation_count = 0;
  return the_target;
}

void print_target_info(const target* const the_target){
  int n_peaks = the_target->n_peaks;
    printf("# n dimensions: %3d,  n peaks: %3d \n", the_target->n_dims, n_peaks);
  for(int ip = 0; ip < n_peaks; ip++){
    print_peak_info(the_target->peaks[ip]);
  }
  printf("# target mean: ");
  for(int id = 0; id < the_target->n_dims; id++){
    printf("%8.5f ", the_target->mean_x[id]);
  }printf("\n");
}
void print_peak_info(const peak* const a_peak){
  int n_dims = a_peak->n_dims;
  printf("# peak position:  ");
  for(int id = 0; id < n_dims; id++){
    printf("%8.5f ", a_peak->position[id]);
  }printf("\n");
  printf("#    width: %8.5f \n", a_peak->width);
  printf("#    height: %8.5f \n", a_peak->height);
  printf("#    shape: %8.5f \n", a_peak->shape);
}

// the density to be sampled:
double pi(const target* const the_target, const double* const x){
  double result = 0.0;
  for(int k = 0; k < the_target->n_peaks; k++){   
    double peak_result;
    peak* a_peak = the_target->peaks[k];
    double rsq = 0;
    for(int i = 0; i < a_peak->n_dims; i++){
      double dx = (x[i] - a_peak->position[i])/a_peak->width;
      rsq += dx*dx;
    }
    if(1){ // gaussian peak
      peak_result = a_peak->height * exp(-0.5*rsq);
    }else{     
      peak_result = a_peak->height * pow((1.0 + rsq), -1.0*a_peak->shape); // -0.5*(n_dims + 1));
    }
    result += peak_result;
  }
  pi_evaluation_count++;
  return result;
}
