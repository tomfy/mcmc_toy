#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "util.h"


vector* new_vector(long capacity){
  vector* v = (vector*)malloc(sizeof(vector));
  v->capacity = capacity;
  v->elements = (double*)calloc(capacity, sizeof(double));
  v->count = 0;
  return v;
}

int push(vector* v, double x){
  if(v->count >= v->capacity){ // if full, increase capacity using realloc.
    long new_capacity = (long)(1.2*(v->capacity + 2));
    //    printf("realloc in push. new_capacity: %ld \n", new_capacity);
    double* new_elements = (double*)realloc(v->elements, new_capacity * sizeof(double*));
    if(! new_elements){ printf("realloc in push failed. \n"); free(v->elements); return -1; } // realloc failed.
    v->capacity = new_capacity;
    v->elements = new_elements;
  }
  v->elements[v->count] = x;
  v->count++;
  return 0; // success
}

int pop(vector* v, double* x){
  if(v->count <= 0){
    return -1;
  }else{
    *x = v->elements[v->count-1];
    v->count--;
    return 0;
  }
}

void free_vector(vector* v){
  free(v->elements);
  free(v);
}

// sort-related

// for use with qsort for sorting arrays of doubles.
int compare_doubles(const void * a, const void * b)
{
  double aa = *(double*)a;
  double bb = *(double*)b;
  int result = 0;
  if(aa < bb){
    result = -1;
  }else if(aa > bb){
    result = 1;
  }
  //  printf("%g %g %i \n", aa, bb, result);
  return result;
}





double Dsquared(int n_dims, double* x, double* y){
  double dsq = 0;
  for(int i = 0; i < n_dims; i++){
    dsq += (x[i] - y[i])*(x[i] - y[i]);
  }
  return dsq;
}

double error(int n_dims, double* mean_x, double* sum_x, long n){ 
  double* est_mean_x = (double*)malloc(n_dims*sizeof(double));
  double e = 0;
  for(int i = 0; i < n_dims; i++){
    est_mean_x[i] = sum_x[i] / n;
    //  printf("i: %d est_mean[i]: %8.5g ", i, est_mean_x[i]);
    e += pow(est_mean_x[i] - mean_x[i], 2.0);
  }// printf("\n");
  free(est_mean_x);
  return sqrt(e);
}

double* sum_arrays(int n_dims, double* x1, double* x2){ // get a new array which is the element by element sum of two arrays (of same size)
  double* sum_array = (double*)malloc(n_dims*sizeof(double));
  for(int i=0; i<n_dims; i++){
    sum_array[i] = x1[i] + x2[i];
  }
  return sum_array;
}
double* add_array_to_array(int n_dims, double* x1, double* x2){ // to each element of x1, add the corresponding elem. of x2
  for(int i=0; i<n_dims; i++){
    x1[i] += x2[i];
  }
  return x1;
}
double* mult_array_by_scalar(double A, int n_dims, double* x){ // multiply each element of x by A
 for(int i=0; i<n_dims; i++){
    x[i] *= A;
  }
  return x;
}

double Kolmogorov_Smirnov_D_statistic_2_sample(const int size1, const double* a1, const int size2, const double* a2){
  // calculates the 2-sample KSD statistic for the the two samples stored in arrays a1 and a2
  double d1 = 1.0/(double)size1;
  double d2 = 1.0/(double)size2;
  double c1 = 0.0;
  double c2 = 0.0;
  double D = 0.0;
  int i=0, i1=0, i2=0;
  while(1){
    if(a1[i1] <= a2[i2]){
      double next_c1 = (i1+1)*d1; // c1 + d1;
      if(next_c1 > c2){
	if( next_c1-c2 > D ){ D = next_c1 - c2; }
      }
      c1 = next_c1;
      i1++;
      if(i1 >= size1){ break; }
    }else{
      double next_c2 = (i2+1)*d2;
      if(next_c2 > c1){
	if( next_c2-c1 > D ){ D = next_c2 - c1; }
      }
      c2 = next_c2;
      i2++;
      if(i2 >= size2){ break; }
    }
  }
  return D;
}

