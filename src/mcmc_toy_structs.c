#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "mcmc_toy_structs.h"
#include "mcmc_toy.h"


// these are the functions that go with the structs in structs.h

// construct a state to a point chosen u.a.r. among the (Ngrid_max+1)^n_dimensions points
State* construct_state(int n_dimensions, int Ngrid_max){ 
  State* the_state = (State*)malloc(sizeof(State));
  the_state->n_dimensions = n_dimensions;
  the_state->point = (int*)malloc(n_dimensions*sizeof(int));
  for(int i=0; i<n_dimensions; i++){
    the_state->point[i] = (int)( drand() * (Ngrid_max+1) );
  }
  the_state->prob = F(n_dimensions, the_state->point, Ngrid_max);
  return the_state;
}

void free_state(State* s){
  free(s->point);
  free(s);
}

// Ndim Histogram
Ndim_histogram* construct_histogram_ndim(int Ndim, int Ngrid_max){
  Ndim_histogram* histogram = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  histogram->Ndim = Ndim;
  histogram->Ngrid_max = Ngrid_max;
  histogram->total_weight = 0;
  histogram->weights = construct_ndim_array_of_double(Ndim, Ngrid_max+1, 0.0);
  return histogram;
}

void free_histogram_ndim(Ndim_histogram* h){
  free_ndim_array_of_double(h->weights);
  free(h);
}

// Ndim_array_of_double
Ndim_array_of_double* construct_ndim_array_of_double(int Ndim, int Nsize, double init_value){
  Ndim_array_of_double* array_struct = (Ndim_array_of_double*)malloc(sizeof(Ndim_array_of_double));
  array_struct->Ndim = Ndim;
  array_struct->Nsize = Nsize;
  if(Ndim == 1){
    double* the_array = (double*)malloc(Nsize*sizeof(double));
    for(int i=0; i<Nsize; i++){
      the_array[i] = init_value;
    }
    array_struct->array = the_array; 
  }else{
    Ndim_array_of_double** the_array = (Ndim_array_of_double**)malloc(Nsize*sizeof(Ndim_array_of_double*));
    for(int i=0; i<Nsize; i++){
      the_array[i] = construct_ndim_array_of_double(Ndim-1, Nsize, init_value);
    }
    array_struct->array = the_array;
  }
  return array_struct;
}

void free_ndim_array_of_double(Ndim_array_of_double* A){
  if(A->Ndim == 1){
    free((double*)A->array);
  }else{
    Ndim_array_of_double** a = ((Ndim_array_of_double**)A->array);
    for(int i=0; i<A->Nsize; i++){      
      free_ndim_array_of_double( a[i] );
    }  
    free(a);
  }
  free(A);
}

double set_ndim_array_of_double_with_function(Ndim_array_of_double* array_struct,
					      double outer_value, double_func_ptr function){
  // for probs. of form Prod(i=1,ndim)(f(x_i)), give it a 1-dim function of x (in [0,1]) as the last argument
  int Ndim = array_struct->Ndim;
  int Nsize = array_struct->Nsize;
  //  printf("ZZZ %i %i \n", Ndim, Nsize);
  double sum = 0.0;
  if(Ndim == 1){
    for(int i=0; i<Nsize; i++){
      double x = (double)i/(double)(Nsize-1);
      //    printf("%8i %8.5g  %8.5g \n", i, outer_value, function(x));
      double innermost_value =  outer_value*function(x);
      ((double*)array_struct->array)[i] = innermost_value;
      sum += innermost_value;
      //   printf("%8.5f ", innermost_value);
    }//printf("\n");
  }else{    
    for(int i=0; i<Nsize; i++){   
      double x = (double)i/(double)(Nsize-1);
      Ndim_array_of_double* inner_array_struct = ((Ndim_array_of_double**)array_struct->array)[i];
      sum +=  set_ndim_array_of_double_with_function(inner_array_struct, outer_value*function(x), function);
    }
  }
  return sum;
}

void normalize_ndim_histogram(Ndim_histogram* pdf){
  multiply_ndim_array_of_double_by_scalar(pdf->weights, 1.0/(pdf->total_weight));
  pdf->total_weight = 1.0; // normalized now, so sum is 1.0
}

void multiply_ndim_array_of_double_by_scalar(Ndim_array_of_double* array_struct, double multiplier){
  int Ndim = array_struct->Ndim;
  int Nsize = array_struct->Nsize;
  if(Ndim == 1){
    for(int i=0; i<Nsize; i++){
      ((double*)array_struct->array)[i] *= multiplier;
    }
  }else{    
    for(int i=0; i<Nsize; i++){   
      Ndim_array_of_double* inner_array_struct = ((Ndim_array_of_double**)array_struct->array)[i];
      multiply_ndim_array_of_double_by_scalar(inner_array_struct, multiplier);
    }
  }
}

double sum_ndim_array_of_double(Ndim_array_of_double* array_struct){
  int Ndim = array_struct->Ndim;
  int Nsize = array_struct->Nsize;
  double sum = 0.0;
  if(Ndim == 1){
    double* the_array = (double*)array_struct->array;
    for(int i=0; i<Nsize; i++){
      sum += the_array[i];
    }
    return sum; 
  }else{
    Ndim_array_of_double** the_array = (Ndim_array_of_double**)array_struct->array;
    for(int i=0; i<Nsize; i++){
      sum += sum_ndim_array_of_double(the_array[i]);
    }
  }
  return sum;
}

double sum_abs_difference_ndim_arrays_of_double(Ndim_array_of_double* a1, Ndim_array_of_double* a2){
  int Ndim = a1->Ndim;
  int Nsize = a1->Nsize;
  assert(a2->Ndim == Ndim   &&  a2->Nsize == Nsize);
  double sum_abs_difference = 0.0;
  if(Ndim == 1){
       double* A1 = (double*)a1->array;
    double* A2 = (double*)a2->array;
    for(int i=0; i<Nsize; i++){
         sum_abs_difference += fabs(A1[i] - A2[i]);
    }
    //   return sum_abs_difference; 
  }else{
    Ndim_array_of_double** A1 = (Ndim_array_of_double**)a1->array;
    Ndim_array_of_double** A2 = (Ndim_array_of_double**)a2->array;
    for(int i=0; i<Nsize; i++){
      sum_abs_difference += sum_abs_difference_ndim_arrays_of_double(A1[i], A2[i]);
    }
  }
  return sum_abs_difference;
}

double* get_pointer_to_element(Ndim_array_of_double* A, int* index_array){
  int Ndim = A->Ndim;
  Ndim_array_of_double* B = A;
  for(int j=0; j<Ndim-1; j++){
      B = ((Ndim_array_of_double**)B->array)[index_array[j]];
  }
  // at this point, B should point to 1 dim 'Ndim_array_of_int', whose array field is a regular array of double:
  double* array_of_double = (double*)B->array;
  return array_of_double + index_array[Ndim-1];
}

void print_ndim_array_of_double(Ndim_array_of_double* A){
  int Ndim = A->Ndim;
  int Nsize = A->Nsize;
  //  printf("In print_ndim_array_of_double, Ndim, Nsize: %i %i \n", Ndim, Nsize);
  if(Ndim == 1){
    printf("%4i %4i   ", Ndim, Nsize);
    for(int i=0; i<Nsize; i++){
      printf("%10.8f ", ((double*)A->array)[i]);
    }printf("\n"); 
  }else{
    printf("%4i %4i \n", Ndim, Nsize);
    for(int i=0; i<Nsize; i++){
      print_ndim_array_of_double( ((Ndim_array_of_double**)A->array)[i] );
    }
  }
}
