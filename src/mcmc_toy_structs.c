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

// Single_T_chain
/* typedef struct */
/* { */
/*   State current_state; */
/*   double temperature; */
/*   int generation; */
/*   Proposal proposal; // each T can have its own proposal */
/* } Single_T_chain; */
Single_T_chain* construct_single_T_chain(double T, Proposal P, State* S){
 Single_T_chain* single_T_chain = (Single_T_chain*)malloc(sizeof(Single_T_chain));
    single_T_chain->current_state = S;
    single_T_chain->temperature = T;
    single_T_chain->generation = 0;
    single_T_chain->proposal = P;
    return single_T_chain;
}
State* single_T_chain_mcmc_step(Single_T_chain* chain){
  State* the_state = chain->current_state;
  Proposal prop = chain->proposal;
  int Ndim = the_state->n_dimensions;
  int* index_array = the_state->point;
  int* prop_index_array = propose(Ndim, index_array, &prop, is_ball);
  double inverse_T = 1.0/chain->temperature;
  double p_of_proposed_state = pow( F(Ndim, prop_index_array, Ngrid_max), inverse_T );
  double p_of_current_state = pow( the_state->prob, inverse_T );
  int n_jumps = 0;
  int dsq = 0;
  //  printf("ZZZ\n");
// prop ratio is 1 for symmetric proposal
  if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept
    // printf("accept dsq: %i \n", dsq);
    int h0max = Ngrid_max/2;   // count if jumped from < 0.5 to > 0.5 in x and y
    for(int i=0; i<Ndim; i++){
      if( ( the_state->point[i] <= h0max  && prop_index_array[i] > h0max ) 
	  || ( the_state->point[i] > h0max  &&  prop_index_array[i] <= h0max ) ){ 
	n_jumps++;
      }
	int d = the_state->point[i] - prop_index_array[i];
	dsq += d*d;
    }
    // set state to the proposed state:
    free(the_state->point); 
    the_state->point = prop_index_array;
    the_state->prob = p_of_proposed_state;
    n_accept++;
  }else{ 
    // reject proposed move.
    //     printf("reject, dsq: %i \n", dsq);
    n_reject++;
    free(prop_index_array);
    // dsq += 0 
  }
  sum_dsq += (double)dsq;
  //  printf("AAAA \n");
  if(mcmc_out_hist != NULL){
    double* bp = get_pointer_to_element(mcmc_out_hist->weights, the_state->point); 
    (*bp) += 1.0;
    mcmc_out_hist->total_weight += 1.0;
  }
  chain->generation++;
  //  printf("BBBBB \n");
  return the_state;
}
void free_single_T_chain(Single_T_chain* chain){
  free_state(chain->current_state);
  free(chain);
}

// Multi_T_chain
/* typedef struct */
/* { */
/*   int n_temperatures; */
/*   Single_T_chain** coupled_chains; // array of metropolis-coupled chains at different T's */
/* } Multi_T_chain; */
Multi_T_chain* construct_multi_T_chain(int n_temperatures, double* temperatures, Proposal* proposals, State** states){
  Multi_T_chain* multi_T_chain = (Multi_T_chain*)malloc(sizeof(Multi_T_chain));
  multi_T_chain->n_temperatures = n_temperatures;
  multi_T_chain->coupled_chains = (Single_T_chain**)malloc(n_temperatures*sizeof(Single_T_chain*));
  for(int i=0; i<n_temperatures; i++){
    // construct single-T chain:
    /* Single_T_chain* single_T_chain = (Single_T_chain*)malloc(sizeof(Single_T_chain)); */
    /* single_T_chain->current_state = states[i]; */
    /* single_T_chain->temperature = temperatures[i]; */
    /* single_T_chain->generation = 0; */
    /* single_T_chain->proposal = proposals[i]; */
    //
    multi_T_chain->coupled_chains[i] = construct_single_T_chain(temperatures[i], proposals[i], states[i]);
  }
  return multi_T_chain;
}

void multi_T_chain_within_T_mcmc_step(Multi_T_chain* multi_T_chain){ 
  State * state;
  printf("%i  ", multi_T_chain->coupled_chains[0]->generation);
  for(int i=0; i<multi_T_chain->n_temperatures; i++){
    //    printf("in multi T mcmc step. i: %i \n", i);
    state = single_T_chain_mcmc_step( multi_T_chain->coupled_chains[i]);
    print_array_of_int(state->n_dimensions, state->point);
  } printf("\n");
}
void free_multi_T_chain(Multi_T_chain* chain){
  for(int i = 0; i<chain->n_temperatures; i++){
    free_single_T_chain(chain->coupled_chains[i]);
  }
  free(chain);
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

Ndim_histogram* construct_copy_histogram_ndim(Ndim_histogram* A){
  Ndim_histogram* histogram = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  histogram->Ndim = A->Ndim;
  histogram->Ngrid_max = A->Ngrid_max;
  histogram->total_weight = A->total_weight;
  histogram->weights = construct_copy_ndim_array_of_double(A->weights);
  return histogram;
}

void* free_histogram_ndim(Ndim_histogram* h){
  free_ndim_array_of_double(h->weights);
  free(h);
  return NULL;
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

Ndim_array_of_double* construct_copy_ndim_array_of_double(Ndim_array_of_double* A){
  Ndim_array_of_double* copy = (Ndim_array_of_double*)malloc(sizeof(Ndim_array_of_double));
  copy->Ndim = A->Ndim;
  copy->Nsize = A->Nsize;
  if(A->Ndim == 1){
    double* the_array = (double*)malloc(A->Nsize*sizeof(double));
    for(int i=0; i<A->Nsize; i++){
      the_array[i] = ((double*)A->array)[i];
    }
    copy->array = the_array; 
  }else{
    Ndim_array_of_double** the_array = (Ndim_array_of_double**)malloc(A->Nsize*sizeof(Ndim_array_of_double*));
    for(int i=0; i<A->Nsize; i++){
      the_array[i] = construct_copy_ndim_array_of_double( ((Ndim_array_of_double**)A->array)[i] );
    }
    copy->array = the_array;
  }
  return copy;
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
