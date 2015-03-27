#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "mcmc_toy_structs.h"
#include "mcmc_toy.h"


// these are the functions that go with the structs in structs.h

// construct a state to a point chosen u.a.r. among the (Ngrid_max+1)^n_dimensions points
State* construct_state(int n_dimensions, int Ngrid_max){ 
  State* the_state = (State*)malloc(sizeof(State));
  the_state->n_dimensions = n_dimensions;
  the_state->ipoint = (int*)malloc(n_dimensions*sizeof(int));
the_state->point = (double*)malloc(n_dimensions*sizeof(double));
  for(int i=0; i<n_dimensions; i++){
    the_state->ipoint[i] = (int)( drand() * (Ngrid_max+1) );
    the_state->point[i] = (double)( drand() * 1.0);
  }
  the_state->prob = F(n_dimensions, the_state->point);
  return the_state;
}

void free_state(State* s){
  free(s->ipoint);
  free(s);
}

// Single_T_chain
Single_T_chain* construct_single_T_chain(double T, Proposal P, State* S){
 Single_T_chain* single_T_chain = (Single_T_chain*)malloc(sizeof(Single_T_chain));
    single_T_chain->current_state = S;
    single_T_chain->temperature = T;
    single_T_chain->generation = 0;
    single_T_chain->proposal = P;
    single_T_chain->mcmc_out_hist = construct_ndim_histogram(S->n_dimensions, Ngrid_max);
    return single_T_chain;
}

State* single_T_chain_mcmc_step(Single_T_chain* chain){
  //  printf("top of mcmc_step.\n");
  State* the_state = chain->current_state;
  Proposal prop = chain->proposal;
  //printf("proposal: %12.5f %12.5f %12.5f\n", prop.W1, prop.W2, prop.p1);
  int Ndim = the_state->n_dimensions;
  double* x_array = the_state->point;
  double* prop_x_array = propose(Ndim, x_array, &prop, is_ball);
  double inverse_T = 1.0/chain->temperature;
  double p_of_proposed_state = pow( F(Ndim, prop_x_array), inverse_T );
  double p_of_current_state = pow( the_state->prob, inverse_T );
  int n_jumps = 0;
  double dsq = 0;
  //   printf("ZZZ\n");
// prop ratio is 1 for symmetric proposal
//  print_array_of_double(Ndim, prop_x_array); printf("  %8.4g  %8.4g  \n", p_of_proposed_state, p_of_current_state);
  if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept
    // printf("accept dsq: %i \n", dsq);
    double h0max = 0.5;   // count if jumped from < 0.5 to > 0.5 in x and y
    for(int i=0; i<Ndim; i++){
      if( ( the_state->point[i] <= h0max  && prop_x_array[i] > h0max ) 
	  || ( the_state->point[i] > h0max  &&  prop_x_array[i] <= h0max ) ){ 
	n_jumps++;
      }
	double d = the_state->point[i] - prop_x_array[i];
	dsq += d*d;
    }
    // set state to the proposed state:
    free(the_state->point); 
    the_state->point = prop_x_array;
    the_state->prob = p_of_proposed_state;
    n_accept++;
  }else{ 
    // reject proposed move.
    //     printf("reject, dsq: %i \n", dsq);
    n_reject++;
    free(prop_x_array);
    // dsq += 0 
  }
  //  printf("%8i ", chain->generation); print_array_of_double(Ndim, the_state->point); printf("\n");
  //  print_array_of_double(Ndim, the_state->point); printf("\n");
  sum_dsq += (double)dsq;
  //    printf("AAAA \n");
  int* index_array = i_array_from_x_array(Ndim, the_state->point); //
  //  print_array_of_double(Ndim, the_state->point); printf("\n");
  //  printf(" "); 
  //   print_array_of_int(Ndim, index_array); printf("\n");
  if(chain->mcmc_out_hist != NULL){
    double* bp = get_pointer_to_element(chain->mcmc_out_hist->weights, index_array); 
    (*bp) += 1.0;
    chain->mcmc_out_hist->total_weight += 1.0;
    //    printf("# *bp, total: %8.5f  %8.5f \n", *bp, chain->mcmc_out_hist->total_weight);
  }
  chain->generation++;
  printf("%8i ", chain->generation); print_array_of_double(Ndim, the_state->point); 
  printf("  "); print_array_of_int(Ndim, index_array); printf("\n");

  if(OUTPUT_TVD_VS_N){
	if(chain->generation%1000 == 0){
	 
	  Ndim_histogram* hist_copy = construct_copy_ndim_histogram(chain->mcmc_out_hist);
	  normalize_ndim_histogram(hist_copy);
	  //printf("ABCD\n");
	  double tvd = total_variation_distance(targp, hist_copy);
	  // printf("AABCD\n");
	  fprintf(tvd_vs_gen_fstream, "%8i %12.6g\n", chain->generation, tvd);
	  free_ndim_histogram(hist_copy);
	}
      }
  // printf("BBBBB \n");
  return the_state;
}
void free_single_T_chain(Single_T_chain* chain){
  free_state(chain->current_state);
  free(chain);
}

// Multi_T_chain
Multi_T_chain* construct_multi_T_chain(int n_temperatures, double* temperatures, Proposal* proposals, State** states){
  Multi_T_chain* multi_T_chain = (Multi_T_chain*)malloc(sizeof(Multi_T_chain));
  multi_T_chain->n_temperatures = n_temperatures;
  multi_T_chain->coupled_chains = (Single_T_chain**)malloc(n_temperatures*sizeof(Single_T_chain*));
  for(int i=0; i<n_temperatures; i++){
    multi_T_chain->coupled_chains[i] = construct_single_T_chain(temperatures[i], proposals[i], states[i]);
  }
  return multi_T_chain;
}

void multi_T_chain_within_T_mcmc_step(Multi_T_chain* multi_T_chain){ 
  State * state;
  if(OUTPUT_SAMPLES){ printf("%i  ", multi_T_chain->coupled_chains[0]->generation); }
  for(int i=0; i<multi_T_chain->n_temperatures; i++){
    state = single_T_chain_mcmc_step( multi_T_chain->coupled_chains[i]);
    if(OUTPUT_SAMPLES){ print_array_of_int(state->n_dimensions, state->ipoint); }
  } if(OUTPUT_SAMPLES){ printf("\n"); }
}
void free_multi_T_chain(Multi_T_chain* chain){
  for(int i = 0; i<chain->n_temperatures; i++){
    free_single_T_chain(chain->coupled_chains[i]);
  }
  free(chain);
}


// Ndim Histogram
Ndim_histogram* construct_ndim_histogram(int Ndim, int n_bins){ // there are bin boundaries x_0, .. x_Ngrid_max, and
  // bins 0 through Ngrid_max+1 ; so Ngrid_max+2 bins!
  Ndim_histogram* histogram = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  histogram->Ndim = Ndim;
  histogram->Ngrid_max = Ngrid_max;
  histogram->total_weight = 0;
  histogram->weights = construct_ndim_array_of_double(Ndim, n_bins, 0.0);
  return histogram;
}

Ndim_histogram* construct_copy_ndim_histogram(Ndim_histogram* A){
  Ndim_histogram* histogram = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  histogram->Ndim = A->Ndim;
  histogram->Ngrid_max = A->Ngrid_max;
  histogram->total_weight = A->total_weight;
  histogram->weights = construct_copy_ndim_array_of_double(A->weights);
  return histogram;
}

void* free_ndim_histogram(Ndim_histogram* h){
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
      double x = Xmin + (double)i/(double)(Nsize-1)*(Xmax-Xmin);
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

double set_ndim_array_of_double_to_target(Ndim_array_of_double* array_struct,
					  double outer_value, Targ_peak_1dim* peaks, Binning_spec* bins){
  /* double peak0_prob = gsl_cdf_gaussian(bin_boundaries[n_bins]-peaks[0].position, peaks[0].sigma)  */
  /*   - gsl_cdf_gaussian(bin_boundaries[0]-peaks[0].position, peaks[0].sigma); */
  /* double peak1_prob = gsl_cdf_gaussian(bin_boundaries[n_bins]-peaks[1].position, peaks[1].sigma)  */
  /*   - gsl_cdf_gaussian(bin_boundaries[1]-peaks[1].position, peaks[1].sigma); */
  // for probs. of form Prod(i=1,ndim)(f(x_i)), give it a 1-dim function of x (in [0,1]) as the last argument

  int Ndim = array_struct->Ndim;
  int Nsize = array_struct->Nsize;
  assert (Nsize == bins->n_bins);
  double* bin_boundaries = bins->bin_edges;
  double sum = 0.0;
  if(Ndim == 1){
    for(int i=0; i<Nsize; i++){
      //double x = Xmin + (double)i/(double)(Nsize-1)*(Xmax-Xmin);
      //    printf("%8i %8.5g  %8.5g \n", i, outer_value, function(x));
      double bin_prob = gsl_cdf_gaussian_P(bin_boundaries[i+1]-peaks[0].position, peaks[0].sigma) 
	- gsl_cdf_gaussian_P(bin_boundaries[i]-peaks[0].position, peaks[0].sigma)
	+ gsl_cdf_gaussian_P(bin_boundaries[i+1]-peaks[1].position, peaks[1].sigma) 
	- gsl_cdf_gaussian_P(bin_boundaries[i]-peaks[1].position, peaks[1].sigma);
      //printf("i, binboundaries[i],[i+1], bin_prob: %i %12.5g  %12.5g  %12.5g \n", i, bin_boundaries[i], bin_boundaries[i+1], bin_prob);

     double innermost_value =  outer_value*bin_prob;
      ((double*)array_struct->array)[i] = innermost_value;
      sum += innermost_value;
      //   printf("%8.5f ", innermost_value);
    }//printf("\n");
  }else{    
    //    printf("Ndim: %i \n", Ndim);
    for(int i=0; i<Nsize; i++){  
  double bin_prob = gsl_cdf_gaussian_P(bin_boundaries[i+1]-peaks[0].position, peaks[0].sigma) 
	- gsl_cdf_gaussian_P(bin_boundaries[i]-peaks[0].position, peaks[0].sigma)
	+ gsl_cdf_gaussian_P(bin_boundaries[i+1]-peaks[1].position, peaks[1].sigma) 
	- gsl_cdf_gaussian_P(bin_boundaries[i]-peaks[1].position, peaks[1].sigma);
  // double bin_prob = gsl_cdf_gaussian_P(bin_boundaries[i+1], peaks[0].sigma) - gsl_cdf_gaussian_P(bin_boundaries[i], peaks[0].sigma)
  //	+ gsl_cdf_gaussian_P(bin_boundaries[i+1], peaks[1].sigma) - gsl_cdf_gaussian_P(bin_boundaries[i], peaks[1].sigma);
  //qprintf("i, binboundaries[i],[i+1], bin_prob: %i %12.5g  %12.5g  %12.5g \n", i, bin_boundaries[i], bin_boundaries[i+1], bin_prob);

      Ndim_array_of_double* inner_array_struct = ((Ndim_array_of_double**)array_struct->array)[i];
      sum +=  set_ndim_array_of_double_to_target(inner_array_struct, outer_value*bin_prob, peaks, bins);
    }
  }
  return sum;
}

/* double set_ndim_array_of_double_with_cdf(Ndim_array_of_double* array_struct, */
/* 					      double outer_value){ */
/*   // for probs. of form Prod(i=1,ndim)(f(x_i)), give it a 1-dim function of x (in [0,1]) as the last argument */
/*   int Ndim = array_struct->Ndim; */
/*   int Nsize = array_struct->Nsize; */
/*   //  printf("ZZZ %i %i \n", Ndim, Nsize); */
/*   double sum = 0.0; */
/*   int Nmiddle = (Ngrid_max-1)/2;  */
/*   if(Ndim == 1){ */
/*     for(int i=0; i<Nmiddle; i++){ */
/*       double x = Xmin + (double)i/(double)(Nsize-1)*(Xmax-Xmin); */
/*       //    printf("%8i %8.5g  %8.5g \n", i, outer_value, function(x)); */
/*       double innermost_value =  outer_value*function(x); */
/*       ((double*)array_struct->array)[i] = innermost_value; */
/*       sum += innermost_value; */
/*       //   printf("%8.5f ", innermost_value); */
/*     }//printf("\n"); */
/*   }else{     */
/*     for(int i=0; i<Nsize; i++){    */
/*       double x = (double)i/(double)(Nsize-1); */
/*       Ndim_array_of_double* inner_array_struct = ((Ndim_array_of_double**)array_struct->array)[i]; */
/*       sum +=  set_ndim_array_of_double_with_function(inner_array_struct, outer_value*function(x), function); */
/*     } */
/*   } */
/*   return sum; */
/* } */


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
  //  printf("a1, a2->Ndim: %i %i \n", a1->Ndim, a2->Ndim);
  // printf("a1, a2->Nsize: %i %i \n", a1->Nsize, a2->Nsize);
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

// Ndim_array_of_int
Ndim_array_of_int* construct_ndim_array_of_int(int Ndim, int Nsize, int init_value){
  Ndim_array_of_int* array_struct = (Ndim_array_of_int*)malloc(sizeof(Ndim_array_of_int));
  array_struct->Ndim = Ndim;
  array_struct->Nsize = Nsize;
  if(Ndim == 1){
    int* the_array = (int*)malloc(Nsize*sizeof(int));
    for(int i=0; i<Nsize; i++){
      the_array[i] = init_value;
    }
    array_struct->array = the_array; 
  }else{
    Ndim_array_of_int** the_array = (Ndim_array_of_int**)malloc(Nsize*sizeof(Ndim_array_of_int*));
    for(int i=0; i<Nsize; i++){
      the_array[i] = construct_ndim_array_of_int(Ndim-1, Nsize, init_value);
    }
    array_struct->array = the_array;
  }
  return array_struct;
}

// target 1dim
void normalize_targ_1dim(int n_modes, Targ_peak_1dim* peaks){
 double w_sum = 0.0;
  for(int i = 0; i<n_modes; i++){
    w_sum += peaks[i].weight;
  }
  for(int i=0; i<n_modes; i++){
    peaks[i].weight /= w_sum;
  }
}

Binning_spec* construct_binning_spec(int n_bins, Targ_peak_1dim* peaks){
  double* bin_boundaries = (double*)malloc((n_bins+1)*sizeof(double));
  printf("n_bins: %i \n", n_bins);

  // assume 2 peaks with same sigma and weight
  double xmid = 0.5*(peaks[0].position + peaks[1].position);
  int Nmiddle = n_bins/2; //n_bins should be even
  double p_this_bin;
  // first peak:
  double position = peaks[0].position;
  double sig = peaks[0].sigma;
  double weight = peaks[0].weight;
  double peak1_prob = gsl_cdf_gaussian_P(xmid,sig) - 0; 
  double pinc = peak1_prob/(n_bins-1);
  double p = pinc;
  bin_boundaries[0] = -INFINITY;
  bin_boundaries[n_bins] = INFINITY;
  for(int i=1; i<=Nmiddle; i++){
    bin_boundaries[i] = gsl_cdf_gaussian_Pinv(p,sig) + position;
    assert(!( i<Nmiddle && bin_boundaries[i]>= xmid) );
    p_this_bin = gsl_cdf_gaussian_P(bin_boundaries[i]-position, sig)
      - gsl_cdf_gaussian_P(bin_boundaries[i-1]-position, sig); 
    //  printf("# %6i  %12.7g  %12.7g  %12.7g  %12.7g\n", i, bin_boundaries[i], p, p_this_bin, f_1dim(bin_boundaries[i]));
 printf("# %6i  %12.7g  %12.7g  %12.7g %12.7g  %12.7g\n", i-1, bin_boundaries[i-1], bin_boundaries[i], 
	   p+peak1_prob, p_this_bin, f_1dim(bin_boundaries[i]) );
    p += 2*pinc;
  }
  bin_boundaries[Nmiddle] = xmid;
  /* p_this_bin =  gsl_cdf_gaussian_P(bin_boundaries[Nmiddle]-position, sig)  */
  /*   - gsl_cdf_gaussian_P(bin_boundaries[Nmiddle-1]-position, sig); */
  /* printf("# %6i  %12.7g  %12.7g  %12.7g  %12.7g\n", Nmiddle-1, bin_boundaries[Nmiddle], p, p_this_bin, f_1dim(bin_boundaries[Nmiddle])); */

  position = peaks[1].position;
  sig = peaks[1].sigma;
  weight = peaks[1].weight;
  double peak2_prob = 1.0 - gsl_cdf_gaussian_P(xmid-position,sig); 
  // printf("peak2 prob: %12.7g \n",  peak2_prob);
  pinc = peak2_prob/(n_bins-1);
  p = 2*pinc;
  for(int i=Nmiddle+1; i<n_bins; i++){
    double p2mid = 
      bin_boundaries[i] = gsl_cdf_gaussian_Pinv(p + gsl_cdf_gaussian_P(xmid-position, sig), sig) + position;
    p_this_bin = gsl_cdf_gaussian_P(bin_boundaries[i]-position, sig);
    p_this_bin -=  gsl_cdf_gaussian_P(bin_boundaries[i-1]-position, sig) ;
    printf("# %6i  %12.7g  %12.7g  %12.7g %12.7g  %12.7g\n", i-1, bin_boundaries[i-1], bin_boundaries[i], 
	   p+peak1_prob, p_this_bin, f_1dim(bin_boundaries[i]) );
    p += 2*pinc;
  }
  double p_top_bin = 1.0 - gsl_cdf_gaussian_P(bin_boundaries[n_bins-1]-position, sig);
  printf("# %6i  %12.7g  %12.7g %12.7g %12.7g  %12.7g\n", n_bins-1, bin_boundaries[n_bins-1], INFINITY, p+peak1_prob-pinc, p_top_bin, 0.0 );
  Binning_spec* binning_spec = (Binning_spec*)malloc(sizeof(Binning_spec));
  binning_spec->n_bins = n_bins;
  binning_spec->x_lo = -INFINITY;
  binning_spec->x_hi = INFINITY;
  binning_spec->bin_edges = bin_boundaries;
  return binning_spec;
}
