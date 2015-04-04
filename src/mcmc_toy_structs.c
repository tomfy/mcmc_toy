#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "mcmc_toy_structs.h"
#include "mcmc_toy.h"


// these are the functions that go with the structs in structs.h

// construct a state to a point chosen u.a.r. among the (Ngrid_max+1)^n_dimensions points
State* construct_state(int n_dimensions, const Target_1dim* targ_1d){ 
  State* the_state = (State*)malloc(sizeof(State));
  the_state->n_dimensions = n_dimensions;
  //  the_state->ipoint = (int*)malloc(n_dimensions*sizeof(int));
  the_state->point = (double*)malloc(n_dimensions*sizeof(double));
  for(int i=0; i<n_dimensions; i++){
    //  the_state->ipoint[i] = (int)( drand() * (Ngrid_max+1) ); 
    the_state->point = draw_ndim(n_dimensions, targ_1d); // draw from target distribution
    // [i] = (double)( drand() * 1.0);
  }
  the_state->prob = F(targ_1d, n_dimensions, the_state->point);
  return the_state;
}

void free_state(State* s){
  free(s->point);
  free(s);
}

// Single_T_chain
Single_T_chain* construct_single_T_chain(double T, Proposal P, State* S, const Binning_spec_set* bins){
  Single_T_chain* single_T_chain = (Single_T_chain*)malloc(sizeof(Single_T_chain));
 
  single_T_chain->current_state = S;
  single_T_chain->temperature = T;
  single_T_chain->generation = 0;
  single_T_chain->proposal = P;
  single_T_chain->dsq_sum = 0.0;
  single_T_chain->n_jumps = 1;
  //   single_T_chain->bins = bins;
  if(DO_TVD){
    single_T_chain->mcmc_out_hist = construct_ndim_histogram(S->n_dimensions, bins->ndim_bins);
    //  single_T_chain->mcmc_out_hist->bins = bins;
    single_T_chain->mcmc_out_1d_hists = (Ndim_histogram**)malloc(S->n_dimensions*sizeof(Ndim_histogram*));
    for(int i=0; i<S->n_dimensions; i++){
      single_T_chain->mcmc_out_1d_hists[i] = construct_ndim_histogram(1, bins->onedim_bins);
    }
    single_T_chain->mcmc_out_orthants_hist = construct_ndim_histogram(S->n_dimensions, bins->orthants_bins);
    single_T_chain->mcmc_out_reflected_hist = construct_ndim_histogram(S->n_dimensions, bins->positive_bins);
    if(DO_EXACT){
      single_T_chain->exact_draw_hist = construct_ndim_histogram(S->n_dimensions, bins->ndim_bins);
   single_T_chain->exact_draw_1d_hists = (Ndim_histogram**)malloc(S->n_dimensions*sizeof(Ndim_histogram*));
    for(int i=0; i<S->n_dimensions; i++){
      single_T_chain->exact_draw_1d_hists[i] = construct_ndim_histogram(1, bins->onedim_bins);
    }
    single_T_chain->exact_draw_orthants_hist = construct_ndim_histogram(S->n_dimensions, bins->orthants_bins);
    single_T_chain->exact_draw_reflected_hist = construct_ndim_histogram(S->n_dimensions, bins->positive_bins);
      //   single_T_chain->exact_draw_hist->bins = bins;
    }
  }
  return single_T_chain;
}

State* single_T_chain_mcmc_step(Single_T_chain* chain){
  State* the_state = chain->current_state;
  int Ndim = the_state->n_dimensions;
  double inverse_T = 1.0/chain->temperature;
  
  double* x_array = the_state->point;
  double p_of_current_state = pow( the_state->prob, inverse_T );

  Proposal prop = chain->proposal;
  double* prop_x_array = propose(Ndim, x_array, &prop);
  double p_of_proposed_state = pow( F(g_targ_1d, Ndim, prop_x_array), inverse_T );

  int n_jumps = 0;
  double dsq = 0;
  // prop ratio is 1 for symmetric proposal
  if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept
    double h0max = 0.0;   // count if jumped from < 0.5 to > 0.5 in x and y
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
    chain->n_accept++;
  }else{ 
    // reject proposed move. // dsq += 0 , n_jumps += 0
    chain->n_reject++;
    free(prop_x_array);
  }
  chain->dsq_sum += dsq;
  chain->n_jumps += n_jumps;
  chain->generation++;
  chain->current_state = the_state;
  return the_state;
}

void single_T_chain_histogram_current_state(Single_T_chain* chain){
  State* the_state = chain->current_state;
  int n_dim = the_state->n_dimensions;
  if(chain->mcmc_out_hist != NULL){
    add_data_pt_to_ndim_histogram(chain->mcmc_out_hist, n_dim, the_state->point);
    add_data_pt_to_ndim_histogram(chain->mcmc_out_orthants_hist, n_dim, the_state->point);
    double* reflected_point = (double*)malloc(n_dim*sizeof(double));
    for(int i=0; i<n_dim; i++){
      reflected_point[i] = fabs(the_state->point[i]);
    }
    add_data_pt_to_ndim_histogram(chain->mcmc_out_reflected_hist, n_dim, reflected_point);
    free(reflected_point);
    for(int i=0; i<n_dim; i++){
      add_data_pt_to_ndim_histogram(chain->mcmc_out_1d_hists[i], 1, the_state->point+i);
    }
  }
}

void single_T_chain_histogram_exact_draw(Single_T_chain* chain, double* exact_draw){
  State* the_state = chain->current_state;
  int n_dim = the_state->n_dimensions;
  if(chain->exact_draw_hist != NULL){
    add_data_pt_to_ndim_histogram(chain->exact_draw_hist, n_dim, exact_draw);
    add_data_pt_to_ndim_histogram(chain->exact_draw_orthants_hist, n_dim, exact_draw);
    double* reflected_point = (double*)malloc(n_dim*sizeof(double));
    for(int i=0; i<n_dim; i++){
      reflected_point[i] = fabs(exact_draw[i]);
    }
    add_data_pt_to_ndim_histogram(chain->exact_draw_reflected_hist, n_dim, reflected_point);
    free(reflected_point);
    for(int i=0; i<n_dim; i++){
      add_data_pt_to_ndim_histogram(chain->exact_draw_1d_hists[i], 1, exact_draw+i);
    }
  }
}

void single_T_chain_output_tvd(Single_T_chain* chain){
  fprintf(g_tvd_vs_gen_fstream, "%8i  ", chain->generation); 
  if(chain->mcmc_out_hist != NULL){
    double tvd = total_variation_distance(g_targprobs_ndim, chain->mcmc_out_hist);
    fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
    tvd = total_variation_distance(g_targprobs_orthants, chain->mcmc_out_orthants_hist);
    fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
    tvd = total_variation_distance(g_targprobs_one_orthant, chain->mcmc_out_reflected_hist);
    fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
    for(int i=0; i<chain->current_state->n_dimensions; i++){
      //   printf("g_targprobs_1dim, Ndim: %i \n", g_targprobs_1dim->Ndim);
      tvd = total_variation_distance(g_targprobs_1dim, chain->mcmc_out_1d_hists[i]);
      fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
    }

    /* Ndim_histogram* hist_copy = construct_copy_ndim_histogram(chain->mcmc_out_hist); */
    /* normalize_ndim_histogram(hist_copy); */
    /* tvd = total_variation_distance(g_targp, hist_copy); */
    /* fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd); */
    /* free_ndim_histogram(hist_copy); */
  }
  if(chain->exact_draw_hist != NULL){
    double tvd = total_variation_distance(g_targprobs_ndim, chain->exact_draw_hist);
    fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
    tvd = total_variation_distance(g_targprobs_orthants, chain->exact_draw_orthants_hist);
    fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
    tvd = total_variation_distance(g_targprobs_one_orthant, chain->exact_draw_reflected_hist);
    fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
    for(int i=0; i<chain->current_state->n_dimensions; i++){
      //   printf("g_targprobs_1dim, Ndim: %i \n", g_targprobs_1dim->Ndim);
      tvd = total_variation_distance(g_targprobs_1dim, chain->exact_draw_1d_hists[i]);
      fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
    }
    /* Ndim_histogram* hist_copy = construct_copy_ndim_histogram(chain->exact_draw_hist); */
    /* normalize_ndim_histogram(hist_copy); */
    /* double tvd = total_variation_distance(g_targp, hist_copy); */
    /* fprintf(g_tvd_vs_gen_fstream, "%12.6g  ", tvd); */
    /* free_ndim_histogram(hist_copy); */
  }
  fprintf(g_tvd_vs_gen_fstream, "\n");
}



void free_single_T_chain(Single_T_chain* chain){
  free_state(chain->current_state);
  free_ndim_histogram(chain->mcmc_out_hist);
  free_ndim_histogram(chain->exact_draw_hist);
  free(chain);
  
}

// Multi_T_chain
Multi_T_chain* construct_multi_T_chain(int n_temperatures, double* temperatures, Proposal* proposals, State** states, const Binning_spec_set* bins){
  Multi_T_chain* multi_T_chain = (Multi_T_chain*)malloc(sizeof(Multi_T_chain));
  multi_T_chain->n_temperatures = n_temperatures;
  multi_T_chain->coupled_chains = (Single_T_chain**)malloc(n_temperatures*sizeof(Single_T_chain*));
  for(int i=0; i<n_temperatures; i++){
    printf("i: %i;  before construct single T chain. \n", i);
    multi_T_chain->coupled_chains[i] = construct_single_T_chain(temperatures[i], proposals[i], states[i], bins);
  }
  return multi_T_chain;
}

void multi_T_chain_within_T_mcmc_step(Multi_T_chain* multi_T_chain){ 
  int n_dim = multi_T_chain->coupled_chains[0]->current_state->n_dimensions;
  for(int i=0; i<multi_T_chain->n_temperatures; i++){
    Single_T_chain* chain = multi_T_chain->coupled_chains[i];
    single_T_chain_mcmc_step( chain );
    double* exact_xs = NULL;
    if(DO_EXACT){
      exact_xs = draw_ndim(n_dim, g_targ_1d);
    }
    if(OUTPUT_SAMPLES){
      printf("%6i  %8.5f %12.8g ", chain->generation, chain->temperature, chain->current_state->prob); 
      print_array_of_double(n_dim, chain->current_state->point); printf(" ");
      print_array_of_int(n_dim, i_array_from_x_array(chain->mcmc_out_hist->bins, n_dim, chain->current_state->point));
      if(DO_EXACT){ print_array_of_double(n_dim, exact_xs); 
	print_array_of_int(n_dim, i_array_from_x_array(chain->exact_draw_hist->bins, n_dim, exact_xs));
      }

      printf("\n");
    }
   
    if(DO_TVD){
      single_T_chain_histogram_current_state(multi_T_chain->coupled_chains[i]);
      if(DO_EXACT){ 
	single_T_chain_histogram_exact_draw(multi_T_chain->coupled_chains[i], exact_xs);
	//add_data_pt_to_ndim_histogram(chain->exact_draw_hist, n_dim,  exact_xs); }
	if(chain->generation%TVD_EVERY == 0){	 
	  single_T_chain_output_tvd(chain);
	}
      }
    }
    free(exact_xs);
  } 
}

void free_multi_T_chain(Multi_T_chain* chain){
  for(int i = 0; i<chain->n_temperatures; i++){
    free_single_T_chain(chain->coupled_chains[i]);
  }
  free(chain);
}

// Ndim Histogram
Ndim_histogram* construct_ndim_histogram
// (int Ndim, int n_bins){ // there are bin boundaries x_0, .. x_Ngrid_max, and
(int n_dim, const Binning_spec* bins){
 
  Ndim_histogram* histogram = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  histogram->Ndim = n_dim;
  histogram->Ngrid_max = bins->n_bins;
  histogram->total_weight = 0;
  histogram->weights = construct_ndim_array_of_double(n_dim, bins->n_bins, 0.0);
  histogram->bins = bins;
  return histogram;
}

Ndim_histogram* construct_copy_ndim_histogram(const Ndim_histogram* A){
  Ndim_histogram* histogram = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  histogram->Ndim = A->Ndim;
  histogram->Ngrid_max = A->Ngrid_max;
  histogram->total_weight = A->total_weight;
  histogram->weights = construct_copy_ndim_array_of_double(A->weights);
  return histogram;
}
void add_data_pt_to_ndim_histogram(Ndim_histogram* A, int n_dim, double* xs){
  assert(n_dim == A->Ndim);
  int* i_array = i_array_from_x_array(A->bins, n_dim, xs); //
  double* bp = get_pointer_to_element(A->weights, i_array); 
  (*bp) += 1.0;
  A->total_weight += 1.0;
  free(i_array);
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

double set_ndim_array_of_double_to_target(Ndim_array_of_double* array_struct,
					  double outer_value, Target_peak_1dim* peaks, const Binning_spec* bins){
 
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

double sum_abs_difference_ndim_arrays_of_double(const Ndim_array_of_double* a1, const Ndim_array_of_double* a2){
  int Ndim = a1->Ndim;
  int Nsize = a1->Nsize;
    printf("a1, a2->Ndim: %i %i \n", a1->Ndim, a2->Ndim);
   printf("a1, a2->Nsize: %i %i \n", a1->Nsize, a2->Nsize);
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

void print_ndim_array_of_double(const Ndim_array_of_double* A){
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
void normalize_targ_1dim(Target_1dim* targ){
  double w_sum = 0.0;
  for(int i = 0; i<targ->n_modes; i++){
    w_sum += targ->peaks[i].weight;
  }
  for(int i=0; i<targ->n_modes; i++){
    targ->peaks[i].weight /= w_sum;
  }
}

// more general version - can be used with > 2 peaks, with non-equal sigmas and weights
Binning_spec* construct_binning_spec(int n_bins, const Target_1dim* targ_1d, double xlo, double xhi){ 
  double* bin_boundaries = (double*)malloc((n_bins+1)*sizeof(double));
  int n_modes = targ_1d->n_modes;
  Target_peak_1dim* peaks = targ_1d->peaks;
  double integral = integral_f_1dim(targ_1d, xlo, xhi);
  double p_bin = integral/(double)n_bins;
  //  printf("integral: %g, p_bin: %g \n", integral, p_bin);
  bin_boundaries[0] = xlo;
  for(int i=1; i<n_bins; i++){   
    bin_boundaries[i] = find_bin_upper_edge( targ_1d, xlo, i*p_bin);
  }
  bin_boundaries[n_bins] = xhi;
  Binning_spec* binning_spec = (Binning_spec*)malloc(sizeof(Binning_spec));
  binning_spec->n_bins = n_bins;
  binning_spec->x_lo = xlo;
  binning_spec->x_hi = xhi;
  binning_spec->bin_edges = bin_boundaries;
  //  print_array_of_double(n_bins+1, bin_boundaries); printf("\n");
  return binning_spec;
}

Binning_spec* construct_binning_spec_old(int n_bins, const Target_1dim* targ_1d){
  double* bin_boundaries = (double*)malloc((n_bins+1)*sizeof(double));
  // printf("n_bins: %i \n", n_bins);
  int n_modes = targ_1d->n_modes;
  Target_peak_1dim* peaks = targ_1d->peaks;
  // assume 2 peaks with same sigma and weight
  double xmid = 0.5*(peaks[0].position + peaks[1].position);
  //  printf("xmid: %12.5g\n", xmid);
  int Nmiddle = n_bins/2; //n_bins should be even
  double p_this_bin;
  //  printf("peak1,2 probs: %12.5g  %12.5g \n", peak1_prob, peak2_prob); 
  // first peak:
  double position = peaks[0].position;
  double sig = peaks[0].sigma;
  double weight = peaks[0].weight;
  double peak1_prob = gsl_cdf_gaussian_P(xmid-position,sig) - 0; 
  double pinc = peak1_prob/(n_bins/2); // (n_bins-1);
  double p = pinc;
  bin_boundaries[0] = -INFINITY;
  bin_boundaries[n_bins] = INFINITY;
  for(int i=1; i<=Nmiddle; i++){
    bin_boundaries[i] = gsl_cdf_gaussian_Pinv(p,sig) + position;
    assert(!( i<Nmiddle && bin_boundaries[i]>= xmid) );
    p_this_bin = gsl_cdf_gaussian_P(bin_boundaries[i]-position, sig)
      - gsl_cdf_gaussian_P(bin_boundaries[i-1]-position, sig); 
    //  printf("# %6i  %12.7g  %12.7g  %12.7g  %12.7g\n", i, bin_boundaries[i], p, p_this_bin, f_1dim(bin_boundaries[i]));
    printf("# %6i  %12.7g  %12.7g  %12.7g  %12.7g\n", i-1, bin_boundaries[i-1], bin_boundaries[i], 
	    p_this_bin, f_1dim(targ_1d, bin_boundaries[i]) );
    p += pinc; // 2*pinc;
  }
  bin_boundaries[Nmiddle] = xmid;
  /* p_this_bin =  gsl_cdf_gaussian_P(bin_boundaries[Nmiddle]-position, sig)  */
  /*   - gsl_cdf_gaussian_P(bin_boundaries[Nmiddle-1]-position, sig); */
  /* printf("# %6i  %12.7g  %12.7g  %12.7g  %12.7g\n", Nmiddle-1, bin_boundaries[Nmiddle], p, p_this_bin, f_1dim(bin_boundaries[Nmiddle])); */

  position = peaks[1].position;
  sig = peaks[1].sigma;
  weight = peaks[1].weight;
  double peak2_prob = 1.0 - gsl_cdf_gaussian_P(xmid-position,sig); 
  // printf("peak1,2 probs: %12.5g  %12.5g \n", peak1_prob, peak2_prob);
  // printf("peak2 prob: %12.7g \n",  peak2_prob);
  pinc = peak2_prob/(n_bins/2); // (n_bins-1);
  p = pinc; // 2*pinc;
  for(int i=Nmiddle+1; i<n_bins; i++){
    double p2mid = 
      bin_boundaries[i] = gsl_cdf_gaussian_Pinv(p + gsl_cdf_gaussian_P(xmid-position, sig), sig) + position;
    p_this_bin = gsl_cdf_gaussian_P(bin_boundaries[i]-position, sig);
    p_this_bin -=  gsl_cdf_gaussian_P(bin_boundaries[i-1]-position, sig) ;
    printf("# %6i  %12.7g  %12.7g  %12.7g %12.7g  %12.7g\n", i-1, bin_boundaries[i-1], bin_boundaries[i], 
	   p+peak1_prob, p_this_bin, f_1dim(targ_1d, bin_boundaries[i]) );
    p += pinc; // 2*pinc;
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

int x_to_bin(const Binning_spec* bin_spec, double x){ // 1-dim
  int n_bins = bin_spec->n_bins;
  double* edges = bin_spec->bin_edges;
  int min_bin = 0;
  int max_bin = n_bins-1;
  int try_bin = max_bin/2;
  while(max_bin > min_bin){
    if(x < edges[try_bin+1]){
      max_bin = try_bin;
    }else{
      min_bin = try_bin+1;
    }
    try_bin = (min_bin + max_bin)/2;
  }  
  return min_bin;
}

int* i_array_from_x_array(const Binning_spec* bins, int Ndim, double* x_array){
  //int* i_array_from_x_array(int Ndim, double* x_array){
  int* i_array = (int*)malloc(Ndim*sizeof(int));
  for(int i=0; i<Ndim; i++){
    //   printf("i, x[i]: %i %g \n", i, x_array[i]);
    i_array[i] = x_to_bin(bins, x_array[i]); // i_from_x(x_array[i]);
    //  printf("AFTER\n");
  }
  return i_array;
}

// obsolete/unused stuff ...

/* double set_ndim_array_of_double_with_function(Ndim_array_of_double* array_struct, */
/* 					      double outer_value, double_func_ptr function){ */
/*   // for probs. of form Prod(i=1,ndim)(f(x_i)), give it a 1-dim function of x (in [0,1]) as the last argument */
/*   int Ndim = array_struct->Ndim; */
/*   int Nsize = array_struct->Nsize; */
/*   //  printf("ZZZ %i %i \n", Ndim, Nsize); */
/*   double sum = 0.0; */
/*   if(Ndim == 1){ */
/*     for(int i=0; i<Nsize; i++){ */
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
