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

// construct a state
State* construct_state(int n_dimensions, const Target_1dim* targ_1d, int t_index, double* temperatures){ 
  State* the_state = (State*)calloc(1, sizeof(State));
  the_state->n_dimensions = n_dimensions;
  the_state->walker_index = t_index; // initially walker_index is = t_index
  the_state->most_recent_end = -1; // 
  the_state->hottest_since_cold = t_index; // index of hottest since being cold 
  the_state->point = (double*)calloc(n_dimensions, sizeof(double));
  for(int i=0; i<n_dimensions; i++){
    the_state->point = draw_ndim(n_dimensions, targ_1d, temperatures[t_index]); // draw from target distribution
  }
  the_state->prob = (USELOGPROB)? -1.0 : f_ndim(targ_1d, n_dimensions, the_state->point);
  the_state->log_prob = (USELOGPROB)? log_f_ndim(targ_1d, n_dimensions, the_state->point) : 0.0;
  return the_state;
}

void free_state(State* s){
  free(s->point);
  free(s);
}

// Single_T_chain
Single_T_chain* construct_single_T_chain(int chain_number, double T, double rate, Proposal P, State* S, const Binning_spec_set* bins, const char* type){
  Single_T_chain* single_T_chain = (Single_T_chain*)calloc(1, sizeof(Single_T_chain));
  single_T_chain->chain_number = chain_number;
  single_T_chain->current_state = S;
  single_T_chain->sum_x = (double*)calloc(S->n_dimensions, sizeof(double));
  single_T_chain->temperature = T;
  single_T_chain->rate = rate;
  single_T_chain->generation = 0;
  single_T_chain->within_T_steps_taken = 0;
  single_T_chain->n_try = 0;
  single_T_chain->n_accept = 0;
  single_T_chain->n_try_1 = 0;
  single_T_chain->n_try_2 = 0;
  single_T_chain->n_accept_1 = 0;
  single_T_chain->n_accept_2 = 0;
  single_T_chain->proposal = P;
  single_T_chain->dsq_sum = 0.0;
  single_T_chain->n_jumps = 1;
  single_T_chain->type = type;
  single_T_chain->first_tvd_etc_out_done = 0;
 
    single_T_chain->mcmc_out_hist = construct_ndim_histogram(S->n_dimensions, bins->ndim_bins);
   single_T_chain->mcmc_out_1d0_hist = construct_ndim_histogram(1, bins->onedim_bins);
single_T_chain->mcmc_out_1dall_hist = construct_ndim_histogram(1, bins->onedim_bins);

    single_T_chain->mcmc_out_orthants_hist = construct_ndim_histogram(S->n_dimensions, bins->orthants_bins);
    single_T_chain->mcmc_out_reflected_hist = construct_ndim_histogram(S->n_dimensions, bins->positive_bins);
  
    // TVD and K-S D statistic stuff
    single_T_chain->n_old_gees = 0;
    single_T_chain->old_gees = NULL;
    single_T_chain->n_next_block_gees = FIRST_SUMMARY_GEN;
    single_T_chain->next_block_gees = (double*)calloc(FIRST_SUMMARY_GEN, sizeof(double));
  
  return single_T_chain;
}

State* single_T_chain_mcmc_step(Single_T_chain* chain){
  State* the_state = chain->current_state;
  int Ndim = the_state->n_dimensions;
  double inverse_T = 1.0/chain->temperature;
  
  double* x_array = the_state->point;
  double p_of_current_state_T = pow( the_state->prob, inverse_T ); // the_state->prob is the T=1 prob.
  
  Proposal prop = chain->proposal;
  int which_prop; // 1 or 2
  double* prop_x_array = propose(Ndim, x_array, &prop, &which_prop);
  double p_of_proposed_state = (USELOGPROB)? -1.0 : f_ndim(g_targ_1d, Ndim, prop_x_array);
  double log_p_of_proposed_state = (USELOGPROB)? log_f_ndim(g_targ_1d, Ndim, prop_x_array) : 0.0;
  //  printf("log(p), log_p: %g %g %g \n\n", log(p_of_proposed_state), log_p_of_proposed_state, fabs(log(p_of_proposed_state) - log_p_of_proposed_state) );
  /* if(fabs(log(p_of_proposed_state) - log_p_of_proposed_state) > 1e-10){ */
  /*   fprintf(stderr, "log(p), log_p: %g %g %g \n", log(p_of_proposed_state), log_p_of_proposed_state, fabs(log(p_of_proposed_state) - log_p_of_proposed_state) ); */
  /*   //  exit(1); */
  /* } */

  double p_of_proposed_state_T = (USELOGPROB)? -1.0 : pow( p_of_proposed_state, inverse_T );
  double log_p_of_proposed_state_T = (USELOGPROB)? log_p_of_proposed_state*inverse_T : 0.0;
  int n_jumps = 0;
  double dsq = 0;

  // prop ratio is 1 for symmetric proposal
  if(strcmp(chain->type, "mcmc") == 0){ // mcmc step
    int accept = (USELOGPROB)?
      ((log_p_of_proposed_state > the_state->log_prob) || (log(drand()) < inverse_T*(log_p_of_proposed_state - the_state->log_prob))) :
      ((p_of_proposed_state_T > p_of_current_state_T) || (p_of_proposed_state_T > drand() * p_of_current_state_T));
    if(accept){ // accept
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
      the_state->log_prob = log_p_of_proposed_state;
      chain->n_accept++;
      //      fprintf(stderr, "which_prop, naccept 1,2: %i %i %i   %g\n", which_prop, chain->n_accept_1, chain->n_accept_2, chain->temperature);
      if(which_prop == 1){ chain->n_accept_1++; }else{ chain->n_accept_2++; }
    }else{ 
      // reject proposed move. 
      free(prop_x_array);
    }
  }else{ // do iid draw
    double* exact_xs = draw_ndim(Ndim, g_targ_1d, chain->temperature);
    for(int i=0; i<Ndim; i++){
      double d = the_state->point[i] - exact_xs[i];
      dsq += d*d;
    }
	
    free(the_state->point);
    the_state->point = exact_xs;
    the_state->prob = f_ndim(g_targ_1d, Ndim, exact_xs);
    the_state->log_prob = log_f_ndim(g_targ_1d, Ndim, exact_xs);
    chain->n_accept++;
  }  

  chain->dsq_sum += dsq; 
  chain->n_jumps += n_jumps;
  chain->within_T_steps_taken++;
  chain->n_try++; 
  if(which_prop == 1){ chain->n_try_1++; }else{ chain->n_try_2++; }
  chain->current_state = the_state;
  return the_state;
}



void single_T_chain_histogram_current_state(Single_T_chain* chain){
  State* the_state = chain->current_state;
  int n_dim = the_state->n_dimensions;
  if(chain->mcmc_out_hist != NULL){
    add_data_pt_to_ndim_histogram(chain->mcmc_out_hist, n_dim, the_state->point);
    add_data_pt_to_ndim_histogram(chain->mcmc_out_orthants_hist, n_dim, the_state->point);
    double* reflected_point = (double*)calloc(n_dim, sizeof(double));
    for(int i=0; i<n_dim; i++){
      reflected_point[i] = fabs(the_state->point[i]);
    }
    add_data_pt_to_ndim_histogram(chain->mcmc_out_reflected_hist, n_dim, reflected_point);
    free(reflected_point);
    add_data_pt_to_ndim_histogram(chain->mcmc_out_1d0_hist, 1, the_state->point); // zeroth component
    for(int i=0; i<n_dim; i++){
      add_data_pt_to_ndim_histogram(chain->mcmc_out_1dall_hist, 1, the_state->point+i);
    }
  }
}

void single_T_chain_output_tvd(Single_T_chain* chain){

  double* all_gees = (double*)calloc(chain->generation, sizeof(double));
  int n_recent_gees = chain->generation - chain->n_old_gees; // number of states stored in next_block (i.e. after old_gees)
  assert(n_recent_gees == chain->n_next_block_gees || g_n_pi_evaluations >= g_max_pi_evaluations );
  qsort(chain->next_block_gees, n_recent_gees, sizeof(double), cmpfunc);
  if(chain->old_gees == NULL){
    all_gees = copy_array(n_recent_gees, chain->next_block_gees);
  }else{
    all_gees = merge_sorted_arrays(chain->n_old_gees, chain->old_gees, n_recent_gees, chain->next_block_gees);
  }
  if(chain->chain_number == 0){ // only output cold-chain info.
 
    fprintf(g_tvd_vs_gen_fstream, "%8ld %8i %8i %8i  ", g_n_pi_evaluations, chain->generation, chain->n_try, chain->n_accept);

  // tvds cols 5-9
  if(chain->mcmc_out_hist != NULL){
    double tvd = total_variation_distance(g_targprobs_ndim, chain->mcmc_out_hist);
    fprintf(g_tvd_vs_gen_fstream, "%10.7g ", tvd);
    tvd = total_variation_distance(g_targprobs_orthants, chain->mcmc_out_orthants_hist);
    fprintf(g_tvd_vs_gen_fstream, "%10.7g ", tvd);
    tvd = total_variation_distance(g_targprobs_one_orthant, chain->mcmc_out_reflected_hist);
    fprintf(g_tvd_vs_gen_fstream, "%10.7g ", tvd);
    tvd = total_variation_distance(g_targprobs_1dim, chain->mcmc_out_1d0_hist);
    fprintf(g_tvd_vs_gen_fstream, "%10.7g ", tvd);
    tvd = total_variation_distance(g_targprobs_1dim, chain->mcmc_out_1dall_hist);
    fprintf(g_tvd_vs_gen_fstream, "%10.7g  ", tvd);
  }

  double dmusq = 0.0;
  for(int i=0; i<chain->current_state->n_dimensions; i++){ 
    double dmu = chain->sum_x[i]/chain->generation - g_targ_1d->mean; //  muhat - mu
    dmusq += dmu*dmu; // get square of (muhat - mu) vector.
  } 
  double KSD =   // anderson_darling_statistic2(chain->within_T_steps_taken, all_gees, cdf);
    Kolmogorov_smirnov_D_statistic_1_sample(chain->generation, all_gees, cdf);
  //  double ADS = anderson_darling_statistic(chain->within_T_steps_taken, all_gees, cdf);
  double EDFS = edf_statistic(chain->generation, all_gees, cdf);
  // cols 10-14
  fprintf(g_tvd_vs_gen_fstream, "%10.7g %10.7g %10.7g %10.7g %10.7g %8.4g  ",  
	  dmusq, KSD, EDFS, pow(chain->sum_q/chain->generation, 2.0), chain->sum_p/chain->generation, 
	  chain->temperature);
  fprintf(g_tvd_vs_gen_fstream, "\n");
  }
 
  if(chain->old_gees != NULL){ free(chain->old_gees); }
  chain->old_gees = all_gees;
  chain->n_old_gees = chain->generation; // chain->n_old_gees + chain->n_next_block_gees;
  //  printf("chain->n_old_gees: %i \n", chain->n_old_gees);
  free(chain->next_block_gees);
  int next_n_next_block_gees = (int)((SUMMARY_GEN_FACTOR - 1.0)*chain->n_old_gees);
  chain->n_next_block_gees = next_n_next_block_gees;
  chain->next_block_gees = (double*)calloc(next_n_next_block_gees, sizeof(double));
}

void single_T_chain_output_state(Single_T_chain* chain){
  
  fprintf(g_state_vs_gen_fstream, "%2i %2i %2i  ", chain->current_state->walker_index, chain->current_state->most_recent_end, chain->current_state->hottest_since_cold);
  print_array_of_double(g_state_vs_gen_fstream, 
                        chain->current_state->n_dimensions, chain->current_state->point);
  // fprintf(g_state_vs_gen_fstream, "\n");
}

void free_single_T_chain(Single_T_chain* chain){
  free_state(chain->current_state);
  free_ndim_histogram(chain->mcmc_out_hist);
  free(chain); 
}

// Multi_T_chain
Multi_T_chain* construct_multi_T_chain(int n_temperatures, double* temperatures, double* rates, Proposal* proposals, State** states, const Binning_spec_set* bins, const char* type){
  Multi_T_chain* multi_T_chain = (Multi_T_chain*)calloc(1, sizeof(Multi_T_chain));
  multi_T_chain->n_temperatures = n_temperatures;
  multi_T_chain->generation = 0;
  multi_T_chain->count_within_T_updates = 0;
  multi_T_chain->coupled_chains = (Single_T_chain**)calloc(n_temperatures, sizeof(Single_T_chain*));
  multi_T_chain->walker_index_to_T_index = (int*)malloc(n_temperatures*sizeof(int));
  multi_T_chain->temperature_swap_tries_accepts = construct_ndim_array_of_int(2, multi_T_chain->n_temperatures, 0);
  for(int i=0; i<n_temperatures; i++){
    multi_T_chain->coupled_chains[i] = construct_single_T_chain(i, temperatures[i], rates[i], proposals[i], states[i], bins, type);
    multi_T_chain->walker_index_to_T_index[i] = i;
  }
  return multi_T_chain;
}

void multi_T_chain_within_T_mcmc_step(Multi_T_chain* multi_T_chain){ 
  int n_dim = multi_T_chain->coupled_chains[0]->current_state->n_dimensions;
  for(int i=0; i<multi_T_chain->n_temperatures; i++){
    // printf("i: %i \n", i);
    Single_T_chain* chain = multi_T_chain->coupled_chains[i];
    if(drand() < chain->rate){
      single_T_chain_mcmc_step( chain );
    }
    add_arrays(n_dim, chain->sum_x, chain->current_state->point);
    chain->sum_q += g_shortrange(chain->current_state);
    chain->sum_p += chain->current_state->prob; 
    chain->generation++; // generation increments whether a particular chain does a within_T step or not.
    if(i==0 /* || (chain->generation % (i+1) == 0 */){ single_T_chain_histogram_current_state(chain); 
    int j = chain->generation - chain->n_old_gees - 1;
    //  printf("XXXXXXXXXXXXXXXXXXX\n");
    //  printf("# i, j: %i %i  n_next_block: %i %i %i \n", i, j, chain->n_next_block_gees, chain->generation, chain->n_old_gees);
        //    printf("\n");
    assert((j < chain->n_next_block_gees) && (j >= 0));
    chain->next_block_gees[j] = g(chain->current_state);
  }
  }
  multi_T_chain->generation++;

  // output
  for(int i=0; i<multi_T_chain->n_temperatures; i++){    
    if(OUTPUT_SAMPLES){
      Single_T_chain* chain = multi_T_chain->coupled_chains[i];
      printf("%6i  %11.8f %11.8g ", chain->generation, chain->temperature, chain->current_state->prob); 
      print_array_of_double(stdout, n_dim, chain->current_state->point); printf(" ");
      printf("\n");
    }  
  }
  Single_T_chain* chain = multi_T_chain->coupled_chains[0];
  if(chain->generation == chain->n_old_gees + chain->n_next_block_gees){
    multi_T_chain_output_tvd(multi_T_chain);
    multi_T_chain->summary_generation = multi_T_chain->generation;
  }

}

void multi_T_chain_T_swap_mcmc_step(Multi_T_chain* multi_T_chain, int i_c, int i_h){
  assert(i_c < i_h);
  Single_T_chain* cold = multi_T_chain->coupled_chains[i_c];
  Single_T_chain* hot = multi_T_chain->coupled_chains[i_h];
  int cold_walker_index = cold->current_state->walker_index;
  int hot_walker_index = hot->current_state->walker_index;
  double T_c = cold->temperature;
  double T_h = hot->temperature;
  double pr_c =  cold->current_state->prob;
  double pr_h =  hot->current_state->prob;
  double log_pr_c = cold->current_state->log_prob;
  double log_pr_h = hot->current_state->log_prob;
  int ich[2] = {i_c, i_h};
  int ihc[2] = {i_h, i_c};
  int* pch = get_pointer_to_int_element(multi_T_chain->temperature_swap_tries_accepts, ich);
  (*pch)++;
  int accept = (USELOGPROB)?
    ( (log_pr_h > log_pr_c) || (log(drand()) <  (1.0/T_c - 1.0/T_h)*(log_pr_h - log_pr_c) ) ) :
    ( (pr_h > pr_c)  ||  (drand() < pow( (pr_h/pr_c), (1.0/T_c - 1.0/T_h) ) ) ); 
    if(accept){ // accept
    State* swap_state = cold->current_state;
    cold->current_state = hot->current_state;
    hot->current_state = swap_state;
 int* phc = get_pointer_to_int_element(multi_T_chain->temperature_swap_tries_accepts, ihc);
 (*phc)++;
 multi_T_chain->walker_index_to_T_index[cold_walker_index] = i_h;
 multi_T_chain->walker_index_to_T_index[hot_walker_index] = i_c;
 if(i_h == multi_T_chain->n_temperatures-1){
   hot->current_state->most_recent_end = 1;
 }
if(i_c == 0){
   cold->current_state->most_recent_end = 0;
   cold->current_state->hottest_since_cold = 0; 
 }
 if(i_h > hot->current_state->hottest_since_cold){
   hot->current_state->hottest_since_cold = i_h; // this walker has a new hottest temperature is
   // has experienced since being cold.
 }
    } // end of accept 
}

void multi_T_chain_output_tvd(Multi_T_chain* multi_T_chain){
  int n_temperatures = multi_T_chain->n_temperatures;
  double T_hot = multi_T_chain->coupled_chains[n_temperatures-1]->temperature;
  double cold_rate =  multi_T_chain->coupled_chains[0]->rate;
  //  for(int i=0; i<multi_T_chain->n_temperatures; i++){  
    int i=0;
    Single_T_chain* chain = multi_T_chain->coupled_chains[i];
 if(! chain->first_tvd_etc_out_done){
    //  fprintf(g_tvd_vs_gen_fstream, "# n_T T_hot rate n_pi_evals   gen   n_try  n_accept  ");
    fprintf(g_tvd_vs_gen_fstream, "#  nT     T_hot    cold_rate  n_pi_eval  gen    n_try    n_acc    tvd_ndim   tvd_orth   tvd_refl   tvd_1d1   tvd_1dall ");
    fprintf(g_tvd_vs_gen_fstream, "    dmusq      ksd       edfs       sq_avg_q    sq_avg_p      T \n");
    chain->first_tvd_etc_out_done = 1;
  }
    fprintf(g_tvd_vs_gen_fstream, "%5i %9.4f %10.7g ", n_temperatures, T_hot, cold_rate); 
    single_T_chain_output_tvd(chain);
    //  }
    // multi_T_chain->next_summary_generation = (int)(multi_T_chain->summary_generation * SUMMARY_GEN_FACTOR);
}
void multi_T_chain_output_withinT_accept_info(Multi_T_chain* multi_T_chain){
  for(int i=0; i<multi_T_chain->n_temperatures; i++){
    Single_T_chain* chain = multi_T_chain->coupled_chains[i];
    fprintf(g_accept_info_fstream, "# T: %10.6g  %8i %8i   %8i %8i ", chain->temperature, chain->n_try_1, chain->n_accept_1, chain->n_try_2, chain->n_accept_2);
    if (i<multi_T_chain->n_temperatures-1) {
      int x[2] = {i, i+1};
      int y[2] = {i+1, i};
    int* ptries = get_pointer_to_int_element(multi_T_chain->temperature_swap_tries_accepts, x);
    int* paccepts = get_pointer_to_int_element(multi_T_chain->temperature_swap_tries_accepts, y);
    fprintf(g_accept_info_fstream, " %8i %8i \n", *ptries, *paccepts); 
    }else{
        fprintf(g_accept_info_fstream, "\n");
    }
  }
}
void multi_T_chain_output_Tswap_accept_info(Multi_T_chain* multi_T_chain){
  for(int i=0; i<multi_T_chain->n_temperatures; i++){
    fprintf(g_accept_info_fstream, "# %12.5f  ", multi_T_chain->coupled_chains[i]->temperature);
    for(int j=0; j<multi_T_chain->n_temperatures; j++){
      int ij[2] = {j, i};
      int* p = get_pointer_to_int_element(multi_T_chain->temperature_swap_tries_accepts, ij);
      fprintf(g_accept_info_fstream, "%8i ", *p);
    }fprintf(g_accept_info_fstream, "\n");
  }
}

void multi_T_chain_print_state(Multi_T_chain* multi_T_chain){
  fprintf(g_state_vs_gen_fstream, "%6d  ", multi_T_chain->generation);
  for(int i = 0; i < multi_T_chain->n_temperatures; i++){
    single_T_chain_output_state(multi_T_chain->coupled_chains[i]);
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
 
  Ndim_histogram* histogram = (Ndim_histogram*)calloc(1, sizeof(Ndim_histogram));
  histogram->Ndim = n_dim;
  histogram->Ngrid_max = bins->n_bins;
  histogram->total_weight = 0;
  histogram->weights = construct_ndim_array_of_double(n_dim, bins->n_bins, 0.0);
  histogram->bins = bins;
  return histogram;
}

Ndim_histogram* construct_copy_ndim_histogram(const Ndim_histogram* A){
  Ndim_histogram* histogram = (Ndim_histogram*)calloc(1, sizeof(Ndim_histogram));
  histogram->Ndim = A->Ndim;
  histogram->Ngrid_max = A->Ngrid_max;
  histogram->total_weight = A->total_weight;
  histogram->weights = construct_copy_ndim_array_of_double(A->weights);
  return histogram;
}
void add_data_pt_to_ndim_histogram(Ndim_histogram* A, int n_dim, double* xs){
  assert(n_dim == A->Ndim);
  int* i_array = i_array_from_x_array(A->bins, n_dim, xs); //
  //  printf("z %10.7g \n", xs[0]);
  double* bp = get_pointer_to_double_element(A->weights, i_array); 
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
  Ndim_array_of_double* array_struct = (Ndim_array_of_double*)calloc(1, sizeof(Ndim_array_of_double));
  array_struct->Ndim = Ndim;
  array_struct->Nsize = Nsize;
  if(Ndim == 1){
    double* the_array = (double*)calloc(Nsize, sizeof(double));
    for(int i=0; i<Nsize; i++){
      the_array[i] = init_value;
    }
    array_struct->array = the_array; 
  }else{
    Ndim_array_of_double** the_array = (Ndim_array_of_double**)calloc(Nsize, sizeof(Ndim_array_of_double*));
    for(int i=0; i<Nsize; i++){
      the_array[i] = construct_ndim_array_of_double(Ndim-1, Nsize, init_value);
    }
    array_struct->array = the_array;
  }
  return array_struct;
}

Ndim_array_of_double* construct_copy_ndim_array_of_double(Ndim_array_of_double* A){
  Ndim_array_of_double* copy = (Ndim_array_of_double*)calloc(1, sizeof(Ndim_array_of_double));
  copy->Ndim = A->Ndim;
  copy->Nsize = A->Nsize;
  if(A->Ndim == 1){
    double* the_array = (double*)calloc(A->Nsize, sizeof(double));
    for(int i=0; i<A->Nsize; i++){
      the_array[i] = ((double*)A->array)[i];
    }
    copy->array = the_array; 
  }else{
    Ndim_array_of_double** the_array = (Ndim_array_of_double**)calloc(A->Nsize, sizeof(Ndim_array_of_double*));
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
					  double outer_value, const Target_1dim* targ_1d, const Binning_spec* bins){
  Target_peak_1dim* peaks = targ_1d->peaks;
  int Ndim = array_struct->Ndim;
  int Nsize = array_struct->Nsize;
  assert (Nsize == bins->n_bins);
  double* bin_boundaries = bins->bin_edges;
  double sum = 0.0;
  if(Ndim == 1){
    for(int i=0; i<Nsize; i++){
      double bin_prob = 0.0;
      for(int j=0; j<targ_1d->n_modes; j++){
	bin_prob += 
	  peaks[j].weight * (gsl_cdf_gaussian_P(bin_boundaries[i+1]-peaks[j].position, peaks[j].sigma) 
			     - gsl_cdf_gaussian_P(bin_boundaries[i]-peaks[j].position, peaks[j].sigma) );
      };

      double innermost_value =  outer_value*bin_prob;
      ((double*)array_struct->array)[i] = innermost_value;
      sum += innermost_value;
    }
  }else{    
    for(int i=0; i<Nsize; i++){  
      double bin_prob = 0.0; 
      for(int j=0; j<targ_1d->n_modes; j++){
	bin_prob += 
	  peaks[j].weight * (gsl_cdf_gaussian_P(bin_boundaries[i+1]-peaks[j].position, peaks[j].sigma) 
			     - gsl_cdf_gaussian_P(bin_boundaries[i]-peaks[j].position, peaks[j].sigma) );
      };

      Ndim_array_of_double* inner_array_struct = ((Ndim_array_of_double**)array_struct->array)[i];
      sum +=  set_ndim_array_of_double_to_target(inner_array_struct, outer_value*bin_prob, targ_1d, bins);
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
  //  printf("a1, a2->Ndim: %i %i \n", a1->Ndim, a2->Ndim);
  //  printf("a1, a2->Nsize: %i %i \n", a1->Nsize, a2->Nsize);
  assert(a2->Ndim == Ndim   &&  a2->Nsize == Nsize);
  double sum_abs_difference = 0.0;
  if(Ndim == 1){
    double* A1 = (double*)a1->array;
    double* A2 = (double*)a2->array;
    for(int i=0; i<Nsize; i++){
      //   printf("bin probs, targ, data: %i %g %g \n", i, A1[i], A2[i]);
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

double* get_pointer_to_double_element(Ndim_array_of_double* A, int* index_array){
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
  Ndim_array_of_int* array_struct = (Ndim_array_of_int*)calloc(1, sizeof(Ndim_array_of_int));
  array_struct->Ndim = Ndim;
  array_struct->Nsize = Nsize;
  if(Ndim == 1){
    int* the_array = (int*)calloc(Nsize, sizeof(int));
    for(int i=0; i<Nsize; i++){
      the_array[i] = init_value;
    }
    array_struct->array = the_array; 
  }else{
    Ndim_array_of_int** the_array = (Ndim_array_of_int**)calloc(Nsize, sizeof(Ndim_array_of_int*));
    for(int i=0; i<Nsize; i++){
      the_array[i] = construct_ndim_array_of_int(Ndim-1, Nsize, init_value);
    }
    array_struct->array = the_array;
  }
  return array_struct;
}

int* get_pointer_to_int_element(Ndim_array_of_int* A, int* index_array){
  int Ndim = A->Ndim;
  Ndim_array_of_int* B = A;
  //  printf("index array ");
  for(int j=0; j<Ndim-1; j++){
    B = ((Ndim_array_of_int**)B->array)[index_array[j]];
    //  printf(" %i ", index_array[j]);
  }
  // at this point, B should point to 1 dim 'Ndim_array_of_int', whose array field is a regular array of double:
  int* array_of_int = (int*)B->array;
  //  printf(" %i \n", index_array[Ndim-1] ); 
  //  printf(" value: %i \n", array_of_int[index_array[Ndim-1]]);;
  return array_of_int + index_array[Ndim-1];
}

// target 1dim
Target_1dim* construct_target_1dim(int n_peaks, Target_peak_1dim* peaks){
  Target_1dim* targ_1d = (Target_1dim*)calloc(1, sizeof(Target_1dim));
  targ_1d->normalized = 0; // not normalized yet.
  targ_1d->n_modes = n_peaks;
  targ_1d->peaks = peaks; 
  normalize_target_1dim(targ_1d);
  double mean = 0.0;
  double avg_musq = 0.0;
  double sum_sigmasq = 0.0;
  for (int i=0; i<n_peaks; i++){
    double w_i = peaks[i].weight;
    double mu_i = peaks[i].position;
    double sigma_i = peaks[i].sigma;
    mean +=  w_i * mu_i;
    avg_musq += w_i * mu_i * mu_i;
    sum_sigmasq += w_i * sigma_i * sigma_i;
  }
  double variance = (avg_musq - mean*mean) + sum_sigmasq;
  targ_1d->mean = mean;
  targ_1d->variance = variance;
  //  printf("mean, variance: %g %g \n", mean, variance);
  return targ_1d;
} 

void normalize_target_1dim(Target_1dim* targ){
  double w_sum = 0.0;
  for(int i = 0; i<targ->n_modes; i++){
    w_sum += targ->peaks[i].weight;
  }
  for(int i=0; i<targ->n_modes; i++){
    targ->peaks[i].weight /= w_sum;
  }
  targ->normalized = 1; // normalized now.
}

// more general version - can be used with > 2 peaks, with non-equal sigmas and weights
Binning_spec* construct_binning_spec(int n_bins, const Target_1dim* targ_1d, double xlo, double xhi){ 
  double* bin_boundaries = (double*)calloc((n_bins+1), sizeof(double));
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
  Binning_spec* binning_spec = (Binning_spec*)calloc(1, sizeof(Binning_spec));
  binning_spec->n_bins = n_bins;
  binning_spec->x_lo = xlo;
  binning_spec->x_hi = xhi;
  binning_spec->bin_edges = bin_boundaries;
  //   print_array_of_double(n_bins+1, bin_boundaries); printf("\n");
  return binning_spec;
}

Binning_spec* construct_binning_spec_old(int n_bins, const Target_1dim* targ_1d){
  double* bin_boundaries = (double*)calloc((n_bins+1), sizeof(double));
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
  Binning_spec* binning_spec = (Binning_spec*)calloc(1, sizeof(Binning_spec));
  binning_spec->n_bins = n_bins;
  binning_spec->x_lo = -INFINITY;
  binning_spec->x_hi = INFINITY;
  binning_spec->bin_edges = bin_boundaries;
  return binning_spec;
}

void print_binning_spec(const Binning_spec* bins){
  printf("Binning specification: n_bins %i\n", bins->n_bins);
  for(int i=0; i<=bins->n_bins; i++){
    printf("%10.6g \n", bins->bin_edges[i]);
  } printf("\n");
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
  int* i_array = (int*)calloc(Ndim, sizeof(int));
  for(int i=0; i<Ndim; i++){
    //   printf("i, x[i]: %i %g \n", i, x_array[i]);
    i_array[i] = x_to_bin(bins, x_array[i]); // i_from_x(x_array[i]);
    //  printf("AFTER\n");
  }
  return i_array;
} 
