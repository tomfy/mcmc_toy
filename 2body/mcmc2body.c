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
#include "chain_architecture.h"
#include "chain_state.h"
#include "mcmc2body.h"

gsl_rng* g_rng;
long pi_evaluation_count; 
int n_thin;

// ************* main: ******************

int main(int argc, char* argv[]){

  // ***** process command line args *****

  // target arguments:
  int i_arg = 1;
  target* the_target = set_up_target_from_argv(argv, &i_arg);

  // chain_architecture arguments:
  int n_levels = atoi(argv[i_arg++]);
  int n_per_level = atoi(argv[i_arg++]); // 1 ('1-body') or 2 ('2_body')
  int symmetry = atoi(argv[i_arg++]); // 0: asymmetrical; 1: symmetrical.
  double Thot = atof(argv[i_arg++]);
  double kernel_width = atof(argv[i_arg++]);
  double proposal_width_factor = atof(argv[i_arg++]);
  double interpolation_power = atof(argv[i_arg++]); // used in symmetrical 2-body

  // run arguments:
  long seed = atoi(argv[i_arg++]);
  double burn_in_arg = atof(argv[i_arg++]);
  long max_n_updates = atoi(argv[i_arg++]); 
  long max_pi_evals = atoi(argv[i_arg++]); // bail out after this many evaluations of target density.
  long burn_in = (burn_in_arg > 0)? (long)burn_in_arg : (long)(fabs(burn_in_arg) * max_pi_evals); // burn-in defined in terms of numbers of pi evaluations.
  long n_runs = atoi(argv[i_arg++]);

  // output arguments:
  n_thin =  atoi(argv[i_arg++]);
  char* output_format = (char*)malloc(32*sizeof(char)); 
  sprintf(output_format, "%s", argv[i_arg++]);
  int verbose = atoi(argv[i_arg++]);

  // ******  
 
  // ***** output the run parameters *****
  printf("# dimensions: %i \n", the_target->n_dims);
  printf("# n_peaks: %i \n", the_target->n_peaks);

  print_target_info(the_target);


  printf("# n_levels: %i \n", n_levels);
  printf("# n_per_level: %i \n", n_per_level);
  printf("# symmetry: %i \n", symmetry);
  printf("# Thot: %10.6g  kernel_width: %10.8g \n", Thot, kernel_width);
  printf("# proposal_width_factor: %10.8g \n", proposal_width_factor);
  printf("# interpolation_power: %10.8g \n", interpolation_power);

  printf("# seed: %ld \n", seed);
  printf("# max updates: %ld \n", max_n_updates);
  printf("# max pi evals: %ld \n", max_pi_evals);
  printf("# burn-in: %ld \n", burn_in);
  printf("# runs: %ld \n", n_runs);

  printf("# n_thin: %d \n", n_thin);
  printf("# output format %s \n", output_format);
  printf("# verbosity: %d \n", verbose);

  double peak_width = min_peak_width(the_target); // ->peaks[0]->width, the_target->peaks[1]->width);
  // int verbose = 0;

  // ***** Set up RNG *****
  gsl_rng_env_setup();
  const gsl_rng_type* rng_type = gsl_rng_default;
  g_rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(g_rng, seed);


  //*********  Set up target distribution *******************
  // target* the_target = set_up_target(n_dims, n_peaks, 0.5*peak_separation, peak_width, height_ratio, shape_param);

  //********** Set up the chain architecture ****************
  double min_prop_width = proposal_width_factor*peak_width;
  double max_prop_width = proposal_width_factor*kernel_width;
  max_prop_width = min_prop_width;
  double min_kernel_width = kernel_width; // peak_width;
  double max_kernel_width = kernel_width;
  chain_architecture* the_arch = set_up_chain_architecture(n_levels, n_per_level, symmetry, Thot, 
                                                           min_prop_width, max_prop_width,
                                                           // proposal_width_factor*peak_width, proposal_width_factor*kernel_width, 
                                                           min_kernel_width, max_kernel_width,
                                                           // peak_width, kernel_width, 
                                                           interpolation_power);
  print_chain_architecture_info(the_arch);

  // loop over runs: 
  for(int i_run = 0; i_run < n_runs; i_run++){
    // ***** Setup initial states of walkers *********
    double init_width = 0.1*peak_width;
    chain_state* the_chain_state = set_up_chain_state(the_target, the_arch, init_width);

    check_state_consistency(the_target, the_arch, the_chain_state);


    
    // *************************************************************
    // ******************* Loop through updates ********************
    // *************************************************************
    long next_summary = 10000;
    pi_evaluation_count = 0;
    long n_updates_done;
    if(n_per_level == 1){ // standard 1-body heating.

      for(long i = 1; i <= burn_in; i++){ // burn-in
        step_1b(the_target, the_arch, the_chain_state);
        if(pi_evaluation_count >= burn_in){ break; }
      }
      reset_chain_state(the_chain_state);
      pi_evaluation_count = 0;

      for(long i = 1; i <= max_n_updates; i++){ // loop through updates       
        step_1b(the_target, the_arch, the_chain_state);
        accumulate_x_sums(the_chain_state);

        if(i % n_thin == 0){
          // check_state_consistency(the_target, the_arch, the_chain_state);
          if(verbose){
            printf("%4ld %4ld ", i, pi_evaluation_count);
            if(strcmp(output_format, "T") == 0){
              print_states_T_order(the_chain_state);
            }else if(strcmp(output_format, "walker") == 0){
              print_states_walker_order(the_chain_state);
            }else{
              print_states_L0_Rtop_only(the_chain_state);
            }
            printf("\n");
          }
        }
        if(pi_evaluation_count >= next_summary){ // get KSD statistic, etc.
          vector* v = the_chain_state->all_coldx0s;
          qsort(v->elements, v->count, sizeof(double), compare_doubles);
          printf("# KSD: %4ld ", pi_evaluation_count);
          double avg_KSD = 0.0;
          for(int iw = 0; iw < 2*the_chain_state->n_levels; iw++){
            vector* vw = the_chain_state->w_coldx0s[iw];
            qsort(vw->elements, vw->count, sizeof(double), compare_doubles);
            double KSD = Kolmogorov_Smirnov_D_statistic_2_sample(v->count, v->elements, vw->count, vw->elements);
            avg_KSD += KSD;
            printf("# %2d %5ld %6.4f  ", iw, vw->count, KSD); 
          }
          printf(" %6.4f\n", avg_KSD/(2*the_chain_state->n_levels));
          next_summary = (long)(next_summary * 1.2);
          if(next_summary > max_pi_evals){ next_summary = max_pi_evals; }
        }
        if(pi_evaluation_count >= max_pi_evals){ n_updates_done = the_chain_state->updates; break ;}
      }
    }else{ // 2-body auxiliary distributions

      for(long i = 1; i <= burn_in; i++){ // do the burnin
        step_2b(the_target, the_arch, the_chain_state);
        if(pi_evaluation_count >= burn_in){ break; }
      }
      reset_chain_state(the_chain_state);
      pi_evaluation_count = 0;

      for(long i = 1; i <= max_n_updates; i++){ // loop through updates

        step_2b(the_target, the_arch, the_chain_state);
        accumulate_x_sums(the_chain_state);

        if(the_arch->symmetry == 1){
          cold_transition_observe_and_count_sym(the_chain_state);
        }else if(the_arch->symmetry == 0){
          cold_transition_observe_and_count_asym(the_chain_state);
        }
    
        if(i % n_thin == 0){
          check_state_consistency(the_target, the_arch, the_chain_state);
          if(verbose){
            printf("%4ld %4ld ", i, pi_evaluation_count);
            if(strcmp(output_format, "T") == 0){
              print_states_T_order(the_chain_state);
            }else if(strcmp(output_format, "walker") == 0){
              print_states_walker_order(the_chain_state);
            }else{
              print_states_L0_Rtop_only(the_chain_state);
            }
            printf("\n");
          }
        }
        if(pi_evaluation_count >= max_pi_evals){ n_updates_done = the_chain_state->updates; break ;}
      } // end of loop through updates
    }

    // *********************************************

    if(i_run == 0){
      printf("# x-update acceptance rates for each walker: \n# ");
      for(int iw = 0; iw < 2*n_levels; iw++){
        printf("  %2i %5.3f", iw, the_chain_state->w_accepts[iw]/(1.0*n_updates_done));
      }printf("\n");

      printf("# X-move acceptance rates:\n");
      for(int it = 0; it < n_levels; it++){
        if(n_per_level == 2){
          printf("#   T-level: %1i P_A: L, R, avg: %5.3f %5.3f \n", it,          
                 the_chain_state->t_Laccepts[it]/(1.0*n_updates_done),
                 the_chain_state->t_Raccepts[it]/(1.0*n_updates_done)
                 ); 
        }else{
          printf("#   T-level: %1i P_A:  %5.3f %5.3f %5.3f \n", it,   
                 the_chain_state->t_Laccepts[it]/(1.0*n_updates_done),
                 the_chain_state->t_Raccepts[it]/(1.0*n_updates_done),       
                 the_chain_state->t_accepts[it]/(2.0*n_updates_done) );
        }
      }

      printf("# T-swap acceptance rates:\n");
      for(int it = 0; it < n_levels-1; it++){
        printf("#   T-levels: %1i-%1i P_A L,R,both: %5.3f %5.3f %5.3f \n", it, it+1, 
               the_chain_state->t_Tswap_Laccepts[it]/(1.0*n_updates_done),
               the_chain_state->t_Tswap_Raccepts[it]/(1.0*n_updates_done),
               the_chain_state->t_Tswap_accepts[it]/(1.0*n_updates_done)); 
      }

      printf("# LR-swap acceptance rates:\n");
      for(int it = 0; it < n_levels; it++){
        printf("#   T-level: %1i P_A: %5.3f \n", it,          
               the_chain_state->t_LRswap_accepts[it]/(1.0*n_updates_done)
               ); 
      }

      for(int it = 0; it < n_levels; it++){  
        printf("# Level: %2i LR Xavg: ", it);
        double* Lxsum = the_chain_state->t_Lxsums[it];
        for(int i_dim = 0; i_dim < the_target->n_dims; i_dim++){
          printf("%8.5f ", Lxsum[i_dim] / n_updates_done);
        }printf("   ");
        double* Rxsum = the_chain_state->t_Rxsums[it];
        for(int i_dim = 0; i_dim < the_target->n_dims; i_dim++){
          printf("%8.5f ", Rxsum[i_dim] / n_updates_done);
        }printf(" \n");
      }

 printf("# pi evaluation count:  %ld \n", pi_evaluation_count);
    }
    // just output the abs error of estimated mean for each level
    /* for(int i = 0; i < the_target->n_dims; i++){ */
    /*   printf("%8.6g ", the_target->mean_x[i]); */
    /* }printf("\n"); */
    //  printf("L abserrors:  ");
    printf("# Abs Errors: ");
    for(int it = 0; it < n_levels; it++){
      printf("%8.6g ", error(the_target->n_dims, the_target->mean_x, the_chain_state->t_Lxsums[it], n_updates_done));
    }printf("  ");
    //   printf("R abserrors:  ");
    for(int it = 0; it < n_levels; it++){
      printf("%8.6g ", error(the_target->n_dims, the_target->mean_x, the_chain_state->t_Rxsums[it], n_updates_done));
    }printf("\n");
    printf("# Pooled Abs Error: %10.8g \n", error(the_target->n_dims, the_target->mean_x, est_mean_pi(the_arch, the_chain_state), 1)); 

  // ********** clean up **************
  free_chain_state(the_chain_state);

  } // end of loop over runs.


  
} // end of main


// ***** functions *****



double* est_mean_pi(chain_architecture* arch, chain_state* state){ // get the 'best' estimates of the mean of pi
  // i.e. using all the nodes which have pi as their stationary distribution.
  int n_dims = state->n_dims;
  int n_levels = arch->n_levels;
  int n_per_level = arch->n_per_level;
  int symmetry = arch->symmetry;
  long updates = state->updates; // so far ...
  double* xavg;
  if(n_per_level == 1){ // standard 1-body pi^(1/T) heating.
    xavg = sum_arrays(n_dims, state->t_Lxsums[0], state->t_Rxsums[0]);
    xavg = mult_array_by_scalar(0.5/updates, n_dims, xavg);
    printf("# 1body: xavg[0]: %10.8g \n", xavg[0]);
  }else if(n_per_level == 2){ // 2-body
    if(symmetry == 0){ // asymmetric
      xavg = (double*)calloc(n_dims, sizeof(double));
      for(int i=0; i<n_levels; i++){
        xavg = add_array_to_array(n_dims, xavg, state->t_Lxsums[i]);
      }
      xavg = mult_array_by_scalar(1.0/(n_levels * updates), n_dims, xavg);
      printf("# 2bodyA: xavg[0]: %10.8g \n", xavg[0]);
    }else if(symmetry == 1){ // symmetric
      xavg = sum_arrays(n_dims, state->t_Lxsums[0], state->t_Rxsums[n_levels-1]);
      xavg = mult_array_by_scalar(0.5/updates, n_dims, xavg);
      printf("# 2bodyB: xavg[0]: %10.8g \n", xavg[0]);
    }
  }
  printf("# estimated mean x[0]: %10.7g \n", xavg[0]);
  return xavg;
}


// ****************  the end  ***********************
