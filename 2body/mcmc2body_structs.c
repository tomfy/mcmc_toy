#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "mcmc2body.h"
#include "mcmc2body_structs.h"


// function definitions:

void old_print_states_walker_order(int n_Ts, int n_d, double** xs, int* Ts){
  int n_w = 2*n_Ts;
  for(int iw=0; iw < n_w; iw++){
    int it = Ts[iw];
    printf("%2i %2i  ", iw, it);
    for(int j=0; j < n_d; j++){
      printf("%5.3f ", xs[iw][j]);
    }printf("   ");
  } // printf("\n");
}

void print_states_T_order(chain_state* state){
  
  for(int it=0; it < state->n_Ts; it++){
    int iw = state->t_Lws[it];
    printf("%2i %2i  ", iw, it);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
  for(int it=0; it < state->n_Ts; it++){
    int iw = state->t_Rws[it];
    printf("%2i %2i  ", iw, it);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
  // printf("\n");
}


double pi(peaks* the_peaks, double* x){
  double result = 0.0;
  for(int k = 0; k < the_peaks->n_peaks; k++){   
    double peak_result;
    peak* a_peak = the_peaks->peaks[k];
    if(1){
      peak_result = 1.0;
      for(int i = 0; i < a_peak->n_dims; i++){
        peak_result *= exp(-0.5*pow( (x[i] - a_peak->position[i])/a_peak->width, 2) );
      }
    }else{ 
      double rsq = 0;
      for(int i = 0; i < a_peak->n_dims; i++){
        //   printf("peak positioni: %g\n", a_peak->position[i]);
        rsq += pow( (x[i] - a_peak->position[i])/a_peak->width, 2 );
      }
      //     printf("new rsq: %g \n", rsq);
      peak_result = pow((1.0 + rsq), -3); // -0.5*(n_dims + 1));
      //    printf("#  rsq: %12.6g  peak_result: %12.6g  %12.6g\n", rsq, peak_result, pow(rsq, -6) );
    }
    result += a_peak->height * peak_result;
  }
  // printf("#  result: %12.6g \n", result);
  return result;
}


void update_x(int n_dims, peaks* the_peaks,  int iw, chain_state* state, double Tinverse){
  //  printf("top of update_x\n");
  double pi_x = state->w_pis[iw]; // pi(n_dims, n_peaks, peak_positions, heights, sigmas, x);
  //  printf("after pi_x = pis[iw]. n_dims: %i \n", n_dims);

  double* prop_x = (double*)malloc(n_dims*sizeof(double));
  //  printf("after alloc prop_x \n"); fflush(stdout);
  for(int i = 0; i < n_dims; i++){ // propose a new position 
    prop_x[i] = state->w_xs[iw][i] + prop_w * (gsl_rng_uniform(g_rng) - 0.5);
  }

  double pi_prop_x = pi(the_peaks, prop_x);
  //  printf("pi_x: %8.5g  %8.5g   %8.5g\n", pi_x, pi_prop_x, Tinverse);
  if( (pi_prop_x > pi_x)  ||  (gsl_rng_uniform(g_rng)  <  pow(pi_prop_x/pi_x, Tinverse)) ){ // Accept proposal
    state->w_pis[iw] = pi_prop_x;
    // return prop_x;
    // double* old_x = xs[iw];
    free(state->w_xs[iw]);  
    state->w_xs[iw] = prop_x;
    // free(old_x);
  }else{
    free(prop_x);
  }
  //  printf("bottom of update_x\n");
}


void T_swap(chain_state* state, double* inverse_Temperatures){

  for(int i=1; i<state->n_Ts; i++){ // T-swapping among Ls
    double Tc_inverse = // 1.0/Temperatures[i-1]; // 1/T for the colder of two walkers
      inverse_Temperatures[i-1];
    double Th_inverse = // 1.0/Temperatures[i]; // 1/T for the hotter of the two walkers
      inverse_Temperatures[i];
    {
      // first the 'left' set of walkers:
      int iw_hot = state->t_Lws[i]; // index of walker with ith temperature
      int iw_cold = state->t_Lws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
      double pi_cold = state->w_pis[iw_cold];
      double pi_hot = state->w_pis[iw_hot];
      if( 
         (pi_hot > pi_cold) 
         || 
         (gsl_rng_uniform(g_rng)  <  pow(pi_hot/pi_cold, Tc_inverse - Th_inverse)) 
          ){ // Accept proposed T-swap
        int tmp = state->t_Lws[i]; 
        state->t_Lws[i] = state->t_Lws[i-1];
        state->t_Lws[i-1] = tmp;

        tmp = state->w_ts[iw_hot];
        state->w_ts[iw_hot] = state->w_ts[iw_cold];
        state->w_ts[iw_cold] = tmp;
      }
    }
    {
      // then the 'right' set of walkers:
      int iw_hot = state->t_Rws[i]; // index of walker with ith temperature
      int iw_cold = state->t_Rws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
      double pi_cold = state->w_pis[iw_cold];
      double pi_hot = state->w_pis[iw_hot];
      if( 
         (pi_hot > pi_cold) 
         || 
         (gsl_rng_uniform(g_rng)  <  pow(pi_hot/pi_cold, Tc_inverse - Th_inverse)) 
          ){ // Accept proposed T-swap
        int tmp = state->t_Rws[i]; 
        state->t_Rws[i] = state->t_Rws[i-1];
        state->t_Rws[i-1] = tmp;

        tmp = state->w_ts[iw_hot];
        state->w_ts[iw_hot] = state->w_ts[iw_cold];
        state->w_ts[iw_cold] = tmp;
      }
    }
  }
  //  check_state_consistency(n_dims, n_Ts, pis, xs, Ts, LRs, Lws, Rws);
  //  printf("bottom of T_swap\n");
}



int check_state_consistency(int n_dims, int n_Ts, double* pis, double** xs,  int* Ts, int* LRs, int* Lws, int* Rws){
  for(int iw = 0; iw < 2*n_Ts; iw++){
    int iT = Ts[iw]; // should be Tindex of iw walker.
    int iLR = LRs[iw]; // should be LR index (0 for L, 1 for R) of iw walker
    if(iLR == 0){ // Left
      assert(iw == Lws[iT]); // Lws
    }else if(iLR == 1){ // Right
      assert(iw == Rws[iT]);
    }else{
      assert(1 == 0);
    }
  }
}
       
