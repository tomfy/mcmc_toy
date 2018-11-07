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
//#include "mcmc2body.h"

// function definitions:

chain_state* set_up_chain_state(const target* const the_target, const chain_architecture* const arch, double init_width){
  chain_state* state = (chain_state*)malloc(sizeof(chain_state)); 
  int n_dims = the_target->n_dims;
  int n_levels = arch->n_levels;
  int n_walkers = 2*n_levels;
  state->n_dims = n_dims;
  state->n_levels = n_levels;
  state->n_walkers = n_walkers;
  state->w_xs = (double**)malloc(n_walkers*sizeof(double*));
  state->w_pis = (double*)malloc(n_walkers*sizeof(double)); 
  state->w_ts = (int*)malloc(n_walkers*sizeof(int));
  state->w_LRs = (int*)malloc(n_walkers*sizeof(int));
  state->t_Lws = (int*)malloc(n_levels*sizeof(int));
  state->t_Rws = (int*)malloc(n_levels*sizeof(int));
  state->t_PIs = (double*)malloc(n_levels*sizeof(double));
  state->t_dsqrs = (double*)malloc(n_levels*sizeof(double));

  state->updates = 0;
 state->t_Lxsums = (double**)malloc(n_levels*sizeof(double*));
 state->t_Rxsums = (double**)malloc(n_levels*sizeof(double*));
  state->w_near_peak = (int*)calloc(n_walkers, sizeof(int)); 
  state->w_transition_counts = (int*)calloc(n_walkers, sizeof(int));
  state->t_Lnear_peak = (int*)calloc(n_levels, sizeof(int)); 
  state->t_Ltransition_counts = (int*)calloc(n_levels, sizeof(int));
  state->t_Rnear_peak = (int*)calloc(n_levels, sizeof(int)); 
  state->t_Rtransition_counts = (int*)calloc(n_levels, sizeof(int));

  state->w_accepts = (int*)calloc(n_walkers, sizeof(int)); 
  state->t_accepts = (int*)calloc(n_levels, sizeof(int)); 
  state->t_Laccepts = (int*)calloc(n_levels, sizeof(int)); 
  state->t_Raccepts = (int*)calloc(n_levels, sizeof(int)); 
  state->t_Tswap_accepts = (int*)calloc(n_levels, sizeof(int));
  state->t_Tswap_Laccepts = (int*)calloc(n_levels, sizeof(int));
  state->t_Tswap_Raccepts = (int*)calloc(n_levels, sizeof(int));
  state->t_LRswap_accepts = (int*)calloc(n_levels, sizeof(int));
  for(int i=0; i < n_walkers; i++){
    double* pt = (double*)malloc(n_dims*sizeof(double));
    for(int j=0; j < n_dims; j++){
      pt[j] = gsl_ran_gaussian(g_rng, init_width);
    }
    state->w_xs[i] = pt;
    state->w_pis[i] = pi(the_target, pt);
  
    state->w_ts[i] = i / 2; //  initially e.g. n_levels = 3:  Ts = (0,0,1,1,2,2)
    state->w_LRs[i] =  i % 2; //  initially e.g. n_levels = 3:  LRs = (0,1,0,1,0,1)    
    if(state->w_LRs[i] == 0){
      state->t_Lws[state->w_ts[i]] = i;
    }else if(state->w_LRs[i] == 1){
      state->t_Rws[state->w_ts[i]] = i;
    }else{
      printf("LRs[i] is neither 0 nor 1. Bye.\n");
      exit(0);
    }
    
  } // end loop through walkers
  for(int it = 0; it < n_levels; it++){
    int iwx = state->t_Lws[it];
    int iwy = state->t_Rws[it];
    double* x = state->w_xs[iwx];
    double* y = state->w_xs[iwy];
    double pix = state->w_pis[iwx];
    double piy = state->w_pis[iwy];
    double Dsqr = Dsquared(n_dims, x, y);
    double Lsqr = arch->kernel_widths[it];
    Lsqr *= Lsqr;
    state->t_dsqrs[it] = Dsqr;
    state->t_PIs[it] = PI(pix, piy, it, n_levels, Kernel(Dsqr, Lsqr), arch);
    state->t_Lnear_peak[it] = 0; // -1: peak at -1, 1: peak at 1.
    state->t_Rnear_peak[it] = 0; // -1: peak at -1, 1: peak at 1.
    state->t_Lxsums[it] = (double*)calloc(n_dims, sizeof(double));
    state->t_Rxsums[it] = (double*)calloc(n_dims, sizeof(double));
  }
  return state;
}

void reset_chain_state(chain_state* state){ // reset the cumulative variables to initial values
  int n_dims = state->n_dims;
  int n_walkers = state->n_walkers;
  int n_levels = state->n_levels;

  state->updates = 0;

  free(state->w_transition_counts);
  state->w_transition_counts = (int*)calloc(n_walkers, sizeof(int));
  free(state->t_Ltransition_counts);
  state->t_Ltransition_counts = (int*)calloc(n_levels, sizeof(int));
  free(state->t_Rtransition_counts);
  state->t_Rtransition_counts = (int*)calloc(n_levels, sizeof(int));

  free(state->w_accepts);
  state->w_accepts = (int*)calloc(n_walkers, sizeof(int));
  free(state->t_accepts);
  state->t_accepts = (int*)calloc(n_levels, sizeof(int));
  free(state->t_Laccepts);
  state->t_Laccepts = (int*)calloc(n_levels, sizeof(int));
  free(state->t_Raccepts);
  state->t_Raccepts = (int*)calloc(n_levels, sizeof(int));
  free(state->t_Tswap_accepts);
  state->t_Tswap_accepts = (int*)calloc(n_levels, sizeof(int));
  free(state->t_Tswap_Laccepts);
  state->t_Tswap_Laccepts = (int*)calloc(n_levels, sizeof(int));
  free(state->t_Tswap_Raccepts);
  state->t_Tswap_Raccepts = (int*)calloc(n_levels, sizeof(int));
  free(state->t_LRswap_accepts);
  state->t_LRswap_accepts = (int*)calloc(n_levels, sizeof(int));

  for(int it=0; it<n_levels; it++){
    free(state->t_Lxsums[it]);
    state->t_Lxsums[it] = (double*)calloc(n_dims, sizeof(double));
    free(state->t_Rxsums[it]);
    state->t_Rxsums[it] = (double*)calloc(n_dims, sizeof(double));
  }

}

void print_states_walker_order(const chain_state* const state){
  int n_w = 2*state->n_levels;
  for(int iw=0; iw < n_w; iw++){
    int it = state->w_ts[iw];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");
  }
}

void print_states_T_order(const chain_state* const state){
  
  for(int it=0; it < state->n_levels; it++){
    int iw = state->t_Lws[it];
    printf("%2i %2i %3i %3i  ", iw, it, state->w_transition_counts[iw], state->t_Ltransition_counts[it]);
    for(int j=0; j < state->n_dims; j++){
      printf("%7.5f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
  for(int it=0; it < state->n_levels; it++){
    int iw = state->t_Rws[it];
    printf("%2i %2i %3i %3i  ", iw, it, state->w_transition_counts[iw], state->t_Rtransition_counts[it]);
    for(int j=0; j < state->n_dims; j++){
      printf("%7.5f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
}

void print_states_cold_only(const chain_state* const state){ // print L and R T-level 0
  int it = 0;
  {
    int iw = state->t_Lws[it];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
  {
    int iw = state->t_Rws[it];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
}

void print_states_L0_Rtop_only(const chain_state* const state){ // print L level 0, R level n_levels-1
  int it = 0;
  {
    int iw = state->t_Lws[it];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");    
  }
  it = state->n_levels-1;
  {
    int iw = state->t_Rws[it];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
}

int check_state_consistency(const target* const the_target, const chain_architecture* const arch, const chain_state* const state){
  for(int iw = 0; iw < 2*state->n_levels; iw++){
    int it = state->w_ts[iw]; // should be Tindex of iw walker.
    int iLR = state->w_LRs[iw]; // should be LR index (0 for L, 1 for R) of iw walker
    if(iLR == 0){ // Left
      assert(iw == state->t_Lws[it]); // Lws
    }else if(iLR == 1){ // Right
      assert(iw == state->t_Rws[it]);
    }else{
      assert(1 == 0);
    }
  }
  for(int it=0; it < state->n_levels; it++){
    double pix = state->w_pis[state->t_Lws[it]];
    double piy = state->w_pis[state->t_Rws[it]];
    double dsqr = state->t_dsqrs[it];
    double Lsqr = arch->kernel_widths[it];
    Lsqr *= Lsqr;
    double Kvalue = Kernel(dsqr, Lsqr);

    if(0){
      double* x = state->w_xs[state->t_Lws[it]];
      double* y = state->w_xs[state->t_Rws[it]];
      double pixc = pi(the_target, x);
      double piyc = pi(the_target, y);
      double dsqrc = Dsquared(state->n_dims, x, y);
    
      if( fabs(pixc - pix)/(pixc + pix + 1e-300) < 1e-10){
      }else{
        printf("pix, pixc:  %20.12g %20.12g %20.12g \n", pix, pixc, pixc-pix);

        exit(1);
      }
      if( fabs(piyc - piy)/(piyc + piy + 1e-300) < 1e-10){
      }else{
        printf("piy, piyc:  %20.12g %20.12g %20.12g \n", piy, piyc, piyc-piy);

        exit(1);
      }
      assert( fabs(piyc - piy)/(piyc + piy + 1e-300)  < 1e-10);
      assert( fabs(dsqrc - dsqr)/(dsqrc + dsqr)  < 1e-10);
    }

    // double PI(double pix, double piy, int it, int n_levels, double K, double Tinverse){
    double PIc = PI(pix, piy, it, state->n_levels, Kvalue, arch);
    double fracdiff = fabs(PIc - state->t_PIs[it])/(PIc + state->t_PIs[it]);
    if(fracdiff > 1.0e-10){

      printf("In check state consistency. it: %i  %10.7g %10.7g  %10.7g %10.7g \n", it, pix, piy, PIc, state->t_PIs[it]);
      printf("dsqr: %g, Lsqr: %g  K: %g \n", dsqr, Lsqr, Kvalue);
      assert(fracdiff < 1.0e-10);
    }
  }
}
// end  chain_state  functions


// ************************ 1-body updates *************************************

void update_x(const target* const the_target, const chain_architecture* const arch, chain_state* state, int iw){
  double pi_x = state->w_pis[iw]; 
  double* prop_x = (double*)malloc(state->n_dims*sizeof(double));
  int it = state->w_ts[iw];
  double prop_w = (state->w_LRs[iw] == 0)? arch->Lprop_widths[it] : arch->Rprop_widths[it];
  double Tinverse = arch->inverse_Temperatures[it];
  //  printf("%i  %g  %g \n", it, Tinverse, prop_w);

  for(int i = 0; i < state->n_dims; i++){ // propose a new position 
    prop_x[i] = state->w_xs[iw][i] 
      + gsl_ran_gaussian(g_rng, prop_w); // normal (aka gaussian) proposal
  }

  double pi_prop_x = pi(the_target, prop_x);
  if( (pi_prop_x > pi_x)  ||  (gsl_rng_uniform(g_rng)  <  pow(pi_prop_x/pi_x, Tinverse)) ){ // Accept proposal
    {
   int old_near_peak = (state->w_xs[iw][0] <= 0)? -1 : 1;
      int new_near_peak = (prop_x[0] <= 0)? -1 : 1;
      if(new_near_peak != old_near_peak){ state->t_Ltransition_counts[it]++; }
    }
    state->w_pis[iw] = pi_prop_x;
    free(state->w_xs[iw]);  
    state->w_xs[iw] = prop_x;
    state->w_accepts[iw]++;
    state->t_accepts[state->w_ts[iw]]++;
    if(state->w_LRs[iw] == 0){
      state->t_Laccepts[state->w_ts[iw]]++;
    }else{
      state->t_Raccepts[state->w_ts[iw]]++;
    }
    int new_near_peak = (prop_x[0] < 0)? -1 : 1;
    if ((state->w_ts[iw] == 0 /* cold*/ ) && (state->w_near_peak[iw] != new_near_peak)){
      state->w_transition_counts[iw]++;
      state->w_near_peak[iw] = new_near_peak;
    }
  }else{ // Reject proposal
    free(prop_x);
  }
}

void T_swap(chain_state* state, double* inverse_Temperatures){
  for(int i=1; i<state->n_levels; i++){ // T-swapping among Ls
    double Tc_inverse = inverse_Temperatures[i-1];  // 1/T for the colder of two walkers
    double Th_inverse = inverse_Temperatures[i];  // 1/T for the hotter of the two walkers
    {  //  first the 'left' set of walkers:
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
        //   printf("T-swap accepted (L). T-levels: %i %i \n", i-1, i);
        state->t_Tswap_accepts[i-1]++; // indices used are 0 through n_levels-2
      }
    }
    {  //  then the 'right' set of walkers:
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

        state->t_Tswap_accepts[i-1]++;
      }
    }
  }
}

// ****************** 2-body updates ************************************

void update_x_2b(const target* const the_target, const chain_architecture* const arch, chain_state* state, int it){

  double Lprop_w = arch->Lprop_widths[it];
  double Rprop_w = arch->Rprop_widths[it];
  double Lsqr = arch->kernel_widths[it] * arch->kernel_widths[it];
  double Tinverse = arch->inverse_Temperatures[it];

  int iw_x = state->t_Lws[it];
  int iw_y = state->t_Rws[it];
  double pi_x = state->w_pis[iw_x];
  double pi_y = state->w_pis[iw_y];

  // update x
  {
    double* propx = (double*)malloc(state->n_dims*sizeof(double));  
    for(int i = 0; i < state->n_dims; i++){ // propose a new position 
      propx[i] = state->w_xs[iw_x][i] + gsl_ran_gaussian(g_rng, Lprop_w); // normal (aka gaussian) proposal
    }
    double pi_propx = pi(the_target, propx);

    double PI_x_y = state->t_PIs[it];
    double Dsqr_prop = Dsquared(state->n_dims, propx, state->w_xs[iw_y]);
    double PI_propx_y = PI(pi_propx, pi_y, it, state->n_levels, Kernel(Dsqr_prop, Lsqr), arch); 
    if( (PI_propx_y > PI_x_y)  ||  (gsl_rng_uniform(g_rng)*PI_x_y  <  PI_propx_y) ){ // Accept proposal

      int old_near_peak = (state->w_xs[iw_x][0] <= 0)? -1 : 1;
      int new_near_peak = (propx[0] <= 0)? -1 : 1;
      if(new_near_peak != old_near_peak){ state->t_Ltransition_counts[it]++; }

      state->w_pis[iw_x] = pi_propx;
      free(state->w_xs[iw_x]);  
      state->w_xs[iw_x] = propx;

      state->w_accepts[iw_x]++;
      state->t_Laccepts[it]++;

      state->t_dsqrs[it] = Dsqr_prop;
      state->t_PIs[it] = PI_propx_y;
    }else{ // Reject proposal
      free(propx);
    }
  }
  // update y
  {  
    pi_x = state->w_pis[iw_x]; // x and pi_x may have changed above, so get the up-to-date value!
    double* propy = (double*)malloc(state->n_dims*sizeof(double));
    for(int i = 0; i < state->n_dims; i++){ // propose a new position
      propy[i] = state->w_xs[iw_y][i] + gsl_ran_gaussian(g_rng, Rprop_w); // normal (aka gaussian) proposal
    }
    double pi_propy = (arch->symmetry == 1)? pi(the_target, propy) : 1.0; // for asymmetric, don't need to calculate pi(y) here (do later if accepted)

    double PI_x_y = state->t_PIs[it];
    double Dsqr_prop = Dsquared(state->n_dims, state->w_xs[iw_x], propy);
    double PI_x_propy = PI(pi_x, pi_propy, it, state->n_levels, Kernel(Dsqr_prop, Lsqr), arch); 
    if( (PI_x_propy > PI_x_y)  ||  (gsl_rng_uniform(g_rng)*PI_x_y <  PI_x_propy) ){ // Accept proposal

      int old_near_peak = (state->w_xs[iw_y][0] <= 0)? -1 : 1;
      int new_near_peak = (propy[0] <= 0)? -1 : 1;
      if(new_near_peak != old_near_peak){ state->t_Rtransition_counts[it]++; }

      state->w_pis[iw_y] = (arch->symmetry)? pi_propy : pi(the_target, propy); // if asymmetric, didn't get pi(y) before, do now.
      free(state->w_xs[iw_y]);  
      state->w_xs[iw_y] = propy;

      state->w_accepts[iw_y]++;
      //   state->t_accepts[it]++;
      state->t_Raccepts[it]++;

      state->t_dsqrs[it] = Dsqr_prop;
      state->t_PIs[it] = PI_x_propy;
    }else{ // Reject proposal
      free(propy);
    }
  }
}

void cold_transition_observe_and_count_sym(chain_state* state){
  // a walker's nearest peak is only reset when the walker is cold (i.e. x, it=0, or y, it=nTs-1
  // a transition occurs if previous near peak (when walker was cold) is different from latest near peak.
  int it = 0;
  // L (i.e. x)
  int iw = state->t_Lws[it];
  int new_near_peak = (state->w_xs[iw][0] < 0)? -1 : 1;
  if(new_near_peak != state->w_near_peak[iw]){
    state->w_transition_counts[iw]++;
    state->w_near_peak[iw] = new_near_peak;
  }
  // R (i.e. y)
  it = state->n_levels-1;
  iw = state->t_Rws[it];
  new_near_peak = (state->w_xs[iw][0] < 0)? -1 : 1;
  if(new_near_peak != state->w_near_peak[iw]){
    state->w_transition_counts[iw]++;
    state->w_near_peak[iw] = new_near_peak;
  }
}

void cold_transition_observe_and_count_asym(chain_state* state){
  // asymmetrical; L <-> cold.
  // a walker's nearest peak is only reset when the walker is cold (i.e. x, it=0, or y, it=nTs-1
  // a transition occurs if previous near peak (when walker was cold) is different from latest near peak.
  for(int it = 0; it < state->n_levels; it++){
  // L (i.e. x)
  int iw = state->t_Lws[it];
  int new_near_peak = (state->w_xs[iw][0] < 0)? -1 : 1;
  if(new_near_peak != state->w_near_peak[iw]){
    state->w_transition_counts[iw]++;
    state->w_near_peak[iw] = new_near_peak;
  }
  }

}

void T_swap_2b_A(const chain_architecture* const arch, chain_state* state){ // swap both x and y between T levels at same time.

  for(int i=1; i<state->n_levels; i++){ // T-swapping among Ls
    //  double Tc_inverse = inverse_Temperatures[i-1];  // 1/T for the colder of two walkers
    //  double Th_inverse = inverse_Temperatures[i];  // 1/T for the hotter of the two walkers
  
    int iw_xhot = state->t_Lws[i]; // index of walker with ith temperature
    int iw_xcold = state->t_Lws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
    int iw_yhot = state->t_Rws[i]; // index of walker with ith temperature
    int iw_ycold = state->t_Rws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
    double pi_xcold = state->w_pis[iw_xcold];
    double pi_xhot = state->w_pis[iw_xhot];
    double pi_ycold = state->w_pis[iw_ycold];
    double pi_yhot = state->w_pis[iw_yhot]; 
    double Lsqr_cold = arch->kernel_widths[i-1] * arch->kernel_widths[i-1];
    double Lsqr_hot = arch->kernel_widths[i] * arch->kernel_widths[i];
   
    //   printf("PI_prop_cold.  it cold: %i \n", i-1);
    double PI_prop_cold = PI(pi_xhot, pi_yhot, i-1, state->n_levels, 
                             Kernel(state->t_dsqrs[i], Lsqr_cold), 
                             arch);
    double PI_cold_ratio = PI_prop_cold / state->t_PIs[i-1];

    //   double PI_current_cold =   PI(pi_xcold, pi_ycold, i-1, state->n_levels,
    //                              Kernel(Dsquared(state->n_dims, state->w_xs[iw_xcold], state->w_xs[iw_ycold]), Lsq),
    //                              Tc_inverse);
    //   printf("current PI cold:  %8.6g %8.6g \n", PI_current_cold, state->t_PIs[i-1]);
    //    double fracdiff = (fabs(PI_current_cold - state->t_PIs[i-1])/(PI_current_cold + state->t_PIs[i-1]));
    //   printf("fracdiff: %10.8g \n", fracdiff);

    //    printf("PI_prop_hot. it hot: %i\n", i);
    double PI_prop_hot = PI(pi_xcold, pi_ycold, i, state->n_levels,  
                            Kernel(state->t_dsqrs[i-1], Lsqr_hot),
                            arch); 
      
    double PI_hot_ratio = PI_prop_hot / state->t_PIs[i];
    // double PI_current_hot =    PI(pi_xhot, pi_yhot, i, state->n_levels,
    //                              Kernel(Dsquared(state->n_dims, state->w_xs[iw_xhot], state->w_xs[iw_yhot]), Lsq),
    //                              Th_inverse);
    
    double PI_ratio = PI_cold_ratio * PI_hot_ratio;

    if( (PI_ratio > 1.0)  ||  (gsl_rng_uniform(g_rng)  <  PI_ratio) ) { // Accept proposed T-swap

      // swap colder and hotter L walkers:
      int tmp = state->t_Lws[i]; 
      state->t_Lws[i] = state->t_Lws[i-1];
      state->t_Lws[i-1] = tmp;

      tmp = state->w_ts[iw_xhot];
      state->w_ts[iw_xhot] = state->w_ts[iw_xcold];
      state->w_ts[iw_xcold] = tmp;

      // swap colder and hotter R walkers
      tmp = state->t_Rws[i]; 
      state->t_Rws[i] = state->t_Rws[i-1];
      state->t_Rws[i-1] = tmp;

      tmp = state->w_ts[iw_yhot];
      state->w_ts[iw_yhot] = state->w_ts[iw_ycold];
      state->w_ts[iw_ycold] = tmp;

      state->t_Tswap_accepts[i-1]++; // indices used are 0 through n_levels-2
      //  state->t_Tswap_Laccepts[i-1]++;

      double tmp_dsqr = state->t_dsqrs[i-1];
      state->t_dsqrs[i-1] = state->t_dsqrs[i];
      state->t_dsqrs[i] = tmp_dsqr;

      state->t_PIs[i-1] = PI_prop_cold;
      state->t_PIs[i] = PI_prop_hot;
    } // end of if accepted
  } // loop through T levels
} // end T_swap_2b


void T_swap_2b_B(const chain_architecture* arch, chain_state* state){ // swap xs, then ys
  for(int i=1; i<state->n_levels; i++){ // T-swapping among Ls
    //  double Tc_inverse = inverse_Temperatures[i-1];  // 1/T for the colder of two walkers
    //  double Th_inverse = inverse_Temperatures[i];  // 1/T for the hotter of the two walkers

      double Lsqr_cold = arch->kernel_widths[i-1] * arch->kernel_widths[i-1];
      double Lsqr_hot = arch->kernel_widths[i] * arch->kernel_widths[i];

    {  //  first the 'left' set of walkers, i.e. x's:
      int iw_xhot = state->t_Lws[i]; // index of walker with ith temperature
      int iw_xcold = state->t_Lws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
      int iw_yhot = state->t_Rws[i]; // index of walker with ith temperature
      int iw_ycold = state->t_Rws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
      double pi_xcold = state->w_pis[iw_xcold];
      double pi_xhot = state->w_pis[iw_xhot];
      double pi_ycold = state->w_pis[iw_ycold];
      double pi_yhot = state->w_pis[iw_yhot];
      double Lsqr_cold = arch->kernel_widths[i-1] * arch->kernel_widths[i-1];
      double Lsqr_hot = arch->kernel_widths[i] * arch->kernel_widths[i];

      // proposed colder state, with xhot, ycold, i.e x swapped
      double dsqr_prop_cold = Dsquared(state->n_dims, state->w_xs[iw_xhot], state->w_xs[iw_ycold]);
      double PI_prop_cold = PI(pi_xhot, pi_ycold, i-1, state->n_levels, 
                               Kernel(dsqr_prop_cold, Lsqr_cold), 
                               arch);
      double PI_cold_ratio = PI_prop_cold / state->t_PIs[i-1];
      /* PI(pi_xcold, pi_ycold, i-1, state->n_levels,  */
      /*    Kernel(Dsquared(state->n_dims, state->w_xs[iw_xcold], state->w_xs[iw_ycold]), Lsq), */
      /*    Tc_inverse); */

   // proposed hotter state, with xcold, yhot, i.e x swapped
      double dsqr_prop_hot = Dsquared(state->n_dims, state->w_xs[iw_xcold], state->w_xs[iw_yhot]);
      double PI_prop_hot = PI(pi_xcold, pi_yhot, i, state->n_levels,  
                              Kernel(dsqr_prop_hot, Lsqr_hot),
                              arch); 

    
      
   
      double PI_hot_ratio = PI_prop_hot / state->t_PIs[i];
      /* PI(pi_xhot, pi_yhot, i, state->n_levels,  */
      /*    Kernel(Dsquared(state->n_dims, state->w_xs[iw_xhot], state->w_xs[iw_yhot]), Lsq), */
      /*    Th_inverse); */
    
      double PI_ratio = PI_cold_ratio * PI_hot_ratio;

      //   printf("%5i  %7.5g  %7.5g  %7.5g  %7.5g   %7.5g  %7.5g  %7.5g \n", i-1, state->t_PIs[i-1], state->t_PIs[i], PI_prop_cold, PI_prop_hot,  PI_cold_ratio, PI_hot_ratio, PI_ratio);
      
      if( (PI_ratio > 1.0)  ||  (gsl_rng_uniform(g_rng)  <  PI_ratio) ) { // Accept proposed T-swap
        int tmp = state->t_Lws[i]; 
        state->t_Lws[i] = state->t_Lws[i-1];
        state->t_Lws[i-1] = tmp;

        tmp = state->w_ts[iw_xhot];
        state->w_ts[iw_xhot] = state->w_ts[iw_xcold];
        state->w_ts[iw_xcold] = tmp;
        //   state->t_Tswap_accepts[i-1]++; // indices used are 0 through n_levels-2
        state->t_Tswap_Laccepts[i-1]++;

        state->t_dsqrs[i-1] = dsqr_prop_cold;
        state->t_dsqrs[i] = dsqr_prop_hot;
        state->t_PIs[i-1] = PI_prop_cold;
        state->t_PIs[i] = PI_prop_hot;
      }
    }
    {  //  then the 'right' set of walkers (i.e. swap y's):
      int iw_xhot = state->t_Lws[i]; // index of walker with ith temperature
      int iw_xcold = state->t_Lws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
      int iw_yhot = state->t_Rws[i]; // index of walker with ith temperature
      int iw_ycold = state->t_Rws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
      double pi_xcold = state->w_pis[iw_xcold];
      double pi_xhot = state->w_pis[iw_xhot];
      double pi_ycold = state->w_pis[iw_ycold];
      double pi_yhot = state->w_pis[iw_yhot];

      // proposed colder state, with xcold, yhot, i.e y swapped
      double dsqr_prop_cold = Dsquared(state->n_dims, state->w_xs[iw_xcold], state->w_xs[iw_yhot]);
      double PI_prop_cold = PI(pi_xcold, pi_yhot, i-1, state->n_levels, 
                               Kernel(dsqr_prop_cold, Lsqr_cold), 
                               arch);
      double PI_cold_ratio = PI_prop_cold / state->t_PIs[i-1];

      // proposed hotter state, with xhot, ycold, i.e y swapped
      double dsqr_prop_hot = Dsquared(state->n_dims, state->w_xs[iw_xhot], state->w_xs[iw_ycold]);
      double PI_prop_hot = PI(pi_xhot, pi_ycold, i, state->n_levels,  
                              Kernel(dsqr_prop_hot, Lsqr_hot),
                              arch); 
      double PI_hot_ratio = PI_prop_hot / state->t_PIs[i];

      double PI_ratio = PI_cold_ratio * PI_hot_ratio;
      
      if( (PI_ratio > 1.0)  ||  (gsl_rng_uniform(g_rng)  <  PI_ratio)  ) { // Accept proposed T-swap
        int tmp = state->t_Rws[i]; 
        state->t_Rws[i] = state->t_Rws[i-1];
        state->t_Rws[i-1] = tmp;

        tmp = state->w_ts[iw_yhot];
        state->w_ts[iw_yhot] = state->w_ts[iw_ycold];
        state->w_ts[iw_ycold] = tmp;
        //      state->t_Tswap_accepts[i-1]++; // indices used are 0 through n_levels-2
        state->t_Tswap_Raccepts[i-1]++;

        state->t_dsqrs[i-1] = dsqr_prop_cold;
        state->t_dsqrs[i] = dsqr_prop_hot;
        state->t_PIs[i-1] = PI_prop_cold;
        state->t_PIs[i] = PI_prop_hot;
      }
    } // end y swapping block   
  } // loop through T levels
} // end T_swap_2b

void LR_swap(const chain_architecture* arch, chain_state* state){
  for(int it = 0; it < state->n_levels; it++){
    //   printf("it: %8i \n", it);
    int iw_x = state->t_Lws[it];
    int iw_y = state->t_Rws[it];
    double Lsqr = arch->kernel_widths[it]*arch->kernel_widths[it];
    double PI_prop = PI(state->w_pis[iw_y], state->w_pis[iw_x], it, state->n_levels, Kernel(state->t_dsqrs[it], Lsqr), arch);
    if( (PI_prop > state->t_PIs[it])  ||  (gsl_rng_uniform(g_rng)*state->t_PIs[it] < PI_prop) ){ // Accept proposed swap
      assert( state->w_ts[iw_x] == it);
      assert( state->w_ts[iw_y] == it);
      state->w_LRs[iw_x] = 1; // i.e. the walker that was L is now R
      state->w_LRs[iw_y] = 0; // i.e. the walker that was R is now L
      state->t_Lws[it] = iw_y;
      state->t_Rws[it] = iw_x;
      state->t_PIs[it] = PI_prop;
      state->t_LRswap_accepts[it]++;
    }
  }
}


void step_1b(const target* const target, const chain_architecture* const arch, chain_state* state){    
  for(int iw = 0; iw < state->n_walkers; iw++){ // update positions of walkers
    update_x(target, arch, state, iw); 
  }
  T_swap(state, arch->inverse_Temperatures);
  state->updates++;
}


void step_2b(const target* const target, const chain_architecture* const arch, chain_state* state){    
  for(int level = 0; level < arch->n_levels; level++){
    update_x_2b(target, arch, state, level);
  }
  T_swap_2b_A(arch, state); // for each pair of levels: propose swapping both x and y at once, then accept/reject.
  T_swap_2b_B(arch, state); // for each pair of levels: propose swapping x, accept/reject; then propose swapping y, accept/reject.
  if(arch->symmetry == 0){
    LR_swap(arch, state);
  }
  state->updates++;
}

void accumulate_x_sums(chain_state* state){ 
  for(int level = 0; level < state->n_levels; level++){
    int iw_x = state->t_Lws[level];
    int iw_y = state->t_Rws[level];
    for(int i_dim = 0; i_dim < state->n_dims; i_dim++){
      state->t_Lxsums[level][i_dim] += state->w_xs[iw_x][i_dim];
      state->t_Rxsums[level][i_dim] += state->w_xs[iw_y][i_dim];
    }
  }
}

// ************************************************************

// the end





