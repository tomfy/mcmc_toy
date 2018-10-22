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

// peaks  function definitions
peaks* set_up_peaks(int n_dims, int n_peaks, double spacing, double width, double height_ratio, double shape_param){
  peaks* the_peaks = (peaks*)malloc(sizeof(peaks));
  the_peaks->n_peaks = n_peaks;
  the_peaks->peaks = (peak**)malloc(n_peaks*sizeof(peak*)); // array of pointers to peaks.
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
    the_peaks->peaks[i] = a_peak; 
    peak_height *= height_ratio;
  }
  return the_peaks;
}
// end  peaks  functions

//  chain_state  function definitions: 
chain_state* set_up_chain_state(int n_dims, int n_Ts, peaks* the_peaks, double init_width, double* inverse_Temperatures){
  chain_state* the_chain_state = (chain_state*)malloc(sizeof(chain_state));
  int n_walkers = 2*n_Ts;
  the_chain_state->n_dims = n_dims;
  the_chain_state->n_Ts = n_Ts;
  the_chain_state->n_walkers = n_walkers;
  the_chain_state->w_xs = (double**)malloc(n_walkers*sizeof(double*));
  the_chain_state->w_pis = (double*)malloc(n_walkers*sizeof(double)); 
  the_chain_state->w_ts = (int*)malloc(n_walkers*sizeof(int));
  the_chain_state->w_LRs = (int*)malloc(n_walkers*sizeof(int));
  the_chain_state->t_Lws = (int*)malloc(n_Ts*sizeof(int));
  the_chain_state->t_Rws = (int*)malloc(n_Ts*sizeof(int));
  the_chain_state->t_PIs = (double*)malloc(n_Ts*sizeof(double));

  the_chain_state->w_near_peak = (int*)malloc(n_walkers*sizeof(int)); 
  the_chain_state->w_transition_counts = (int*)malloc(n_walkers*sizeof(int));

  the_chain_state->w_accepts = (int*)calloc(n_walkers, sizeof(int)); 
  the_chain_state->t_accepts = (int*)calloc(n_Ts, sizeof(int)); 
  the_chain_state->t_Laccepts = (int*)calloc(n_Ts, sizeof(int)); 
  the_chain_state->t_Raccepts = (int*)calloc(n_Ts, sizeof(int)); 
  the_chain_state->t_Tswap_accepts = (int*)calloc(n_Ts, sizeof(int));
  the_chain_state->t_Tswap_Laccepts = (int*)calloc(n_Ts, sizeof(int));
  the_chain_state->t_Tswap_Raccepts = (int*)calloc(n_Ts, sizeof(int));
  for(int i=0; i < n_walkers; i++){
    double* pt = (double*)malloc(n_dims*sizeof(double));
    for(int j=0; j < n_dims; j++){
        pt[j] = gsl_ran_gaussian(g_rng, init_width);
    }
    the_chain_state->w_xs[i] = pt;
    the_chain_state->w_pis[i] = pi(the_peaks, pt);
  
    the_chain_state->w_ts[i] = i / 2; //  initially e.g. n_Ts = 3:  Ts = (0,0,1,1,2,2)
    the_chain_state->w_LRs[i] =  i % 2; //  initially e.g. n_Ts = 3:  LRs = (0,1,0,1,0,1)    
    if(the_chain_state->w_LRs[i] == 0){
      the_chain_state->t_Lws[the_chain_state->w_ts[i]] = i;
    }else if(the_chain_state->w_LRs[i] == 1){
      the_chain_state->t_Rws[the_chain_state->w_ts[i]] = i;
    }else{
      printf("LRs[i] is neither 0 nor 1. Bye.\n");
      exit(0);
    }

    the_chain_state->w_near_peak[i] = 0; // (pt[0] < 0)? 0 : 1; 
    the_chain_state->w_transition_counts[i] = 0;
    the_chain_state->w_accepts[i] = 0;
    
  } // end loop through walkers
  for(int it = 0; it < n_Ts; it++){
    int iwx = the_chain_state->t_Lws[it];
    int iwy = the_chain_state->t_Rws[it];
    double* x = the_chain_state->w_xs[iwx];
    double* y = the_chain_state->w_xs[iwy];
    double pix = the_chain_state->w_pis[iwx];
    double piy = the_chain_state->w_pis[iwy];
    double Tinverse = inverse_Temperatures[it];
    the_chain_state->t_PIs[it] = PI(pix, piy, it, n_Ts, Kernel(Dsquared(n_dims, x, y), Lsq), Tinverse);
  }
  return the_chain_state;
}

void print_states_walker_order(chain_state* state){
  int n_w = 2*state->n_Ts;
  for(int iw=0; iw < n_w; iw++){
    int it = state->w_ts[iw];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");
  }
}

void print_states_T_order(chain_state* state){
  
  for(int it=0; it < state->n_Ts; it++){
    int iw = state->t_Lws[it];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%7.5f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
  for(int it=0; it < state->n_Ts; it++){
    int iw = state->t_Rws[it];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%7.5f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
}

void print_states_cold_only(chain_state* state){ 
  for(int it=0; it < 1; it++){
    int iw = state->t_Lws[it];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
  for(int it=0; it < 1; it++){
    int iw = state->t_Rws[it];
    printf("%2i %2i %3i  ", iw, it, state->w_transition_counts[iw]);
    for(int j=0; j < state->n_dims; j++){
      printf("%5.3f ", state->w_xs[iw][j]);
    }printf("   ");    
  } 
}

int check_state_consistency(chain_state* state){
  for(int iw = 0; iw < 2*state->n_Ts; iw++){
    int iT = state->w_ts[iw]; // should be Tindex of iw walker.
    int iLR = state->w_LRs[iw]; // should be LR index (0 for L, 1 for R) of iw walker
    if(iLR == 0){ // Left
      assert(iw == state->t_Lws[iT]); // Lws
    }else if(iLR == 1){ // Right
      assert(iw == state->t_Rws[iT]);
    }else{
      assert(1 == 0);
    }
  }
}
// end  chain_state  functions


double pi(peaks* the_peaks, double* x){
  double result = 0.0;
  for(int k = 0; k < the_peaks->n_peaks; k++){   
    double peak_result;
    peak* a_peak = the_peaks->peaks[k];
    double rsq = 0;
      for(int i = 0; i < a_peak->n_dims; i++){
        double dx = (x[i] - a_peak->position[i])/a_peak->width;
        rsq += dx*dx;
      }
    if(1){
      peak_result = exp(-0.5*rsq);
      /* for(int i = 0; i < a_peak->n_dims; i++){ */
      /*   peak_result *= exp(-0.5*pow( (x[i] - a_peak->position[i])/a_peak->width, 2) ); */
      /* } */
    }else{     
      peak_result = pow((1.0 + rsq), -1.0*a_peak->shape); // -0.5*(n_dims + 1));
    }
    result += a_peak->height * peak_result;
  }
  return result;
}


void update_x(peaks* the_peaks, chain_state* state, double Tinverse, double prop_w, int iw){
  double pi_x = state->w_pis[iw]; 
  double* prop_x = (double*)malloc(state->n_dims*sizeof(double));
  prop_w /= sqrt(Tinverse + 1e-4);
  
  for(int i = 0; i < state->n_dims; i++){ // propose a new position 
    prop_x[i] = state->w_xs[iw][i] 
      + gsl_ran_gaussian(g_rng, prop_w); // normal (aka gaussian) proposal
  }

  double pi_prop_x = pi(the_peaks, prop_x);
  if( (pi_prop_x > pi_x)  ||  (gsl_rng_uniform(g_rng)  <  pow(pi_prop_x/pi_x, Tinverse)) ){ // Accept proposal
    state->w_pis[iw] = pi_prop_x;
    free(state->w_xs[iw]);  
    state->w_xs[iw] = prop_x;
    state->w_accepts[iw]++;
    state->t_accepts[state->w_ts[iw]]++;
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
  for(int i=1; i<state->n_Ts; i++){ // T-swapping among Ls
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
        state->t_Tswap_accepts[i-1]++; // indices used are 0 through n_Ts-2
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


// ***************** 2-body updates ***********************

void update_x_2b(peaks* the_peaks, chain_state* state, double Tinverse, double prop_w, int it){
  int iw_x = state->t_Lws[it];
  int iw_y = state->t_Rws[it];
  double pi_x = state->w_pis[iw_x];
  double pi_y = state->w_pis[iw_y];

  // update x
  double* propx = (double*)malloc(state->n_dims*sizeof(double));  
  for(int i = 0; i < state->n_dims; i++){ // propose a new position 
    propx[i] = state->w_xs[iw_x][i] + gsl_ran_gaussian(g_rng, prop_w); // normal (aka gaussian) proposal
  }
  double pi_propx = pi(the_peaks, propx);

  double PI_x_y = state->t_PIs[it];
  double PI_propx_y = PI(pi_propx, pi_y, it, state->n_Ts, Kernel(Dsquared(state->n_dims, propx, state->w_xs[iw_y]), Lsq), Tinverse); 
  if( (PI_propx_y > PI_x_y)  ||  (gsl_rng_uniform(g_rng)*PI_x_y  <  PI_propx_y) ){ // Accept proposal

    state->w_pis[iw_x] = pi_propx;
    free(state->w_xs[iw_x]);  
    state->w_xs[iw_x] = propx;
    state->w_accepts[iw_x]++;
    //    assert(it == state->w_ts[iw_x]);
    state->t_accepts[it]++; // can't we just use it?
    state->t_Laccepts[it]++;
    state->t_PIs[it] = PI_propx_y;
    int new_near_peak = (propx[0] < 0)? -1 : 1;
    if ( ( it == 0 /* cold*/ ) && (state->w_near_peak[iw_x] != new_near_peak)){
      state->w_transition_counts[iw_x]++;
      state->w_near_peak[iw_x] = new_near_peak;
      pi_x = pi_propx;
    }
  }else{ // Reject proposal
    free(propx);
  }

  // update y
  double* propy = (double*)malloc(state->n_dims*sizeof(double));
  for(int i = 0; i < state->n_dims; i++){ // propose a new position
    propy[i] = state->w_xs[iw_y][i] + gsl_ran_gaussian(g_rng, prop_w); // normal (aka gaussian) proposal
  }
  double pi_propy = pi(the_peaks, propy);

  PI_x_y = state->t_PIs[it];
  double PI_x_propy = PI(pi_x, pi_propy, it, state->n_Ts, Kernel(Dsquared(state->n_dims, state->w_xs[iw_x], propy), Lsq), Tinverse); 
  if( (PI_x_propy > PI_x_y)  ||  (gsl_rng_uniform(g_rng)*PI_x_y <  PI_x_propy) ){ // Accept proposal

    state->w_pis[iw_y] = pi_propy;
    free(state->w_xs[iw_y]);  
    state->w_xs[iw_y] = propy;
    state->w_accepts[iw_y]++;
    state->t_accepts[it]++;
    state->t_Raccepts[it]++;
    state->t_PIs[it] = PI_x_propy;
    int new_near_peak = (propy[0] < 0)? -1 : 1;
    if ( ( it == 0 /* cold*/ ) && (state->w_near_peak[iw_y] != new_near_peak)){
      state->w_transition_counts[iw_y]++;
      state->w_near_peak[iw_y] = new_near_peak;
    }
  }else{ // Reject proposal
    free(propy);
  }
}

void T_swap_2b(chain_state* state, double* inverse_Temperatures){
  for(int i=1; i<state->n_Ts; i++){ // T-swapping among Ls
    double Tc_inverse = inverse_Temperatures[i-1];  // 1/T for the colder of two walkers
    double Th_inverse = inverse_Temperatures[i];  // 1/T for the hotter of the two walkers
    {  //  first the 'left' set of walkers, i.e. x's:
      int iw_xhot = state->t_Lws[i]; // index of walker with ith temperature
      int iw_xcold = state->t_Lws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
      int iw_yhot = state->t_Rws[i]; // index of walker with ith temperature
      int iw_ycold = state->t_Rws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature
      double pi_xcold = state->w_pis[iw_xcold];
      double pi_xhot = state->w_pis[iw_xhot];
      double pi_ycold = state->w_pis[iw_ycold];
      double pi_yhot = state->w_pis[iw_yhot];
   
      double PI_prop_cold = PI(pi_xhot, pi_ycold, i-1, state->n_Ts, 
                           Kernel(Dsquared(state->n_dims, state->w_xs[iw_xhot], state->w_xs[iw_ycold]), Lsq), 
                           Tc_inverse);
      double PI_cold_ratio = PI_prop_cold / state->t_PIs[i-1];
        /* PI(pi_xcold, pi_ycold, i-1, state->n_Ts,  */
        /*    Kernel(Dsquared(state->n_dims, state->w_xs[iw_xcold], state->w_xs[iw_ycold]), Lsq), */
        /*    Tc_inverse); */
      double PI_prop_hot = PI(pi_xcold, pi_yhot, i, state->n_Ts,  
                               Kernel(Dsquared(state->n_dims, state->w_xs[iw_xcold], state->w_xs[iw_yhot]), Lsq),
                              Th_inverse); 
      
        double PI_hot_ratio = PI_prop_hot / state->t_PIs[i];
        /* PI(pi_xhot, pi_yhot, i, state->n_Ts,  */
        /*    Kernel(Dsquared(state->n_dims, state->w_xs[iw_xhot], state->w_xs[iw_yhot]), Lsq), */
        /*    Th_inverse); */
    
      double PI_ratio = PI_cold_ratio * PI_hot_ratio;
      
      if( (PI_ratio >= 1.0)  ||  (gsl_rng_uniform(g_rng)  <  PI_ratio) ) { // Accept proposed T-swap
        int tmp = state->t_Lws[i]; 
        state->t_Lws[i] = state->t_Lws[i-1];
        state->t_Lws[i-1] = tmp;

        tmp = state->w_ts[iw_xhot];
        state->w_ts[iw_xhot] = state->w_ts[iw_xcold];
        state->w_ts[iw_xcold] = tmp;
        state->t_Tswap_accepts[i-1]++; // indices used are 0 through n_Ts-2
        state->t_Tswap_Laccepts[i-1]++;
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
    
       double PI_prop_cold = PI(pi_xcold, pi_yhot, i-1, state->n_Ts, 
                                Kernel(Dsquared(state->n_dims, state->w_xs[iw_xcold], state->w_xs[iw_yhot]), Lsq), 
                                Tc_inverse);
      double PI_cold_ratio = PI_prop_cold / state->t_PIs[i-1];
      //   PI(pi_xcold, pi_ycold, i-1, state->n_Ts, 
      //     Kernel(Dsquared(state->n_dims, state->w_xs[iw_xcold], state->w_xs[iw_ycold]), Lsq),
      //     Tc_inverse);
      double PI_prop_hot = PI(pi_xhot, pi_ycold, i, state->n_Ts,  
                               Kernel(Dsquared(state->n_dims, state->w_xs[iw_xhot], state->w_xs[iw_ycold]), Lsq),
                              Th_inverse); 
        double PI_hot_ratio = PI_prop_hot / state->t_PIs[i];
        /* PI(pi_xhot, pi_yhot, i, state->n_Ts,  */
        /*    Kernel(Dsquared(state->n_dims, state->w_xs[iw_xhot], state->w_xs[iw_yhot]), Lsq), */
        /*    Th_inverse); */

      double PI_ratio = PI_cold_ratio * PI_hot_ratio;
      
      if( (PI_ratio >= 1.0)  ||  (gsl_rng_uniform(g_rng)  <  PI_ratio)  ) { // Accept proposed T-swap
        int tmp = state->t_Rws[i]; 
        state->t_Rws[i] = state->t_Rws[i-1];
        state->t_Rws[i-1] = tmp;

        tmp = state->w_ts[iw_yhot];
        state->w_ts[iw_yhot] = state->w_ts[iw_ycold];
        state->w_ts[iw_ycold] = tmp;
        state->t_Tswap_accepts[i-1]++; // indices used are 0 through n_Ts-2
        state->t_Tswap_Raccepts[i-1]++;
        state->t_PIs[i-1] = PI_prop_cold;
        state->t_PIs[i] = PI_prop_hot;
      }
    } // end y swapping block   
  } // loop through T levels
} // end T_swap_2b

// ************************************************************

double Dsquared(int n_dims, double* x, double* y){
  double dsq = 0;
  for(int i = 0; i < n_dims; i++){
    dsq += (x[i] - y[i])*(x[i] - y[i]);
  }
  return dsq;
}

double Kernel(double Dsq, double Lsq){ // convolution kernel
  return exp(-0.5*Dsq/Lsq);
  //  return (Dsq >= Lsq)? 0.0 : 1.0 - Dsq/Lsq;
}

double PI(double pix, double piy, int it, int n_Ts, double K, double Tinverse){
  //  return pow(pix*piy, Tinverse);
  if(it == 0){
    return pix*K;
  }else{
    int itop = n_Ts - 1;
    if(it == itop){
      return piy*K;
    }else{
      //   return pow( ( (1 - 1.0*it/itop)*pix + (1.0*it/itop)*piy ), Tinverse) * K;
      return ( (1 - 1.0*it/itop)*pix + (1.0*it/itop)*piy ) * K;

    }
  }
}
    
// the end





