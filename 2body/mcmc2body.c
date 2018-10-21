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

gsl_rng* g_rng;
double Lsq;

// ************* main: ******************

int main(int argc, char* argv[]){

  // ***** process command line args *****
  int n_updates = (argc >= 2)? atoi(argv[1]) : 100000; // number of updates to do
  int n_dims = (argc >= 3)? atoi(argv[2]) : 3; // number of dimensions

  int n_Ts = (argc >= 4)? atoi(argv[3]) : 3; // number of T levels
  if(n_Ts < 1){ n_Ts = 1; } 
  double Thot = (argc >= 5)? atof(argv[4]) : 2.0;
  double Tratio = (n_Ts > 1)? pow(Thot, 1.0/(n_Ts-1)) : 1.0; // Tratio for geometric Temperatures
 

  int n_peaks = (argc >= 6)? atoi(argv[5]) : 2;
  double peak_separation = (argc >= 7)? atof(argv[6]) : 1.0; // spacing parameter; peaks are at x[0] = +-peak_separation/2, other coordinates all zero.
  double peak_width = (argc >= 8)? atof(argv[7]) : 0.2; // width parameter; stddev if peak is Normal.
  double height_ratio = (argc >= 9)? atof(argv[8]) : 0.25; // height of peak i+1 rel to peak i.
  double shape_param = (argc >= 10)? atof(argv[9]) : 1.0; 

  double prop_w = (argc >= 11)? atof(argv[10]) : 0.2;

  int seed = (argc >= 12)? atoi(argv[11])  : 12345; // rng seed

  int n_thin = (argc >= 13)? atoi(argv[12])  : 128; // only output every n_thinth

  char* output_order = (char*)malloc(32*sizeof(char));
  sprintf(output_order, "%s", (argc >= 14)? argv[13] : "unknown"); // 1: walkerwise, 0: T-wise, -1 cold only

  double kernel_scale =  (argc >= 15)? atof(argv[14])  : 0.5*peak_separation;
  // ***** 

  Lsq = pow(kernel_scale, 2.0);

  double* inverse_Temperatures = (double*)malloc(n_Ts*sizeof(double)); // the inverse temperatures 
  inverse_Temperatures[0] = 1.0; // cold
  for(int i = 1; i < n_Ts; i++){
    inverse_Temperatures[i] = inverse_Temperatures[i-1] / Tratio;
  }

  // ***** output the run parameters *****
  printf("# n updates: %5i, n_dims: %2i, n_Ts: %2i, Thot: %6.3f, Tratio: %6.3f Kernel_scale: %5.3f\n", 
         n_updates, n_dims, n_Ts, Thot, Tratio, kernel_scale);
  printf("# n peaks: %3i, peak separation: %5.3f, peak width: %5.3f, height ratio: %5.3f shape param: %5.3f  \n", 
         n_peaks, peak_separation, peak_width, height_ratio, shape_param);

  printf("# proposal width: %5.3f, rng seed: %8i,  thin: %5i, output order: %12s \n", prop_w, seed, n_thin, output_order);


  // ***** Set up RNG *****
  gsl_rng_env_setup();
  const gsl_rng_type* rng_type = gsl_rng_default;
  g_rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(g_rng, seed);


  //*********  Set up target distribution *******************
  peaks* the_peaks = set_up_peaks(n_dims, n_peaks, 0.5*peak_separation, peak_width, height_ratio, shape_param);


  // ***** Setup initial states of walkers *********
  double init_width = 2.0*peak_width;
  chain_state* the_chain_state = set_up_chain_state(n_dims, n_Ts, the_peaks, init_width, inverse_Temperatures);

  // ******************************************
  // ***** Loop through n_updates updates *****
  for(int i = 0; i <= n_updates; i++){ // loop through updates
    if(0){
      for(int j = 0; j < 2*n_Ts; j++){ // update positions of walkers
        update_x(the_peaks, the_chain_state, inverse_Temperatures[the_chain_state->w_ts[j]], prop_w, j); 
      }
        T_swap(the_chain_state, inverse_Temperatures);
    }else{
      for(int j = 0; j < n_Ts; j++){
        update_x_2b(the_peaks, the_chain_state, inverse_Temperatures[j], prop_w, j);
      }
      T_swap_2b(the_chain_state, inverse_Temperatures);
    }
    check_state_consistency(the_chain_state);
    if(i % n_thin == 0){
      printf("%5i ", i);
      if(strcmp(output_order, "T") == 0){
        print_states_T_order(the_chain_state);
      }else if(strcmp(output_order, "walker") == 0){
        print_states_walker_order(the_chain_state);
      }else{
        print_states_cold_only(the_chain_state);
      }
      printf("\n");
    }
  } // end of loop through updates
  // *********************************************

  printf("# ");
  for(int iw = 0; iw < 2*n_Ts; iw++){
    printf("  %2i %5.3f", iw, the_chain_state->w_accepts[iw]/(1.0*n_updates));
  }printf("\n");
  printf("# X-move acceptance rates:");
  for(int it = 0; it < n_Ts; it++){
    printf("  T-level: %1i P_A: %5.3f", it, the_chain_state->t_accepts[it]/(2.0*n_updates)); 
  }printf("\n");
  printf("# T-swap acceptance rates:");
  for(int it = 0; it < n_Ts-1; it++){
    printf("  levels: %1i-%1i P_A: %5.3f", it, it+1, the_chain_state->t_Tswap_accepts[it]/(2.0*n_updates)); 
  }printf("\n");
}

// end of main


// ***** functions *****


// the end
