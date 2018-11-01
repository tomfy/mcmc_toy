#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "target.h"
#include "chain_architecture.h"
#include "chain_state.h"
#include "mcmc2body.h"

gsl_rng* g_rng;
// long pi_evaluation_count;

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

 double prop_w_factor = (argc >= 11)? atof(argv[10]) : 2.0;

 int seed = (argc >= 12)? atoi(argv[11])  : 12345; // rng seed

 int n_thin = (argc >= 13)? atoi(argv[12])  : 128; // only output every n_thinth

 char* output_order = (char*)malloc(32*sizeof(char));
 sprintf(output_order, "%s", (argc >= 14)? argv[13] : "unknown"); // 1: walkerwise, 0: T-wise, -1 cold only

 double kernel_scale =  (argc >= 15)? atof(argv[14])  : 0.5*peak_separation;

 double two_body_interpolation_power = (argc >= 16)? atof(argv[15]) : 0; // 0 => geometric, 1 => linear
 int n_per_level = (argc >= 17)? atoi(argv[16]) : 2;
 int symmetry = (argc >= 18)? atoi(argv[17]) : 1;

 long burn_in = 10000;

 // ******  

 double* inverse_Temperatures = (double*)malloc(n_Ts*sizeof(double)); // the inverse temperatures 
 inverse_Temperatures[0] = 1.0; // cold
 for(int i = 1; i < n_Ts; i++){
   inverse_Temperatures[i] = inverse_Temperatures[i-1] / Tratio;
 }

 
 // ***** output the run parameters *****
 printf("# n_per_level: %i \n", n_per_level);
 printf("# n updates: %5i, n_dims: %2i, n_Ts: %2i, Thot: %6.3f, Tratio: %6.3f Kernel_scale: %5.3f\n", 
        n_updates, n_dims, n_Ts, Thot, Tratio, kernel_scale);
 printf("# n peaks: %3i, peak separation: %5.3f, peak width: %5.3f, height ratio: %5.3f shape param: %5.3f  \n", 
        n_peaks, peak_separation, peak_width, height_ratio, shape_param);

 printf("# proposal width factor: %5.3f, rng seed: %8i,  thin: %5i, output order: %12s \n", prop_w_factor, seed, n_thin, output_order);

 printf("# interpolation power: %8.6f \n", two_body_interpolation_power);
 printf("# symmetry: %5i \n", symmetry);



 // ***** Set up RNG *****
 gsl_rng_env_setup();
 const gsl_rng_type* rng_type = gsl_rng_default;
 g_rng = gsl_rng_alloc(rng_type);
 gsl_rng_set(g_rng, seed);


 //*********  Set up target distribution *******************
 target* the_target = set_up_target(n_dims, n_peaks, 0.5*peak_separation, peak_width, height_ratio, shape_param);

 //********** Set up the chain architecture ****************
 chain_architecture* the_arch = set_up_chain_architecture(n_per_level, symmetry, n_Ts, Thot, prop_w_factor*peak_width, prop_w_factor*kernel_scale, peak_width, kernel_scale, two_body_interpolation_power);
 print_chain_architecture_info(the_arch);

 // ***** Setup initial states of walkers *********
 double init_width = 2.0*peak_width;
 chain_state* the_chain_state = set_up_chain_state(n_dims, n_Ts, the_target, the_arch, init_width);


 check_state_consistency(the_target, the_arch, the_chain_state);


 // *************************************************************
 // ********* Loop through n_updates updates ********************

 if(n_per_level == 1){ // standard 1-body heating.
   for(int i = 0; i <= n_updates; i++){ // loop through updates

     step_1b(the_target, the_arch, the_chain_state);
 
     if(i % n_thin == 0){
       // check_state_consistency(the_target, the_chain_state, inverse_Temperatures);
       printf("%5i ", i);
       if(strcmp(output_order, "T") == 0){
         print_states_T_order(the_chain_state);
       }else if(strcmp(output_order, "walker") == 0){
         print_states_walker_order(the_chain_state);
       }else{
         print_states_L0_Rtop_only(the_chain_state);
       }
       printf("\n");
     }
   }
 }else{ // 2-body auxiliary distributions
   for(int i = 0; i <= n_updates; i++){ // loop through updates

     step_2b(the_target, the_arch, the_chain_state);
     
     if(the_arch->symmetry == 1){
       cold_transition_observe_and_count_sym(the_chain_state);
     }else if(the_arch->symmetry == 0){
       cold_transition_observe_and_count_asym(the_chain_state);
     }
    
     if(i % n_thin == 0){
       check_state_consistency(the_target, the_arch, the_chain_state);
       printf("%5i ", i);
       if(strcmp(output_order, "T") == 0){
         print_states_T_order(the_chain_state);
       }else if(strcmp(output_order, "walker") == 0){
         print_states_walker_order(the_chain_state);
       }else{
         print_states_L0_Rtop_only(the_chain_state);
       }
       printf("\n");
     }
   }
 } // end of loop through updates
 // *********************************************

 printf("# ");
 for(int iw = 0; iw < 2*n_Ts; iw++){
   printf("  %2i %5.3f", iw, the_chain_state->w_accepts[iw]/(1.0*n_updates));
 }printf("\n");
 printf("# X-move acceptance rates:\n");
 for(int it = 0; it < n_Ts; it++){
   if(n_per_level == 2){
     printf("#   T-level: %1i P_A: L, R, avg: %5.3f %5.3f \n", it,          
            the_chain_state->t_Laccepts[it]/(1.0*n_updates),
            the_chain_state->t_Raccepts[it]/(1.0*n_updates)
            ); 
   }else{
     printf("#   T-level: %1i P_A:  %5.3f %5.3f %5.3f \n", it,   
            the_chain_state->t_Laccepts[it]/(1.0*n_updates),
            the_chain_state->t_Raccepts[it]/(1.0*n_updates),       
            the_chain_state->t_accepts[it]/(2.0*n_updates) );
   }
 }
 printf("# T-swap acceptance rates:\n");
 for(int it = 0; it < n_Ts-1; it++){
   printf("#   T-levels: %1i-%1i P_A L,R,both: %5.3f %5.3f %5.3f \n", it, it+1, 
          the_chain_state->t_Tswap_Laccepts[it]/(1.0*n_updates),
          the_chain_state->t_Tswap_Raccepts[it]/(1.0*n_updates),
          the_chain_state->t_Tswap_accepts[it]/(2.0*n_updates)); 
 }
 printf("# LR-swap acceptance rates:\n");
 for(int it = 0; it < n_Ts; it++){
   printf("#   T-level: %1i P_A: %5.3f \n", it,          
          the_chain_state->t_LRswap_accepts[it]/(1.0*n_updates)
          ); 
 }
 printf("# pi evaluation count:  %ld \n", the_target->pi_evaluation_count);
}

// end of main


// ***** functions *****

double Dsquared(int n_dims, double* x, double* y){
  double dsq = 0;
  for(int i = 0; i < n_dims; i++){
    dsq += (x[i] - y[i])*(x[i] - y[i]);
  }
  return dsq;
}



// the end
