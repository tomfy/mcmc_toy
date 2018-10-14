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

/* /\* typedef struct{ */
/*   int n_dims; */
/*   double* position; */
/*   double width; */
/*   double height; */
/* }peak; */

/* typedef struct{ */
/*   int n_peaks; */
/*   peak** peaks; */
/* }peaks; */

/* typedef struct{ */
/*   int n_dims; */
/*   double* xs; */
/* }point; */

/* typedef struct{ */
/*   int n_dims; */
/*   int n_Ts; */
/*   int n_walkers; */
/*   double** w_xs; */
/*   double* w_pis; */
/*   int* w_ts; //  */
/*   int* w_LRs; */
/*   int* t_Lws; */
/*   int* t_Rws; */
/* }chain_state; */


/* void old_print_state_walker_order(int n_Ts, int n_d, double** xs);  */
/* void print_states_T_order(chain_state* state); */
/* double pi(peaks* the_peaks, double* x); */
/* void update_x(int n_dims, peaks* the_peaks,  int iw, chain_state* state, double Tinverse); */
/* void T_swap(chain_state* state, double* inverse_Temperatures); */
/* int check_state_consistency(int n_dims, int n_Ts, double* pis, double** xs,  int* Ts, int* LRs, int* Lws, int* Rws); */

gsl_rng* g_rng;
double init_width = 2.0; // points are initialized to uniformly within 2*init_width on a side n_dims-cube centered on origin.
double prop_w = 0.25;

int old_style = 0;

// main:
int main(int argc, char* argv[]){

  // ***** process command line args *****
  int n_updates = (argc >= 2)? atoi(argv[1]) : 100000; // number of updates to do
  int n_dims = (argc >= 3)? atoi(argv[2]) : 3; // number of dimensions
  int n_Ts = (argc >= 4)? atoi(argv[3]) : 3; // number of T levels
  if(n_Ts < 1){ n_Ts = 1; } 
  double Thot = (argc >= 5)? atof(argv[4]) : 2.0;
  double Tratio = (n_Ts > 1)? pow(Thot, 1.0/(n_Ts-1)) : 1.0; // Tratio for geometric Temperatures
  int seed = (argc >= 6)? atoi(argv[5])  : 12345; // rng seed

  double* Temperatures = (double*)malloc(n_Ts*sizeof(double)); // Temperatures[0] is temperature of cold level, etc.
  Temperatures[0] = 1.0; // cold
  for(int i = 1; i < n_Ts; i++){
    Temperatures[i] = Temperatures[i-1] * Tratio;
  }
  double* inverse_Temperatures = (double*)malloc(n_Ts*sizeof(double)); // the inverse temperatures 
  for(int i = 0; i < n_Ts; i++){
    inverse_Temperatures[i] = 1.0/Temperatures[i];
  }

  printf("# n updates: %5i, n_dims: %2i, n_Ts: %2i, Thot: %6.3f, Tratio: %6.3f, rng seed: %8i \n",
         n_updates, n_dims, n_Ts, Thot, Tratio, seed);

  // ***** Set up RNG *****
  gsl_rng_env_setup();
  const gsl_rng_type* rng_type = gsl_rng_default;
  g_rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(g_rng, seed);


  //*********  Set up target distribution *******************
  int n_peaks = 2;
  double A = 1.0; // spacing parameter; peaks are at x[0] = +-A, other coordinates all zero.
  double Sigma = 0.2; // width parameter; stddev if peak is Normal.
  double Hsmall = 0.25; // height of secondary peak(s) rel to primary
  //  double** peak_positions = (double**)malloc(n_peaks*sizeof(double*));
  // double* heights = (double*)malloc(n_peaks*sizeof(double));
  // double* sigmas = (double*)malloc(n_peaks*sizeof(double));

  peaks* the_peaks = (peaks*)malloc(sizeof(peaks));
  the_peaks->n_peaks = n_peaks;
  the_peaks->peaks = (peak**)malloc(n_peaks*sizeof(peak*)); // array of pointers to peaks.
  double peak_height = 1.0;
  for(int i=0; i<n_peaks; i++){
    peak* a_peak = (peak*)malloc(sizeof(peak));
    double* peak_position = (double*)malloc(n_dims*sizeof(double));
    for(int j=0; j<n_dims; j++){
      peak_position[j] = (j == 0)? A*(-1 + 2*i) : 0.0; // peaks at -A, A, 3A, ..
      //  printf("j A peak_position[j]: %i %g %g \n", j, A, peak_position[j]);
    }
    a_peak->n_dims = n_dims;
    a_peak->position = peak_position;
    a_peak->width = Sigma;
    a_peak->height = peak_height;
    the_peaks->peaks[i] = a_peak; 
    
    peak_height *= Hsmall;
  }

  /* for(int j = 0; j < n_peaks; j++){ */
  /*   sigmas[j] = Sigma; */
  /*   heights[j] = (j==0)? 1.0 : Hsmall;  */
  /* } */

  /* for(int i = 0; i < n_peaks; i++){ */
  /*   double* apeak = (double*)malloc(n_dims*sizeof(double)); */
  /*   for(int j = 0; j < n_dims; j++){ */
  /*     apeak[j] = 0; */
  /*   } */
  /*   if(i == 0){ */
  /*     apeak[0] = -1.0*A; */
  /*   }else if(i == 1){ */
  /*     apeak[0] =  A; */
  /*   } */
  /*   peak_positions[i] = apeak;    */
  /* } */

  // ***** Setup initial states of walkers *********
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
   for(int i=0; i < n_walkers; i++){
    double* pt = (double*)malloc(n_dims*sizeof(double));
    for(int j=0; j < n_dims; j++){
      double r = gsl_rng_uniform(g_rng);
      pt[j] = (2.0*r - 1.0)*init_width;
    }
    the_chain_state->w_xs[i] = pt;
    the_chain_state->w_pis[i] = pi(the_peaks, pt);
  
    the_chain_state->w_ts[i] = i / 2; //  n_Ts; // initially e.g. (n_Ts = 3:  Ts = (0,0,1,1,2,2)
    the_chain_state->w_LRs[i] =  i % 2; //  n_Ts; // initially e.g. (n_Ts = 3:  LRs = (0,1,0,1,0,1)    
    if(the_chain_state->w_LRs[i] == 0){
      the_chain_state->t_Lws[the_chain_state->w_ts[i]] = i;
    }else if(the_chain_state->w_LRs[i] == 1){
      the_chain_state->t_Rws[the_chain_state->w_ts[i]] = i;
    }else{
      printf("LRs[i] is neither 0 nor 1. Bye.\n");
      exit(0);
    }
   }


  // ***** Loop through n_updates updates *****
  for(int i = 0; i < n_updates; i++){ // loop through updates
 for(int j = 0; j < n_walkers; j++){ // update positions of walkers
     
          update_x(n_dims, the_peaks, j, the_chain_state, inverse_Temperatures[the_chain_state->w_ts[j]]);
    
    }

    T_swap(the_chain_state, inverse_Temperatures);
    //   check_state_consistency(n_dims, n_Ts, pis, xs, Ts, LRs, Lws, Rws);
    printf("%5i ", i);
    print_states_T_order(the_chain_state);
    printf("\n");
    }

}

// end of main


/* void old_print_states_walker_order(int n_Ts, int n_d, double** xs, int* Ts){ */
/*   int n_w = 2*n_Ts; */
/*   for(int iw=0; iw < n_w; iw++){ */
/*     int it = Ts[iw]; */
/*     printf("%2i %2i  ", iw, it); */
/*     for(int j=0; j < n_d; j++){ */
/*       printf("%5.3f ", xs[iw][j]); */
/*     }printf("   "); */
/*   } // printf("\n"); */
/* } */

/* void print_states_T_order(chain_state* state){ */
  
/*   for(int it=0; it < state->n_Ts; it++){ */
/*     int iw = state->t_Lws[it]; */
/*     printf("%2i %2i  ", iw, it); */
/*     for(int j=0; j < state->n_dims; j++){ */
/*       printf("%5.3f ", state->w_xs[iw][j]); */
/*     }printf("   ");     */
/*   }  */
/*   for(int it=0; it < state->n_Ts; it++){ */
/*     int iw = state->t_Rws[it]; */
/*     printf("%2i %2i  ", iw, it); */
/*     for(int j=0; j < state->n_dims; j++){ */
/*       printf("%5.3f ", state->w_xs[iw][j]); */
/*     }printf("   ");     */
/*   }  */
/*   // printf("\n"); */
/* } */


/* double pi(peaks* the_peaks, double* x){ */
/*   double result = 0.0; */
/*   for(int k = 0; k < the_peaks->n_peaks; k++){    */
/*     double peak_result; */
/*     peak* a_peak = the_peaks->peaks[k]; */
/*     if(1){ */
/*       peak_result = 1.0; */
/*       for(int i = 0; i < a_peak->n_dims; i++){ */
/*         peak_result *= exp(-0.5*pow( (x[i] - a_peak->position[i])/a_peak->width, 2) ); */
/*       } */
/*     }else{  */
/*       double rsq = 0; */
/*       for(int i = 0; i < a_peak->n_dims; i++){ */
/*         //   printf("peak positioni: %g\n", a_peak->position[i]); */
/*         rsq += pow( (x[i] - a_peak->position[i])/a_peak->width, 2 ); */
/*       } */
/*       //     printf("new rsq: %g \n", rsq); */
/*       peak_result = pow((1.0 + rsq), -3); // -0.5*(n_dims + 1)); */
/*       //    printf("#  rsq: %12.6g  peak_result: %12.6g  %12.6g\n", rsq, peak_result, pow(rsq, -6) ); */
/*     } */
/*     result += a_peak->height * peak_result; */
/*   } */
/*   // printf("#  result: %12.6g \n", result); */
/*   return result; */
/* } */


/* void update_x(int n_dims, peaks* the_peaks,  int iw, chain_state* state, double Tinverse){ */
/*   //  printf("top of update_x\n"); */
/*   double pi_x = state->w_pis[iw]; // pi(n_dims, n_peaks, peak_positions, heights, sigmas, x); */
/*   //  printf("after pi_x = pis[iw]. n_dims: %i \n", n_dims); */

/*   double* prop_x = (double*)malloc(n_dims*sizeof(double)); */
/*   //  printf("after alloc prop_x \n"); fflush(stdout); */
/*   for(int i = 0; i < n_dims; i++){ // propose a new position  */
/*     prop_x[i] = state->w_xs[iw][i] + prop_w * (gsl_rng_uniform(g_rng) - 0.5); */
/*   } */

/*   double pi_prop_x = pi(the_peaks, prop_x); */
/*   //  printf("pi_x: %8.5g  %8.5g   %8.5g\n", pi_x, pi_prop_x, Tinverse); */
/*   if( (pi_prop_x > pi_x)  ||  (gsl_rng_uniform(g_rng)  <  pow(pi_prop_x/pi_x, Tinverse)) ){ // Accept proposal */
/*     state->w_pis[iw] = pi_prop_x; */
/*     // return prop_x; */
/*     // double* old_x = xs[iw]; */
/*     free(state->w_xs[iw]);   */
/*     state->w_xs[iw] = prop_x; */
/*     // free(old_x); */
/*   }else{ */
/*     free(prop_x); */
/*   } */
/*   //  printf("bottom of update_x\n"); */
/* } */


/* void T_swap(chain_state* state, double* inverse_Temperatures){ */

/*   for(int i=1; i<state->n_Ts; i++){ // T-swapping among Ls */
/*     double Tc_inverse = // 1.0/Temperatures[i-1]; // 1/T for the colder of two walkers */
/*       inverse_Temperatures[i-1]; */
/*     double Th_inverse = // 1.0/Temperatures[i]; // 1/T for the hotter of the two walkers */
/*       inverse_Temperatures[i]; */
/*     { */
/*       // first the 'left' set of walkers: */
/*       int iw_hot = state->t_Lws[i]; // index of walker with ith temperature */
/*       int iw_cold = state->t_Lws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature */
/*       double pi_cold = state->w_pis[iw_cold]; */
/*       double pi_hot = state->w_pis[iw_hot]; */
/*       if(  */
/*          (pi_hot > pi_cold)  */
/*          ||  */
/*          (gsl_rng_uniform(g_rng)  <  pow(pi_hot/pi_cold, Tc_inverse - Th_inverse))  */
/*           ){ // Accept proposed T-swap */
/*         int tmp = state->t_Lws[i];  */
/*         state->t_Lws[i] = state->t_Lws[i-1]; */
/*         state->t_Lws[i-1] = tmp; */

/*         tmp = state->w_ts[iw_hot]; */
/*         state->w_ts[iw_hot] = state->w_ts[iw_cold]; */
/*         state->w_ts[iw_cold] = tmp; */
/*       } */
/*     } */
/*     { */
/*       // then the 'right' set of walkers: */
/*       int iw_hot = state->t_Rws[i]; // index of walker with ith temperature */
/*       int iw_cold = state->t_Rws[i-1]; // index of walker with (i-1)th, i.e. next colder, temperature */
/*       double pi_cold = state->w_pis[iw_cold]; */
/*       double pi_hot = state->w_pis[iw_hot]; */
/*       if(  */
/*          (pi_hot > pi_cold)  */
/*          ||  */
/*          (gsl_rng_uniform(g_rng)  <  pow(pi_hot/pi_cold, Tc_inverse - Th_inverse))  */
/*           ){ // Accept proposed T-swap */
/*         int tmp = state->t_Rws[i];  */
/*         state->t_Rws[i] = state->t_Rws[i-1]; */
/*         state->t_Rws[i-1] = tmp; */

/*         tmp = state->w_ts[iw_hot]; */
/*         state->w_ts[iw_hot] = state->w_ts[iw_cold]; */
/*         state->w_ts[iw_cold] = tmp; */
/*       } */
/*     } */
/*   } */
/*   //  check_state_consistency(n_dims, n_Ts, pis, xs, Ts, LRs, Lws, Rws); */
/*   //  printf("bottom of T_swap\n"); */
/* } */



/* int check_state_consistency(int n_dims, int n_Ts, double* pis, double** xs,  int* Ts, int* LRs, int* Lws, int* Rws){ */
/*   for(int iw = 0; iw < 2*n_Ts; iw++){ */
/*     int iT = Ts[iw]; // should be Tindex of iw walker. */
/*     int iLR = LRs[iw]; // should be LR index (0 for L, 1 for R) of iw walker */
/*     if(iLR == 0){ // Left */
/*       assert(iw == Lws[iT]); // Lws */
/*     }else if(iLR == 1){ // Right */
/*       assert(iw == Rws[iT]); */
/*     }else{ */
/*       assert(1 == 0); */
/*     } */
/*   } */
/* } */
       
      
      


