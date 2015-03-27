#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "mcmc_toy_structs.h"
#include "mcmc_toy.h"

int n_accept = 0;
int n_reject = 0;

// 
double sum_dsq = 0;
Ndim_histogram* mcmc_out_hist;
Ndim_histogram* targp;

// control parameters and default values:
int Ngrid_max = 25;
double sigma = 0.02;
int is_ball = 1;
int n_modes = 2;
double Xmin = -0.5;
double Xmax = 1.5;
Targ_peak_1dim* peaks;
const gsl_rng_type* rng_type;
gsl_rng* the_rng;
Binning_spec* bins;
FILE* tvd_vs_gen_fstream; 

// main:
int main(int argc, char* argv[]){
  // ********** read control parameters from command line: ***************
  printf("# ");
  for(int iarg = 0; iarg<argc; iarg++){
    printf("%s; ", argv[iarg]);
  }printf("\n");
  
  int n_dimensions = atoi(argv[1]);
  printf("# n_dimensions: %i \n", n_dimensions);
  n_modes = atoi(argv[2]);
  printf("# n_modes: %i \n", n_modes);
  //Targ_peak_1dim* 
  peaks = (Targ_peak_1dim*) malloc(n_modes*sizeof(Targ_peak_1dim));
  for(int j=0; j<n_modes; j++){
    peaks[j].position = atof(argv[3 + 3*j]);
    peaks[j].sigma = atof(argv[4 + 3*j]);
    peaks[j].weight = atof(argv[5 + 3*j]);
    printf("# peak %i; position,sigma,weight: %8.5f %8.5f %8.5f \n",j, peaks[j].position, peaks[j].sigma, peaks[j].weight);
  }
  int next_arg = 3 + 3*n_modes;
  int burn_in_steps = atoi(argv[next_arg]);
  int mcmc_steps = atoi(argv[next_arg+1]);
  int Nreps = atoi(argv[next_arg+2]);
  printf("# burn-in, mcmc steps, n reps: %i %i %i \n", burn_in_steps, mcmc_steps, Nreps);
  int n_temperatures = atoi(argv[next_arg+3]);
  printf("# n temperatures: %i \n", n_temperatures);
  next_arg += 4;
  double* temperatures = (double*)malloc(n_temperatures*sizeof(double));
  Proposal* proposals = (Proposal*)malloc(n_temperatures*sizeof(Proposal));
  for(int j=0; j<n_temperatures; j++){
    temperatures[j] = atof(argv[next_arg + 4*j]);
    proposals[j].W1 = atof(argv[next_arg+1 + 4*j]);
    proposals[j].W2 = atof(argv[next_arg+2 + 4*j]);
    proposals[j].p1 = atof(argv[next_arg+3 + 4*j]);
    printf("#  T,W1,W2,p1:  %8.5f %8.5f %8.5f %8.5f \n", temperatures[j], proposals[j].W1, proposals[j].W2, proposals[j].p1);
  }
  next_arg += 4*n_temperatures;
  is_ball = atoi(argv[next_arg]);
  Ngrid_max = atoi(argv[next_arg+1]);
  printf("# is ball, Ngrid_max: %i %i \n", is_ball, Ngrid_max);
  // ******************************************************************
  normalize_targ_1dim(n_modes, peaks); //

  // *************** set up RNG *****************
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  the_rng = gsl_rng_alloc(rng_type);

  // *************** open output files *************
  tvd_vs_gen_fstream = fopen("tvd_vs_gen", "w");

  Ngrid_max += Ngrid_max % 2; // if odd, add one, to make it be even.
  // /*  ********** get bin boundaries ****************
  bins = construct_binning_spec(Ngrid_max, peaks);
  //  printf("ZZZZZ %i \n", (int)sizeof(bins));
/* printf("%p \n", bins); */
/*   bins->n_bins = 13; */
/*     printf("ZZZZZ\n"); */
//  printf("bins->n_bins: %i \n", bins->n_bins);
//  printf("bins->bin_edges: %p \n", bins->bin_edges);
  //printf("bins->x_lo: %g \n", bins->x_lo);
 
  //  exit(0);
  // ********************************
  // */

  time_t t;
  srand( (unsigned) time(&t) );
  //  srand((unsigned)1234567);
 
  // printf("before targp\n");
  // printf("bins->bin_edges: %p \n", bins->bin_edges);
  targp = init_target_distribution(n_dimensions, Ngrid_max, 1);
  // printf("after targp \n");
  double tvd_sum = 0.0;
  double tvdsq_sum = 0.0;
  int n_jumps = 0;
  int Naccept_sum = 0;
  int Nreject_sum = 0;
  sum_dsq = 0.0;

  for(int j=0; j<Nreps; j++){
    State** states = (State**)malloc(n_temperatures*sizeof(State*));
    for(int i=0; i<n_temperatures; i++){
      states[i] = construct_state(n_dimensions, Ngrid_max);
    }
    // printf("after construct states \n");
    Multi_T_chain* multi_T_chain = construct_multi_T_chain(n_temperatures, temperatures, proposals, states);
    // printf("after multi_T_chain \n");
    // ********** do burn-in ***********
    //	  State* the_state = construct_state(n_dimensions, Ngrid_max);
    for(int n=0; n<=burn_in_steps; n++){
      multi_T_chain_within_T_mcmc_step(multi_T_chain);
    } 
    // printf("after burn-in\n");
    // ********** do post-burn-in *********
    n_accept = 0;
    n_reject = 0;
    mcmc_out_hist = construct_ndim_histogram(n_dimensions, Ngrid_max);
    // printf("after construct ndim hist\n");
    double tvd = 0.0;
    for(int n=1; n<=mcmc_steps; n++){   
      // take a mcmc step
      //	  n_jumps += mcmc_step(the_state, &proposal);
      //    printf("before mcmc_step %i \n", n);
      multi_T_chain_within_T_mcmc_step(multi_T_chain);
	
      /* if(OUTPUT_TVD_VS_N){ */
      /* 	if(n%10000 == 0){ */
      /* 	  Ndim_histogram* hist_copy = construct_copy_ndim_histogram(mcmc_out_hist); */
      /* 	  normalize_ndim_histogram(hist_copy); */
      /* 	  tvd = total_variation_distance(targp, hist_copy); */
      /* 	  printf("#%8i %12.6g\n", n, tvd); */
      /* 	  free_ndim_histogram(hist_copy); */
      /* 	} */
      /* } */
    }

    /* Ndim_histogram* hist_copy = construct_copy_ndim_histogram(mcmc_out_hist); */
    /* normalize_ndim_histogram(hist_copy); */
    /* tvd = total_variation_distance(targp, hist_copy); */
    /* printf("#%8i %12.6g\n", mcmc_steps, tvd); */
    /* free_ndim_histogram(hist_copy); */
    //	  printf("\n\n");
    tvd_sum += tvd;
    tvdsq_sum += tvd*tvd;
    Naccept_sum += n_accept;
    Nreject_sum += n_reject;
    //  free_state(the_state);
    mcmc_out_hist = free_ndim_histogram(mcmc_out_hist);
    //	printf("bottom of Nrep loop. j: %i \n", j);
  } // end loop over reps
  //sum_dsq /= Nreps*mcmc_steps;
  double tvd_avg = tvd_sum/Nreps;
  double tvd_variance = tvdsq_sum/Nreps - tvd_avg*tvd_avg;
  double tvd_stderr = sqrt(tvd_variance)/sqrt((double)Nreps);
  //      printf("%4i  %8i  %10.6f +- %10.6f  %10.6f  %10.6f \n", k, n_jumps, tvd_avg, tvd_stderr, (double)Naccept_sum/(Naccept_sum+Nreject_sum), sqrt(sum_dsq/(Nreps*mcmc_steps)));
  //	printf("\n\n");
  fflush(stdout);
  free_ndim_histogram(targp);
} // end of main

int* ipropose(int Ndim, int* index_array, Proposal* prop, int is_ball){
  int* prop_index_array = (int*)malloc(Ndim*sizeof(int)); 
  int Width = (drand() < prop->p1)? prop->W1 : prop->W2;
  while(1){
    int sum_sq = 0;
    for(int i = 0; i<Ndim; i++){
      int delta_i = (int)((2*Width+1)*drand()) - Width;
      prop_index_array[i] = index_array[i] + delta_i;
      sum_sq += delta_i*delta_i;
    }
    if((sum_sq > 0) && (!is_ball  || sum_sq <= Width*Width)) break;
  }
  return prop_index_array;
}

// propose uniformly from cube or ball centered on present point
double* propose(int Ndim, double* x_array, Proposal* prop, int is_ball){
  double* prop_x_array = (double*)malloc(Ndim*sizeof(double));
  //  printf("porposal:: %8.5f %8.5f %8.5f \n", prop->W1, prop->W2, prop->p1);
  double Width = (drand() < prop->p1)? (double)prop->W1 : (double)prop->W2;
  while(1){
    double sum_sq = 0;
    for(int i = 0; i<Ndim; i++){
      double delta_x = (double)(2*Width*drand()) - Width;
      //   printf("width: %12.5g  delta_x: %12.5g \n", Width, delta_x);
      prop_x_array[i] = x_array[i] + delta_x;
      sum_sq += delta_x*delta_x;
    }
    if(!is_ball  || sum_sq <= Width*Width) break;
  }
  return prop_x_array;
}

int* i_array_from_x_array(int Ndim, double* x_array){
  int* i_array = (int*)malloc(Ndim*sizeof(int));
  for(int i=0; i<Ndim; i++){
    //   printf("i, x[i]: %i %g \n", i, x_array[i]);
    i_array[i] = x_to_bin(bins, x_array[i]); // i_from_x(x_array[i]);
    //  printf("AFTER\n");
  }
  return i_array;
}

int i_from_x(double x){ // xmin, xmax, Ngrid_max are global
  int gridpoint;
  if(x <= Xmin){
    gridpoint = 0;
  }else if(x >= Xmax){
    gridpoint = Ngrid_max;
  }else{
    gridpoint = (int)(Ngrid_max * (x-Xmin)/(Xmax-Xmin) + 0.5);
  }
  return gridpoint;
}
double* x_array_from_i_array(int Ndim, int* i_array){
  double* x_array = (double*)malloc(Ndim*sizeof(double));
  for(int i=0; i<Ndim; i++){
    x_array[i] = x_from_i(i_array[i]);
  }
  return x_array;
}
double x_from_i(int i){
  double x;
  return Xmin + (double)i/Ngrid_max * (Xmax - Xmin);
}

double F(int Ndim, double* x_array){
  double result = 1.0;
  for(int i=0; i<Ndim; i++){
    double x = x_array[i];
    result *= f_1dim(x);
  }
  //  printf("x p: %12.6g %12.6g \n", x_array[0], result);
  return result;
}    



double f_1dim(double x){
  double f = 0.0;
  for(int i=0; i<n_modes; i++){
    double X = (x-peaks[i].position)/peaks[i].sigma;
    f += exp(-0.5 * X*X) * peaks[i].weight/peaks[i].sigma;
    //   printf("i, x, X, position, sigma: %i %8.5f %8.5f %8.5f %8.5f \n", i, x, X, peaks[i].position, peaks[i].sigma);
  }
  return f;
}

double draw_1dim(void){
  double p = peaks[0].weight;
  double drnd = drand();
  for(int i=0; i<n_modes; i++){
    p += peaks[i].weight;
    if(drnd < p){ // draw from i_th peak
      return gsl_ran_gaussian(the_rng, peaks[i].sigma) + peaks[i].position;
    }
  }
  // shouldn't get here if normalization is right, but in case ...
  return gsl_ran_gaussian(the_rng, peaks[n_modes-1].sigma) + peaks[n_modes-1].position;
}

double drand(void){
  return (double)rand() / ((double)RAND_MAX + 1);
}

Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize){
  Ndim_histogram* targ_bin_probs = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  targ_bin_probs->Ndim = Ndim;
  targ_bin_probs->Ngrid_max = Ngrid_max;
  targ_bin_probs->weights = construct_ndim_array_of_double(Ndim, Ngrid_max, 0.0);
  //  printf("before set_...\n");
  targ_bin_probs->total_weight = 
    set_ndim_array_of_double_to_target(targ_bin_probs->weights, 1.0, peaks, bins);
  //  printf("after set_...\n");
  //  print_ndim_array_of_double(targ_bin_probs->weights); exit(0);
  if(normalize){
    normalize_ndim_histogram(targ_bin_probs);
  }
  return targ_bin_probs;
}

double total_variation_distance(Ndim_histogram* targp, Ndim_histogram* mcmc_out){
  // both distributions must be normalized ahead of time
  Ndim_array_of_double* a1 = targp->weights;
  Ndim_array_of_double* a2 = mcmc_out->weights;
  printf("# targp\n#");
  print_ndim_array_of_double(a1);
  printf("# hist\n#");
  print_ndim_array_of_double(a2);
   double tvd = 0.5 * sum_abs_difference_ndim_arrays_of_double(targp->weights, mcmc_out->weights); 
  return tvd;
}

void print_array_of_int(int Nsize, int* array){
  for(int i=0; i<Nsize; i++){
    printf("%6i ", array[i]);
  }// printf("\n");
}

void print_array_of_double(int Nsize, double* array){
  for(int i=0; i<Nsize; i++){
    printf("%8.4g ", array[i]);
  }// printf("\n");
}

int x_to_bin(Binning_spec* bin_spec, double x){
  //  printf("top of x_to_bin. x: %g\n", x);
  //  printf("n_bins, x_lo, x_hi: %i  %g  %g\n", bin_spec->n_bins, bin_spec->x_lo, bin_spec->x_hi);
  int r1;
  int n_bins = bin_spec->n_bins;
  double* edges = bin_spec->bin_edges;
  if(1){
    if(x < edges[0]){ r1 = -1; } //return -1; }
    else if(x < edges[n_bins]){
      for(int i=1; i<=n_bins; i++){
  	//   printf("i, x, bin_edges: %i %g %g \n", i, x, bin_spec->bin_edges[i]);
  	if(x < bin_spec->bin_edges[i]){ r1 = i-1; break; } //return i; }
      }
    }
    //   else if(x <= bin_spec->x_hi){ r1 = bin_spec->n_bins-1; } // return bin_spec->n_bins-1; }
    else if(x >= edges[n_bins]){ r1 = -2; } //return -2; }
    else{
      assert(0); // shouldn't get here
    }
  }
  //  assert(0);
  if(1){
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
  //  printf("in x_to_bin. x, r1, r2: %g %i %i \n", x, r1, min_bin); 
  assert (r1 == min_bin);
}
  return r1; //min_bin;
}
