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

const Ndim_histogram* g_targp;

// control parameters and default values:
int g_Ngrid_max = 25;
int g_is_ball = 1;
int g_n_modes = 2;
const Target_1dim* g_targ_1d;
Target_peak_1dim* g_peaks;
const gsl_rng_type* g_rng_type;
gsl_rng* g_rng;
FILE* tvd_vs_gen_fstream; 

// main:
int main(int argc, char* argv[]){
  // ********** read control parameters from command line: ***************
  printf("# ");
  for(int i=0; i<argc; i++){
    printf("%s; ", argv[i]);
  }printf("\n");
  
  int n_dimensions = atoi(argv[1]);
  printf("# n_dimensions: %i \n", n_dimensions);

 Target_1dim* targ_1d = (Target_1dim*)malloc(sizeof(Target_1dim));
  targ_1d->n_modes = atoi(argv[2]);
  printf("# n_modes: %i \n", targ_1d->n_modes);
  //Target_peak_1dim* 
  Target_peak_1dim* peaks = (Target_peak_1dim*) malloc(targ_1d->n_modes*sizeof(Target_peak_1dim));
  for(int j=0; j<targ_1d->n_modes; j++){
    peaks[j].position = atof(argv[3 + 3*j]);
    peaks[j].sigma = atof(argv[4 + 3*j]);
    peaks[j].weight = atof(argv[5 + 3*j]);
    printf("# peak %i; position,sigma,weight: %8.5f %8.5f %8.5f \n",j, peaks[j].position, peaks[j].sigma, peaks[j].weight);
  }
  targ_1d->peaks = peaks; 
  normalize_targ_1dim(targ_1d);
  g_targ_1d = (const Target_1dim*)targ_1d;


  int next_arg = 3 + 3*g_n_modes;
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
  g_is_ball = atoi(argv[next_arg]);
  g_Ngrid_max = atoi(argv[next_arg+1]);
  printf("# is ball, g_Ngrid_max: %i %i \n", g_is_ball, g_Ngrid_max);
  // ******************************************************************
  //  normalize_targ_1dim(g_n_modes, g_peaks); //

  // *************** set up RNG *****************
  gsl_rng_env_setup();
  g_rng_type = gsl_rng_default;
  g_rng = gsl_rng_alloc(g_rng_type);

  // *************** open output files *************
  tvd_vs_gen_fstream = fopen("tvd_vs_gen", "w");

  g_Ngrid_max += g_Ngrid_max % 2; // if odd, add one, to make it be even.

  // /*  ********** get bin boundaries ****************
  //  g_bins = construct_binning_spec(g_Ngrid_max, g_peaks);
  const Binning_spec* the_bins = construct_binning_spec(g_Ngrid_max, targ_1d);
  //  printf("ZZZZZ %i \n", (int)sizeof(g_bins));
/* printf("%p \n", g_bins); */
/*   g_bins->n_g_bins = 13; */
/*     printf("ZZZZZ\n"); */
 /* printf("g_bins->n_bins: %i \n", g_bins->n_bins); */
 /* printf("g_bins->bin_edges: %p \n", g_bins->bin_edges); */
 /*  printf("g_bins->x_lo: %g \n", g_bins->x_lo); */
 /*  print_array_of_double(g_bins->n_bins, g_bins->bin_edges); */
  //  exit(0);
  // ********************************
  // */

  time_t t;
  srand( (unsigned) time(&t) );
  //  srand((unsigned)1234567);
 
  printf("before targp\n");
  // printf("g_bins->bin_edges: %p \n", g_bins->bin_edges);
  Ndim_histogram* targp = init_target_distribution(n_dimensions, g_Ngrid_max, 1, the_bins);
 printf("after targp\n");
  g_targp = targp;
  // printf("after targp \n");
  double tvd_sum = 0.0;
  double tvdsq_sum = 0.0;
  int n_jumps = 0;
  int Naccept_sum = 0;
  int Nreject_sum = 0;
  //sum_dsq = 0.0;

  for(int j=0; j<Nreps; j++){
    State** states = (State**)malloc(n_temperatures*sizeof(State*));
    for(int i=0; i<n_temperatures; i++){
      states[i] = construct_state(n_dimensions, targ_1d);
    }
    printf("after construct states \n");
    Multi_T_chain* multi_T_chain = construct_multi_T_chain(n_temperatures, temperatures, proposals, states, the_bins);
    // printf("after multi_T_chain \n");
    // ********** do burn-in ***********
    //	  State* the_state = construct_state(n_dimensions, Ngrid_max);
    for(int n=0; n<=burn_in_steps; n++){
      multi_T_chain_within_T_mcmc_step(multi_T_chain);
    } 
    // printf("after burn-in\n");
    // ********** do post-burn-in *********
    //  g_n_accept = 0;
    //  g_n_reject = 0;
    //  mcmc_out_hist = construct_ndim_histogram(n_dimensions, Ngrid_max);
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
    //   Naccept_sum += g_n_accept;
    //  Nreject_sum += g_n_reject;
    //  free_state(the_state);
    //  mcmc_out_hist = free_ndim_histogram(mcmc_out_hist);
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



/* int i_from_x(double x){ // xmin, xmax, Ngrid_max are global */
/*   int gridpoint; */
/*   if(x <= Xmin){ */
/*     gridpoint = 0; */
/*   }else if(x >= Xmax){ */
/*     gridpoint = Ngrid_max; */
/*   }else{ */
/*     gridpoint = (int)(Ngrid_max * (x-Xmin)/(Xmax-Xmin) + 0.5); */
/*   } */
/*   return gridpoint; */
/* } */
/* double* x_array_from_i_array(int Ndim, int* i_array){ */
/*   double* x_array = (double*)malloc(Ndim*sizeof(double)); */
/*   for(int i=0; i<Ndim; i++){ */
/*     x_array[i] = x_from_i(i_array[i]); */
/*   } */
/*   return x_array; */
/* } */
/* double x_from_i(int i){ */
/*   double x; */
/*   return Xmin + (double)i/Ngrid_max * (Xmax - Xmin); */
/* } */

double F(const Target_1dim* targ_1d, int n_dim, double* x_array){
  double result = 1.0;
  for(int i=0; i<n_dim; i++){
    double x = x_array[i];
    result *= f_1dim(targ_1d, x);
  }
  //  printf("x p: %12.6g %12.6g \n", x_array[0], result);
  return result;
}    

double f_1dim(const Target_1dim* targ_1d, double x){
  int n_modes = targ_1d->n_modes;
  Target_peak_1dim* peaks = targ_1d->peaks;
  double f = 0.0;
  for(int i=0; i<n_modes; i++){
    double X = (x-peaks[i].position)/peaks[i].sigma;
    f += exp(-0.5 * X*X) * peaks[i].weight/peaks[i].sigma;
    //   printf("i, x, X, position, sigma: %i %8.5f %8.5f %8.5f %8.5f \n", i, x, X, g_peaks[i].position, g_peaks[i].sigma);
  }
  return f*ONEOVERSQRT2PI;
}


double drand(void){
  return (double)rand() / ((double)RAND_MAX + 1);
}

Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize, const Binning_spec* bins){
  Ndim_histogram* targ_bin_probs = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  targ_bin_probs->Ndim = Ndim;
  targ_bin_probs->Ngrid_max = Ngrid_max;
  targ_bin_probs->weights = construct_ndim_array_of_double(Ndim, Ngrid_max, 0.0);
    printf("before set_...\n");
  targ_bin_probs->total_weight = 
    set_ndim_array_of_double_to_target(targ_bin_probs->weights, 1.0, g_targ_1d->peaks, bins);
    printf("after set_...\n");
  //  print_ndim_array_of_double(targ_bin_probs->weights); exit(0);
  if(normalize){
    normalize_ndim_histogram(targ_bin_probs);
  }
  return targ_bin_probs;
}

double total_variation_distance(const Ndim_histogram* targp, const Ndim_histogram* mcmc_out){
  // both distributions must be normalized ahead of time
  Ndim_array_of_double* a1 = targp->weights;
  Ndim_array_of_double* a2 = mcmc_out->weights;
  /* printf("# targp\n#"); */
  /* print_ndim_array_of_double(a1); */
  /* printf("# hist\n#"); */
  /* print_ndim_array_of_double(a2); */
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


double draw_1dim(const Target_1dim* targ_1d){
  double wsum = 0.0;
  int n_modes = targ_1d->n_modes;
  Target_peak_1dim* peaks = targ_1d->peaks;
  for(int i=0; i<n_modes; i++){
    wsum += peaks[i].weight;
  }
  double random_number = drand();
  double w = 0.0;
  for(int i=0; i<n_modes; i++){
    peaks[i].weight /= wsum;
    w += peaks[i].weight;
    if( random_number < w ){
      return gsl_ran_gaussian(g_rng, peaks[i].sigma) + peaks[i].position;
    }
  }
  return gsl_ran_gaussian(g_rng, peaks[n_modes-1].sigma) + peaks[n_modes-1].position;
}

double* draw_ndim(int n_dim, const Target_1dim* targ_1d){
  double* xs = (double*)malloc(n_dim*sizeof(double));
  for(int i=0; i<n_dim; i++){
    xs[i] = draw_1dim(targ_1d);
  }
  return xs;
}

/* double draw_1dim(int n_mode, Target_peak_1dim* peaks){ */
/*   double p = 0.0; */
/*   double drnd = drand(); */
/*   for(int i=0; i<n_modes; i++){ */
/*     p += peaks[i].weight; */
/*     if(drnd < p){ // draw from i_th peak */
/*       return gsl_ran_gaussian(the_rng, peaks[i].sigma) + peaks[i].position; */
/*     } */
/*   } */
/*   // shouldn't get here if normalization is right, but in case ... */
/*   return gsl_ran_gaussian(the_rng, peaks[n_modes-1].sigma) + peaks[n_modes-1].position; */
/* } */


/* int* ipropose(int Ndim, int* index_array, Proposal* prop, int is_ball){ */
/*   int* prop_index_array = (int*)malloc(Ndim*sizeof(int));  */
/*   int Width = (drand() < prop->p1)? prop->W1 : prop->W2; */
/*   while(1){ */
/*     int sum_sq = 0; */
/*     for(int i = 0; i<Ndim; i++){ */
/*       int delta_i = (int)((2*Width+1)*drand()) - Width; */
/*       prop_index_array[i] = index_array[i] + delta_i; */
/*       sum_sq += delta_i*delta_i; */
/*     } */
/*     if((sum_sq > 0) && (!is_ball  || sum_sq <= Width*Width)) break; */
/*   } */
/*   return prop_index_array; */
/* } */
