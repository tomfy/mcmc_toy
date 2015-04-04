#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "mcmc_toy_structs.h"
#include "mcmc_toy.h"

const Target_1dim* g_targ_1d;
const Ndim_histogram* g_targprobs_ndim;
const Ndim_histogram* g_targprobs_1dim;
const Ndim_histogram* g_targprobs_orthants;
const Ndim_histogram* g_targprobs_one_orthant;

gsl_rng* g_rng;

FILE* g_tvd_vs_gen_fstream; 

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


  int next_arg = 3 + 3*targ_1d->n_modes;
  int burn_in_steps = atoi(argv[next_arg]);
  int mcmc_steps = atoi(argv[next_arg+1]);
  int n_replicates = atoi(argv[next_arg+2]);
  printf("# burn-in, mcmc steps, n reps: %i %i %i \n", burn_in_steps, mcmc_steps, n_replicates);
  int n_temperatures = atoi(argv[next_arg+3]);
  printf("# n temperatures: %i \n", n_temperatures);
  next_arg += 4;
  double* temperatures = (double*)malloc(n_temperatures*sizeof(double));
  Proposal* proposals = (Proposal*)malloc(n_temperatures*sizeof(Proposal));
  for(int j=0; j<n_temperatures; j++){
    temperatures[j] = atof(argv[next_arg + 5*j]);
    proposals[j].shape = argv[next_arg+1 + 5*j];
    //   printf("%s \n", proposals[j].shape);
    proposals[j].W1 = atof(argv[next_arg+2 + 5*j]);
    proposals[j].W2 = atof(argv[next_arg+3 + 5*j]);
    proposals[j].p1 = atof(argv[next_arg+4 + 5*j]);
    printf("#  T, Proposal(shape,W1,W2,p1):  %8.5f %12s %8.5f %8.5f %8.5f \n", 
	   temperatures[j], proposals[j].shape, proposals[j].W1, proposals[j].W2, proposals[j].p1);
  }
  next_arg += 5*n_temperatures;
  int n_bins = atoi(argv[next_arg]);
  n_bins += (n_bins % 2); // to make it even
  printf("#  n_bins: %i \n", n_bins);
  // ******************************************************************
  //  normalize_targ_1dim(g_n_modes, g_peaks); //

  // *************** set up RNG *****************

  gsl_rng_env_setup();
  const gsl_rng_type* rng_type = gsl_rng_default;
  g_rng = gsl_rng_alloc(rng_type);

  // *************** open output files *************
  g_tvd_vs_gen_fstream = fopen("tvd_vs_gen", "w");

  //g_Ngrid_max += g_Ngrid_max % 2; // if odd, add one, to make it be even.

  // /*  ********** get bin boundaries ****************
  /* const Binning_spec* the_bins = //construct_binning_spec(n_bins, targ_1d); */
  /*    construct_binning_spec(n_bins, targ_1d, -INFINITY, INFINITY); // Use for Ndim histograms */
  /* const Binning_spec* two_bins = construct_binning_spec(2, targ_1d, -INFINITY, INFINITY); // bin for each orthant */
  /* const Binning_spec* one_orthant_bins = construct_binning_spec(n_bins/2, targ_1d, 0, INFINITY);  */
  /* const Binning_spec* one_dim_bins = construct_binning_spec(2*n_bins, targ_1d, -INFINITY, INFINITY); */

  Binning_spec_set the_bin_set;
  the_bin_set.ndim_bins =  construct_binning_spec(n_bins, targ_1d, -INFINITY, INFINITY); // Use for Ndim histograms
  the_bin_set.onedim_bins = construct_binning_spec(2*n_bins, targ_1d, -INFINITY, INFINITY);
  the_bin_set.orthants_bins = construct_binning_spec(2, targ_1d, -INFINITY, INFINITY); // bin for each orthant
  the_bin_set.positive_bins =  construct_binning_spec(n_bins/2, targ_1d, 0, INFINITY); 

  Ndim_histogram* targprobs_ndim = init_target_distribution(n_dimensions, n_bins, 1, the_bin_set.ndim_bins);
  g_targprobs_ndim = targprobs_ndim;
  Ndim_histogram* targprobs_1dim = init_target_distribution(1, 2*n_bins, 1, the_bin_set.onedim_bins);
  g_targprobs_1dim = targprobs_1dim;
  Ndim_histogram* targprobs_orthants = init_target_distribution(n_dimensions, 2, 1, the_bin_set.orthants_bins);
  g_targprobs_orthants = targprobs_orthants;
  Ndim_histogram* targprobs_one_orthant = init_target_distribution(n_dimensions, n_bins/2, 1, the_bin_set.positive_bins);
  g_targprobs_one_orthant = targprobs_one_orthant;
  printf("XXX\n");
  double tvd_sum = 0.0;
  double tvdsq_sum = 0.0;

  for(int j=0; j<n_replicates; j++){
    State** states = (State**)malloc(n_temperatures*sizeof(State*));
    for(int i=0; i<n_temperatures; i++){
      states[i] = construct_state(n_dimensions, targ_1d);
    }
    Multi_T_chain* multi_T_chain = construct_multi_T_chain(n_temperatures, temperatures, proposals, states, &the_bin_set);
    printf("YYY\n");
  // ********** do burn-in ***********
    for(int n=0; n<=burn_in_steps; n++){
      multi_T_chain_within_T_mcmc_step(multi_T_chain);
    } 

    // ********** do post-burn-in *********
    double tvd = 0.0;
    for(int n=1; n<=mcmc_steps; n++){   
      // take a mcmc step
      multi_T_chain_within_T_mcmc_step(multi_T_chain);   
    }

    tvd_sum += tvd;
    tvdsq_sum += tvd*tvd;
    //  printf("before free_multi...\n");
    free_multi_T_chain(multi_T_chain);
  } // end loop over reps
  // printf("after reps loop...\n");
  double tvd_avg = tvd_sum/n_replicates;
  double tvd_variance = tvdsq_sum/n_replicates - tvd_avg*tvd_avg;
  double tvd_stderr = sqrt(tvd_variance)/sqrt((double)n_replicates);
  //      printf("%4i  %8i  %10.6f +- %10.6f  %10.6f  %10.6f \n", k, n_jumps, tvd_avg, tvd_stderr, (double)Naccept_sum/(Naccept_sum+Nreject_sum), sqrt(sum_dsq/(n_replicates*mcmc_steps)));
  //	printf("\n\n");
  fflush(stdout);
  free_ndim_histogram(targprobs_ndim);
 free_ndim_histogram(targprobs_1dim);
 free_ndim_histogram(targprobs_orthants);
 free_ndim_histogram(targprobs_one_orthant);
} // end of main


// propose uniformly from cube or ball centered on present point
double* propose(int n_dim, double* x_array, Proposal* prop){
  double* prop_x_array = (double*)malloc(n_dim*sizeof(double));
  //  printf("porposal:: %8.5f %8.5f %8.5f \n", prop->W1, prop->W2, prop->p1);
  double Width = (drand() < prop->p1)? (double)prop->W1 : (double)prop->W2;
  if(!strcmp(prop->shape, "gaussian")){
    for(int i=0; i<n_dim; i++){
      prop_x_array[i] = x_array[i] + gsl_ran_gaussian(g_rng, Width);
    }
  }else{
    while(1){
      double sum_sq = 0;
      for(int i = 0; i<n_dim; i++){
	double delta_x = (double)(2*Width*drand()) - Width;
	//   printf("width: %12.5g  delta_x: %12.5g \n", Width, delta_x);
	prop_x_array[i] = x_array[i] + delta_x;
	sum_sq += delta_x*delta_x;
      }
      if(!strcmp(prop->shape, "cube")  || sum_sq <= Width*Width) break;
    }
  }
  return prop_x_array;
}

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
  return f*ONEOVERSQRT2PI; // to normalize such that integral from -inf to +inf is sum(weights)
}

double integral_f_1dim(const Target_1dim* targ_1d, double x, double y){ // integral of f_1dim from x to y
  Target_peak_1dim* peaks = targ_1d->peaks;
  double integral = 0.0;
  for(int i=0; i<targ_1d->n_modes; i++){
    integral += peaks[i].weight*( gsl_cdf_gaussian_P(y-peaks[i].position, peaks[i].sigma) 
				  - gsl_cdf_gaussian_P(x-peaks[i].position, peaks[i].sigma) );
  }
  // printf("A %g %g %g \n", x, y, integral);
  return integral;
}

double find_bin_upper_edge(const Target_1dim* targ_1d, double xlo, double Q){
  // returns x, the upper limit of integration such that integral from xlo to x of f_1dim
  // is Q 
  double epsilon = 1e-8;
  double ll = -20.0;
  double hh = 20.0;
  double xlb = xlo; // current lower bound on result
  double xub = INFINITY; // current upper bound on result
  double xtry =  (xlo > ll)? xlo + Q/f_1dim(targ_1d, xlb): ll;
  //  printf("xtry: %g \n", xtry);
  double q_xtry = integral_f_1dim(targ_1d, xlo, xtry);
  //      printf("B %g  %g %g \n", xtry, q_xtry, Q);

	while(fabs(q_xtry - Q) > epsilon  || (xub - xlb) > epsilon){
    //  q_xtry = integral_f_1dim(targ_1d, xlo, xtry);
    //    printf("B %g  %g %g \n", xtry, q_xtry, Q);
    if(q_xtry > Q){
      xub = xtry;
    }else{
      xlb = xtry;
    }
    xtry = (xub < hh)? 0.5*(xlb + xub) : 0.5*(xlb + hh);
    q_xtry = integral_f_1dim(targ_1d, xlo, xtry);
    //    printf("C %g  %g %g \n", xlb, xub, xtry);
  }
	//  printf("xlo xlb xub xtry: %g %g %g %g \n", xlo, xlb, xub, xtry);
  return xtry;
}

double drand(void){
  //  return (double)rand() / ((double)RAND_MAX + 1);
    return gsl_rng_uniform(g_rng);
}

Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize, const Binning_spec* bins){
  Ndim_histogram* targ_bin_probs = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  targ_bin_probs->Ndim = Ndim;
  targ_bin_probs->Ngrid_max = Ngrid_max;
  targ_bin_probs->weights = construct_ndim_array_of_double(Ndim, Ngrid_max, 0.0);
  targ_bin_probs->total_weight = 
    set_ndim_array_of_double_to_target(targ_bin_probs->weights, 1.0, g_targ_1d->peaks, bins);
  //  print_ndim_array_of_double(targ_bin_probs->weights); exit(0);
  if(normalize){
    normalize_ndim_histogram(targ_bin_probs);
  }
  return targ_bin_probs;
}

/* double total_variation_distance(const Ndim_histogram* targp, const Ndim_histogram* mcmc_out){ */
/*   // both distributions must be normalized ahead of time */
/*    double tvd = 0.5 * sum_abs_difference_ndim_arrays_of_double(targp->weights, mcmc_out->weights);  */
/*   return tvd; */
/* } */

double total_variation_distance(const Ndim_histogram* targprobs, const Ndim_histogram* hist){
   Ndim_histogram* hist_copy = construct_copy_ndim_histogram(hist);
  normalize_ndim_histogram(hist_copy);
  double tvd =  0.5 * sum_abs_difference_ndim_arrays_of_double(targprobs->weights, hist_copy->weights); 
  // total_variation_distance(targprobs, hist_copy);
  //  fprintf(g_tvd_vs_gen_fstream, "%12.6g ", tvd);
  free_ndim_histogram(hist_copy);
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



// ***** unused/obsolete stuff ****

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
