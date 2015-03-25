#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "mcmc_toy_structs.h"
#include "mcmc_toy.h"

int n_accept = 0;
int n_reject = 0;

// 
double sum_dsq = 0;
Ndim_histogram* mcmc_out_hist;

// control parameters and default values:
int Ngrid_max = 25;
double sigma = 0.02;
int is_ball = 1;
int n_modes = 2;

#define OUTPUT_TVD_VS_N (1)
// main:
int main(int argc, char* argv[]){
  for(int iarg = 0; iarg<argc; iarg++){
    printf("%s; ", argv[iarg]);
  }printf("\n");
  time_t t;
  srand( (unsigned) time(&t) );
  int n_dimensions = 3; 
  // int Ngrid_max = 25;
  if(Ngrid_max % 2 == 0){ Ngrid_max++; } // if even, add one, to make it be odd.
  int mcmc_steps  = 100000;
  int Nreps = 1;
  int burn_in_steps = 200;
  double sigma_in_steps = sigma*Ngrid_max;
  int prop_width_1 = -1;
  int prop_width_2 = 3;
  double prop_width_1_prob = 0.35;
  int n_temperatures = 2;
  //  double small_width_prob = 0.85;
  printf("# sigma: %8.4f  is_ball: %4i n_modes: %4i \n", sigma, is_ball, n_modes);
  printf("# n_dimensions: %4i Ngrid_max: %4i \n", n_dimensions, Ngrid_max);
  printf("# proposal W1, W2, p1: %4i %4i %8.4f \n", prop_width_1, prop_width_2, prop_width_1_prob);
  printf("# sigma*Ngrid_max: %12.8f \n", sigma*Ngrid_max);
  printf("# mcmc_steps %i Nreps: %i \n", mcmc_steps, Nreps);
  

  Ndim_histogram* targp = init_target_distribution(n_dimensions, Ngrid_max, 1);

  //int k = 8;
  for(int k=18; k<=18; k += 5) // = (int)(1.1*(k+1)) ) // loop over different proposals
    {    
      prop_width_1 = k;
      double tvd_sum = 0.0;
      double tvdsq_sum = 0.0;
      int n_jumps = 0;
      int Naccept_sum = 0;
      int Nreject_sum = 0;
      sum_dsq = 0.0;
      for(int j=0; j<Nreps; j++){
	Proposal T0_proposal = {k, 3, 0.35};
	Proposal T1_proposal = {k, 3, 0.35}; 
	Proposal proposals[2] = {T0_proposal, T1_proposal};
	double temperatures[2] = {1.0, 2.0};
	State* states[2] = {construct_state(n_dimensions, Ngrid_max), construct_state(n_dimensions, Ngrid_max)};
	//	for(int iT = 0; iT < n_temperatures; iT++){
	// Proposal proposal = // {3, 3, 0.5};  {prop_width_1, prop_width_2, prop_width_1_prob}; //, small_width, (1.0 - small_width_prob)};
	// 	printf("# proposal: %i %i %8.4f \n", proposal.W1, proposal.W2, proposal.p1);
	//	printf("before construct_multi_T_chain ...\n");
	Multi_T_chain* multi_T_chain = construct_multi_T_chain(n_temperatures, temperatures, proposals, states);
	//	printf("after construct_multi_T_chain ...j: %i \n", j);
	// ********** do burn-in ***********
	State* the_state = construct_state(n_dimensions, Ngrid_max);
	for(int n=0; n<=burn_in_steps; n++){
	  //	  printf("before mcmc_step %i j (replicate:) %i \n", n, j);
	  multi_T_chain_within_T_mcmc_step(multi_T_chain);
	  //printf("after mcmc_step %i \n", n);
	} 

	// ********** do post-burn-in *********
	n_accept = 0;
	n_reject = 0;
	mcmc_out_hist = construct_histogram_ndim(n_dimensions, Ngrid_max);
     
	double tvd = 0.0;
	for(int n=1; n<=mcmc_steps; n++){   
	  // take a mcmc step
	  //	  n_jumps += mcmc_step(the_state, &proposal);
	  //	  printf("before mcmc_step %i \n", n);
	  multi_T_chain_within_T_mcmc_step(multi_T_chain);
	
	  if(OUTPUT_TVD_VS_N){
	    if(n%10000 == 0){
	      Ndim_histogram* hist_copy = construct_copy_histogram_ndim(mcmc_out_hist);
	      normalize_ndim_histogram(hist_copy);
	      tvd = total_variation_distance(targp, hist_copy);
	      printf("#%4i %4i %12.6g\n", k, n, tvd);
	      free_histogram_ndim(hist_copy);
	    }
	  }
	}

	Ndim_histogram* hist_copy = construct_copy_histogram_ndim(mcmc_out_hist);
	normalize_ndim_histogram(hist_copy);
	tvd = total_variation_distance(targp, hist_copy);
	    printf("#%4i %4i %12.6g\n", k, mcmc_steps, tvd);
	free_histogram_ndim(hist_copy);
		printf("\n\n");
	tvd_sum += tvd;
	tvdsq_sum += tvd*tvd;
	Naccept_sum += n_accept;
	Nreject_sum += n_reject;
	free_state(the_state);
	mcmc_out_hist = free_histogram_ndim(mcmc_out_hist);
	//	printf("bottom of Nrep loop. j: %i \n", j);
      } // end loop over reps
      //sum_dsq /= Nreps*mcmc_steps;
      double tvd_avg = tvd_sum/Nreps;
      double tvd_variance = tvdsq_sum/Nreps - tvd_avg*tvd_avg;
      double tvd_stderr = sqrt(tvd_variance)/sqrt((double)Nreps);
      //      printf("%4i  %8i  %10.6f +- %10.6f  %10.6f  %10.6f \n", k, n_jumps, tvd_avg, tvd_stderr, (double)Naccept_sum/(Naccept_sum+Nreject_sum), sqrt(sum_dsq/(Nreps*mcmc_steps)));
          printf("\n\n");
      fflush(stdout);
    
    } // end loop over different proposal widths
  free_histogram_ndim(targp);
} // end of main

// ********** function definitions ******************

/* int mcmc_step(State* the_state, Proposal* prop, int Ngrid_max, Ndim_histogram* hist_ndim){ */
/*   int Ndim = the_state->n_dimensions; */
/*   int* index_array = the_state->point; */
/*   int* prop_index_array = propose(Ndim, index_array, prop, IS_BALL); */
/*   double p_of_proposed_state = F(Ndim, prop_index_array, Ngrid_max); */
/*   double p_of_current_state = the_state->prob; */
/*   int n_jumps = 0; */
/*   int dsq = 0; */
/*   // prop ratio is 1 for symmetric proposal */
/*   if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept */
  
/*     int h0max = Ngrid_max/2;   // count if jumped from < 0.5 to > 0.5 in x and y */
/*     for(int i=0; i<Ndim; i++){ */
/*       if( ( the_state->point[i] <= h0max  && prop_index_array[i] > h0max )  */
/* 	  || ( the_state->point[i] > h0max  &&  prop_index_array[i] <= h0max ) ){  */
/* 	n_jumps++; */
/*       } */
/* 	int d = the_state->point[i] - prop_index_array[i]; */
/* 	dsq += d*d; */
/*     } */
/*     //   printf("accept, dsq: %i \n", dsq); */
/*     // set state to the proposed state: */
/*     free(the_state->point);  */
/*     the_state->point = prop_index_array; */
/*     the_state->prob = p_of_proposed_state; */
/*     n_accept++; */
/*   }else{  */
/*     // reject proposed move. */
/*     //   printf("reject, dsq: %i \n", dsq); */
/*     n_reject++; */
/*     free(prop_index_array); */
/*     // dsq += 0  */
/*   } */
/*   sum_dsq += (double)dsq; */
/*   if(hist_ndim != NULL){ */
/*     double* bp = get_pointer_to_element(hist_ndim->weights, the_state->point);  */
/*     (*bp) += 1.0; */
/*     hist_ndim->total_weight += 1; */
/*   } */
/*   return n_jumps; */
/* } */

int* propose(int Ndim, int* index_array, Proposal* prop, int is_ball){
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
  /* print_array_of_int(Ndim, index_array); */
  /* printf("   ;   "); */
  /* print_array_of_int(Ndim, prop_index_array); printf("\n"); */
  return prop_index_array;
}

double F(int Ndim, int* index_array, int Ngrid_max){
  double result = 1.0;
  for(int i=0; i<Ndim; i++){
    int index = index_array[i];
    //  printf("i, index, Ngrid_max: %i %i %i \n", i, index, Ngrid_max);
    if(index < 0  || index > Ngrid_max){ return 0.0; }
    double x = (double)index/(double)Ngrid_max;
    result *= f_1dim(x);
  }
  return result;
}    

double f_1dim(double x){
  if(x<0.0 || x > 1.0){ // out of bounds, return 0
    return 0.0;
  }
  if(n_modes == 1){
    return exp(-0.5*((x-0.5)/sigma)*((x-0.5)/sigma));
  }else{
    double p1 = 0.2;
    double p2 = 0.8;
    double sig1 = sigma;
    double sig2 = sigma;
    double f = exp(-0.5* ((x-p1)/sig1)*((x-p1)/sig1)) 
      + exp(-0.5* ((x-p2)/sig2)*((x-p2)/sig2)) ;
    return f;
  }
}

double drand(void){
  return (double)rand() / ((double)RAND_MAX + 1);
}




Ndim_histogram* init_target_distribution(int Ndim, int Ngrid_max, int normalize){
  Ndim_histogram* targ_pdf = (Ndim_histogram*)malloc(sizeof(Ndim_histogram));
  targ_pdf->Ndim = Ndim;
  targ_pdf->Ngrid_max = Ngrid_max;
  targ_pdf->weights = construct_ndim_array_of_double(Ndim, Ngrid_max+1, 0.0);
  targ_pdf->total_weight = set_ndim_array_of_double_with_function(targ_pdf->weights, 1.0, f_1dim);
  if(normalize){
    normalize_ndim_histogram(targ_pdf);
  }
  return targ_pdf;
}




double total_variation_distance(Ndim_histogram* targp, Ndim_histogram* mcmc_out){
  // both distributions must be normalized ahead of time
  Ndim_array_of_double* a1 = targp->weights;
  Ndim_array_of_double* a2 = mcmc_out->weights;
  double tvd = 0.5 * sum_abs_difference_ndim_arrays_of_double(targp->weights, mcmc_out->weights);
  return tvd;
}


// *********** not used ******************

void print_array_of_int(int Nsize, int* array){
  for(int i=0; i<Nsize; i++){
    printf("%6i ", array[i]);
  }// printf("\n");
}

Ndim_array_of_int* construct_ndim_array_of_int(int Ndim, int Nsize, int init_value){
  Ndim_array_of_int* array_struct = (Ndim_array_of_int*)malloc(sizeof(Ndim_array_of_int));
  array_struct->Ndim = Ndim;
  array_struct->Nsize = Nsize;
  if(Ndim == 1){
    int* the_array = (int*)malloc(Nsize*sizeof(int));
    for(int i=0; i<Nsize; i++){
      the_array[i] = init_value;
    }
    array_struct->array = the_array; 
  }else{
    Ndim_array_of_int** the_array = (Ndim_array_of_int**)malloc(Nsize*sizeof(Ndim_array_of_int*));
    for(int i=0; i<Nsize; i++){
      the_array[i] = construct_ndim_array_of_int(Ndim-1, Nsize, init_value);
    }
    array_struct->array = the_array;
  }
  return array_struct;
}
