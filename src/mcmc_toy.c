#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "mcmc_toy_structs.h"
#include "mcmc_toy.h"

int n_accept = 0;
int n_reject = 0;
double sum_dsq = 0;

#define SIGMA (0.06)
#define IS_BALL (1)
#define N_MODES (2) // per dimension
// main:
int main(){
  int n_dimensions = 3; 
  int Ngrid_max = 95;
  if(Ngrid_max % 2 == 0){ Ngrid_max++; } // if even, add one, to make it be odd.
  int mcmc_steps  = 400000;
  int Nreps = 10;
  int burn_in_steps = 20000;
  double sigma_in_steps = SIGMA*Ngrid_max;
  printf("# sigma*Ngrid_max: %12.8f \n", SIGMA*Ngrid_max);
  printf("# mcmc_steps %i Nreps: %i \n", mcmc_steps, Nreps);
 

  Ndim_histogram* targp = init_target_distribution(n_dimensions, Ngrid_max, 1);

  int k = 1;
  for(int k=1; k<2*Ngrid_max; k = (int)(1.1*(k+1)) ) // loop over different proposals
    {    
      double tvd_sum = 0.0;
      int n_jumps = 0;
      int Naccept_sum = 0;
      int Nreject_sum = 0;
      sum_dsq = 0.0;
      for(int j=0; j<Nreps; j++){
    
	Proposal proposal = {k, 4, 0.3};
    
	// ********** do burn-in ***********
	State* the_state = construct_state(n_dimensions, Ngrid_max);
	for(int n=0; n<=burn_in_steps; n++){
	  int x = mcmc_step_ndim(the_state, &proposal, Ngrid_max, NULL);
	} 

	// ********** do post-burn-in *********
	n_accept = 0;
	n_reject = 0;
	Ndim_histogram* hist_ndim = construct_histogram_ndim(n_dimensions, Ngrid_max);
     
	for(int n=1; n<=mcmc_steps; n++){   
	  // take a mcmc step
	  n_jumps += mcmc_step_ndim(the_state, &proposal, Ngrid_max, hist_ndim);
	}
	normalize_ndim_histogram(hist_ndim);
	tvd_sum += total_variation_distance(targp, hist_ndim);
	Naccept_sum += n_accept;
	Nreject_sum += n_reject;
	free_state(the_state);
	free_histogram_ndim(hist_ndim);
      } // end loop over reps
      sum_dsq /= Nreps*mcmc_steps;
      printf("%4i  %8i  %10.6f   %10.6f  %10.6f \n", k, n_jumps, tvd_sum/Nreps, (double)Naccept_sum/(Naccept_sum+Nreject_sum), sum_dsq);
    } // end loop over different proposal widths
  free_histogram_ndim(targp);
} // end of main

// ********** function definitions ******************

int mcmc_step_ndim(State* the_state, Proposal* prop, int Ngrid_max, Ndim_histogram* hist_ndim){
  int Ndim = the_state->n_dimensions;
  int* index_array = the_state->point;
  int* prop_index_array = propose(Ndim, index_array, prop, IS_BALL);
  double p_of_proposed_state = F(Ndim, prop_index_array, Ngrid_max);
  double p_of_current_state = the_state->prob;
  int n_jumps = 0;
  int dsq = 0;
  // prop ratio is 1 for symmetric proposal
  if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept
  
    int h0max = Ngrid_max/2;   // count if jumped from < 0.5 to > 0.5 in x and y
    for(int i=0; i<Ndim; i++){
      if( ( the_state->point[i] <= h0max  && prop_index_array[i] > h0max ) 
	  || ( the_state->point[i] > h0max  &&  prop_index_array[i] <= h0max ) ){ 
	n_jumps++;
	int d = the_state->point[i] - prop_index_array[i];
	dsq += d*d;
      }
    }
    // set state to the proposed state:
    free(the_state->point); 
    the_state->point = prop_index_array;
    the_state->prob = p_of_proposed_state;
    n_accept++;
  }else{ 
    // reject proposed move.
    n_reject++;
    free(prop_index_array);
  }
  sum_dsq += (double)dsq;
  if(hist_ndim != NULL){
    double* bp = get_pointer_to_element(hist_ndim->weights, the_state->point); 
    (*bp) += 1.0;
    hist_ndim->total_weight += 1;
  }
  return n_jumps;
}

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
    if(!is_ball  || sum_sq <= Width*Width) break;
  }
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
  if(N_MODES == 1){
    return exp(-0.5*((x-0.5)/SIGMA)*((x-0.5)/SIGMA));
  }else{
    double p1 = 0.2;
    double p2 = 0.8;
    double sig1 = SIGMA;
    double sig2 = SIGMA;
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

void print_int_array(int Nsize, int* array){
  for(int i=0; i<Nsize; i++){
    printf("%6i ", array[i]);
  }printf("\n");
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
