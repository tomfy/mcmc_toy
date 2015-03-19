#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "mcmc_toy.h"

#define REFLECT 0

int n_accept = 0;
int n_reject = 0;

// main:
int main(){
  int n_dimensions = 5; 
  int Ngrid_max = 15;
  if(Ngrid_max % 2 == 0){ Ngrid_max++; } // if even add one, to make it be odd, in case it is not already
  //  int gridsize = Ngrid_max + 1; // gridsize points, [0,Ngrid_max], gridsize will be even
  //  int K = 23; // max number of gridpoints to move in 1 step
  int mcmc_steps  = 100000;
  int burn_in_steps = 20000;
  //  Target_distribution* targp = initialize_pi(n_dimensions, Ngrid_max);
  // printf("after initialize_pi\n");
  // int k = 25;
  for(int k=1; k<=2*(Ngrid_max+1); k = (int)(1.1*(k+1)) )
    {
      //     Histogram_2d* hist2d = initialize_histogram_2d(Ngrid_max);
      //    printf("before initialize_histogram_ndim\n");
      Histogram_ndim* hist_ndim = initialize_histogram_ndim(n_dimensions, Ngrid_max);
      //        printf("after initialize_histogram_ndim\n");
      Proposal proposal = {k, 3, 1.7};
    
      n_accept = 0;
      n_reject = 0;
      State* the_state = initialize_state(n_dimensions, Ngrid_max);
      for(int n=0; n<=burn_in_steps; n++){
	//		printf("before and ... \n");
	int x = mcmc_step_ndim(the_state, &proposal, Ngrid_max, hist_ndim);
	//	printf("after mcmc_step_ndim n= %i\n", n);
      } 
      n_accept = 0;
      n_reject = 0;
      //  hist2d = initialize_histogram_2d(Ngrid_max);
      //    printf("before init hist ndim \n");
      hist_ndim = initialize_histogram_ndim(n_dimensions, Ngrid_max);
      //     printf("after init hist ndim \n");
      int xy_jumps = 0;
      for(int n=1; n<=mcmc_steps; n++){   
	// take a mcmc step
	xy_jumps += mcmc_step_ndim(the_state, &proposal, Ngrid_max, hist_ndim);
	//  printf("%8i %8i %8i %10.5g \n", n, the_state->point[0], the_state->point[1], the_state->prob);
      }
      double tvd = -1; //total_variation_distance(targp, hist_ndim);
      printf("%8i  %8i  %8i  %10.6f   %8i  %8i  %8i \n", k, xy_jumps, mcmc_steps, tvd, n_accept, n_reject, n_accept+n_reject);
    }
}
// end of main

// ********** function definitions ******************

// initialize a state to a point chosen u.a.r. among the (Ngrid_max+1)^n_dimensionsoints
State* initialize_state(int n_dimensions, int Ngrid_max){ // , Target_distribution* targp){
  State* the_state = (State*)malloc(1*sizeof(State));
  //assert(n_dimensions == 2); // other dimensions not implemented at present.
  the_state->n_dimensions = n_dimensions;
  the_state->point = (int*)malloc(n_dimensions*sizeof(int));
  for(int i=0; i<n_dimensions; i++){
    the_state->point[i] = (int)( drand() * (Ngrid_max+1) );
  }
  the_state->prob = F(n_dimensions, the_state->point, Ngrid_max); // targp->pi[the_state->point[0]][the_state->point[1]];
  //  the_state->xhalf = (the_state->point[0] <= Ngrid_max/2)? 0 : 1;
  //  the_state->yhalf = (the_state->point[1] <= Ngrid_max/2)? 0 : 1;
  return the_state;
}


int mcmc_step_ndim(State* the_state, Proposal* prop, int Ngrid_max, Histogram_ndim* hist_ndim){
  int Ndim = hist_ndim->Ndim;
  //  int prop_width = (prop->p1 >= 1 || drand() < prop->p1)? prop->W1 : prop->W2; 
  //  int Ngrid_max = targp->Ngrid_max;
  int* index_array = the_state->point;
  // printf("point: ");
  int* prop_index_array = propose(Ndim, index_array, prop);
  //  printf("after propose \n");
  // print_int_array(Ndim, index_array);
  // print_int_array(Ndim, prop_index_array);
  int n_jumps = 0; // 0 or 1
  double p_of_proposed_state = F(Ndim, prop_index_array, Ngrid_max);
  //  printf("prob of proposed state: %12.6g \n", p_of_proposed_state);
  // f( (double)iprop/(double)Ngrid_max ) * f( (double)jprop/(double)Ngrid_max );
  double p_of_current_state = the_state->prob;
  // prop ratio is 1 for symmetric proposal
  if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept
    // count if jumped from < 0.5 to > 0.5 in x and y
    int h0max = Ngrid_max/2;
    for(int i=0; i<Ndim; i++){
      //    printf("i: %i \n", i);
      if( ( the_state->point[i] <= h0max  && prop_index_array[i] > h0max ) || ( the_state->point[i] > h0max  &&  prop_index_array[i] <= h0max ) ){ 
	n_jumps++;
      }
    }
    // set state to the proposed state:
    the_state->point = prop_index_array;
    the_state->prob = p_of_proposed_state;
    n_accept++;
    // printf("%i %i   %i %i   %i %i \n", the_state->point[0], the_state->point[1], iprop, jprop, x_jumps, y_jumps);
  }else{ 
    // reject proposed move.
    n_reject++;
    free(prop_index_array);
  }
  // ((int*)hist_ndim->bin_count[the_state->point[0]].array)[the_state->point[1]]++;
  //  printf("after step, before get_bin_pointer.\n");
  // print_int_array(Ndim, index_array);
  //  print_int_array(Ndim, prop_index_array);
  //printf("before get_bin_pointer \n");
  int* bp = get_bin_pointer(hist_ndim->bin_count, the_state->point); 
  //  printf("after get_bin_pointer \n");
  (*bp)++;
  hist_ndim->total_count++;
  //    printf("%8i %8i %10.5g \n", the_state->point[0], the_state->point[1], the_state->prob);
  return n_jumps;
}

int* propose(int Ndim, int* index_array, Proposal* prop){
  int* prop_index_array = (int*)malloc(Ndim*sizeof(int)); 
  int Width = (drand() < prop->p1)? prop->W1 : prop->W2;
  for(int i = 0; i<Ndim; i++){
    int delta_i = (int)((2*Width+1)*drand()) - Width;
    prop_index_array[i] = index_array[i] + delta_i;
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
    result *= f(x);
  }
  return result;
}    

double f(double x){
  if(x<0.0 || x > 1.0){ // out of bounds, return 0
    return 0.0;
  }
  double p1 = 0.2;
  double p2 = 0.8;
  double sig1 = 0.06;
  double sig2 = 0.06;
  double f = exp(-0.5* ((x-p1)/sig1)*((x-p1)/sig1)) 
    + exp(-0.5* ((x-p2)/sig2)*((x-p2)/sig2))
    ;
  return f;
}

double drand(void){
  return (double)rand() / ((double)RAND_MAX + 1);
}



Histogram_ndim* initialize_histogram_ndim(int Ndim, int Ngrid_max){
  Histogram_ndim* histogram = (Histogram_ndim*)malloc(1*sizeof(Histogram_ndim));
  histogram->Ndim = Ndim;
  histogram->Ngrid_max = Ngrid_max;
  histogram->total_count = 0;
  //  printf("before initialize_ndim_array_of_int ...\n");
  histogram->bin_count = initialize_ndim_array_of_int(Ndim, Ngrid_max+1, 0);
  // printf("after initialize_ndim_array_of_int ...\n");
  return histogram;
}

Ndim_array_of_double* initialize_ndim_array_of_double(int Ndim, int Nsize, double init_value){
  Ndim_array_of_double* array_struct = (Ndim_array_of_double*)malloc(sizeof(Ndim_array_of_double));
  array_struct->Ndim = Ndim;
  array_struct->Nsize = Nsize;
  if(Ndim == 1){
    double* the_array = (double*)malloc(Nsize*sizeof(double));
    for(int i=0; i<Nsize; i++){
      the_array[i] = init_value;
    }
    array_struct->array = the_array; 
  }else{
    Ndim_array_of_double** the_array = (Ndim_array_of_double**)malloc(Nsize*sizeof(Ndim_array_of_double*));
    for(int i=0; i<Nsize; i++){
      the_array[i] = initialize_ndim_array_of_double(Ndim-1, Nsize, init_value);
    }
    array_struct->array = the_array;
  }
  return array_struct;
}

// ((int*)hist_ndim->bin_count[the_state->point[0]].array)[the_state->point[1]]++;
int* get_bin_pointer(Ndim_array_of_int* A, int* index_array){
  int Ndim = A->Ndim;
  Ndim_array_of_int* B = A;
  for(int j=0; j<Ndim-1; j++){
    //    printf("j, index_array[j]: %i %i \n", j, index_array[j]);
    Ndim_array_of_int* C = ((Ndim_array_of_int**)B->array)[index_array[j]];
    B = C;
  }
  // at this point, B should point to 1 dim 'Ndim_array_of_int', whose array field is a regular array of double:
  int* array_of_double = (int*)B->array;
  return array_of_double + index_array[Ndim-1];
}

Ndim_array_of_int* initialize_ndim_array_of_int(int Ndim, int Nsize, int init_value){
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
      the_array[i] = initialize_ndim_array_of_int(Ndim-1, Nsize, init_value);
    }
    array_struct->array = the_array;
  }
  return array_struct;
}

void print_int_array(int Nsize, int* array){
  for(int i=0; i<Nsize; i++){
    printf("%6i ", array[i]);
  }printf("\n");
}


double total_variation_distance(Target_distribution* targp, Histogram_ndim* hist_ndim){
  int Ngrid_max = targp->Ngrid_max;
  double tvd = 0;
  int point[2] = {0, 0};
  for(int i=0; i<= Ngrid_max; i++){
  point[0] = i;
  for(int j=0; j<= Ngrid_max; j++){
  point[1] = j;
  int* bp = get_bin_pointer(hist_ndim->bin_count, point);
  int bin_count = *bp;
  double mcmc_prob = (double)bin_count/(double)hist_ndim->total_count;
  //         printf("%8.4f %8.4f %8.4f \n", mcmc_prob, targp->pi[i][j], tvd);
  tvd += fabs(mcmc_prob - targp->pi[i][j]);
  }
  }
  tvd *= 0.5;
  return tvd;
  }




// *********** obsolete ? stuff ********************

int propose_1dim(int i, int Width, int Ngrid_max){ // 1 dimensional proposal
  int delta_i;
  delta_i = (int)((2*Width+1)*drand()) - Width;
  int iprop = i + delta_i; 
  //  printf("width, deltai, iprop: %8i %8i %8i \n", Width, delta_i, iprop);
  if(REFLECT){
    while(iprop < 0 || iprop > Ngrid_max){
      if(iprop < 0){
	iprop = -iprop;
      }
      if(iprop > Ngrid_max){
	iprop = Ngrid_max - (iprop - Ngrid_max); // so, e.g. Ngrid_max+1 => Ngrid_max-1, etc.
      }
    }
  }else{
    // just go ahead and propose something out of bounds - it will have pi = 0, so won't go there.
  }
  return iprop;
}

// int mcmc_step(State* the_state, Proposal* prop, Target_distribution* targp, Histogram_2d* hist2d){
int mcmc_step_2d(State* the_state, Proposal* prop, int Ngrid_max, Histogram_ndim* hist_ndim){
  int prop_width = (prop->p1 >= 1 || drand() < prop->p1)? prop->W1 : prop->W2; 
  //  int Ngrid_max = targp->Ngrid_max;
  int i = the_state->point[0];
  int j = the_state->point[1];
  int iprop, jprop;
  while(1){ // don't propose same point
    iprop = propose_1dim(i, prop_width, Ngrid_max);
    jprop = propose_1dim(j, prop_width, Ngrid_max);
    if(iprop != i || jprop != j) break;
  }
  int x_jumps = 0; // 0 or 1
  int y_jumps = 0; // 0 or 1
  double p_of_proposed_state = f( (double)iprop/(double)Ngrid_max ) * f( (double)jprop/(double)Ngrid_max );
  double p_of_current_state = the_state->prob;
  // prop ratio is 1 for symmetric proposal
  if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept
    // count if jumped from < 0.5 to > 0.5 in x and y
    int h0max = Ngrid_max/2;
    if( ( the_state->point[0] <= h0max  && iprop > h0max ) || ( the_state->point[0] > h0max  &&  iprop <= h0max ) ){ 
      x_jumps++;
    }
    if( ( the_state->point[1] <= h0max  && jprop > h0max ) || ( the_state->point[1] > h0max  &&  jprop <= h0max ) ){ 
      y_jumps++;
    }
    // set state to the proposed state:
    the_state->point[0] = iprop; 
    the_state->point[1] = jprop;
    the_state->prob = p_of_proposed_state;
    n_accept++;
    // printf("%i %i   %i %i   %i %i \n", the_state->point[0], the_state->point[1], iprop, jprop, x_jumps, y_jumps);
  }else{ 
    // reject proposed move; do nothing
    n_reject++;
  }
  // ((int*)hist_ndim->bin_count[the_state->point[0]].array)[the_state->point[1]]++;
  int* bp = get_bin_pointer(hist_ndim->bin_count, the_state->point); 
  (*bp)++;
  hist_ndim->total_count++;
  //    printf("%8i %8i %10.5g \n", the_state->point[0], the_state->point[1], the_state->prob);
  return x_jumps + y_jumps;
}


/*
  Target_distribution* initialize_pi(int n_dimensions, int Ngrid_max){
  // assert(n_dimensions == 2);
  Target_distribution* target = (Target_distribution*)malloc(1*sizeof(Target_distribution));
  target->n_dimensions = n_dimensions;
  target->Ngrid_max = Ngrid_max;
  Ndim_array_of_double* pi = (Ndim_array_of_double*)malloc( sizeof(Ndim_array_of_double) );
  initialize_ndim_array_of_double(n_dimensions, Ngrid_max+1, 0);

  double sum_pi = 0.0;
  for(int i=0; i<=Ngrid_max; i++){
  double fx = f( (double)i/(double)Ngrid_max );
  for(int j=0; j<=Ngrid_max; j++){
  double fy = f( (double)j/(double)Ngrid_max );
  double point_prob = fx*fy;
  pi[i][j] = point_prob;
  sum_pi += point_prob;
  }
  }
  // normalize pi 
  for(int i=0; i<=Ngrid_max; i++){
  for(int j=0; j<=Ngrid_max; j++){
  pi[i][j] /= sum_pi;
  }
  }
  target->pi = pi;
  return target;
  }
*/


Histogram_2d* initialize_histogram_2d(int Ngrid_max){
  Histogram_2d* histogram = (Histogram_2d*)malloc(1*sizeof(Histogram_2d));
  histogram->Ngrid_max = Ngrid_max;
  histogram->total_count = 0;
  int** bin_count = (int**)malloc( (Ngrid_max+1) * sizeof(double*) );
  for(int i=0; i<=Ngrid_max; i++){
    bin_count[i] = (int*)malloc( (Ngrid_max+1) * sizeof(double) );
  }

  for(int i=0; i<=Ngrid_max; i++){
    for(int j=0; j<=Ngrid_max; j++){
      bin_count[i][j] = 0;
    }
  }
  histogram->bin_count = bin_count;
  return histogram;
}
