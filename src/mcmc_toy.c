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
  int n_dimensions = 2; // other dimensionalities not implemented.
  int gridsize = 25; // gridsize+1 points, [0,gridsize]
  //  int K = 23; // max number of gridpoints to move in 1 step
  int mcmc_steps  = 4000000;
  int burn_in_steps = 50000;
  Target_distribution* targp = initialize_pi(n_dimensions, gridsize);
 
  // int k = 25;
  for(int k=1; k<=2*gridsize; k += 2)
    {
      Histogram_2d* hist2d = initialize_histogram_2d(gridsize);
      Proposal proposal = {k, 3, 1.7};
      /*
	for(int i=0; i<=gridsize; i++){
	for(int j=0; j<=gridsize; j++){
	printf("%12.5g ", targp[i][j]);
	}printf("\n");
	}
	// */
      n_accept = 0;
      n_reject = 0;
      State* the_state = initialize_state(n_dimensions, gridsize, targp);
      for(int n=0; n<=burn_in_steps; n++){
	int x = mcmc_step(the_state, &proposal, targp, hist2d);
      } 
      n_accept = 0;
      n_reject = 0;
      hist2d = initialize_histogram_2d(gridsize);
      int xy_jumps = 0;
      for(int n=0; n<=mcmc_steps; n++){   
	// take a mcmc step
	xy_jumps += mcmc_step(the_state, &proposal, targp, hist2d);
	//  printf("%8i %8i %8i %10.5g \n", n, the_state->point[0], the_state->point[1], the_state->prob);
      }
      double tvd = total_variation_distance(targp, hist2d);
      printf("%8i  %8i  %8i  %10.6f   %8i  %8i  %8i \n", k, xy_jumps, mcmc_steps, tvd, n_accept, n_reject, n_accept+n_reject);
    }
}
// end of main

// initialize a state to
// a point chosen u.a.r. among the (Ngrid+1)^n_dimensions points
State* initialize_state(int n_dimensions, int Ngrid, Target_distribution* targp){
  State* the_state = (State*)malloc(1*sizeof(State));
  assert(n_dimensions == 2); // other dimesions not implemented at present.
  the_state->n_dimensions = n_dimensions;
  the_state->point = (int*)malloc(n_dimensions*sizeof(int));
  for(int i=0; i<n_dimensions; i++){
    the_state->point[i] = (int)( drand() * (Ngrid+1) );
  }
  the_state->prob = targp->pi[the_state->point[0]][the_state->point[1]];
  the_state->xhalf = (the_state->point[0] <= Ngrid/2)? 0 : 1;
  the_state->yhalf = (the_state->point[1] <= Ngrid/2)? 0 : 1;
  return the_state;
}

int mcmc_step(State* the_state, Proposal* prop, Target_distribution* targp, Histogram_2d* hist2d){
  int prop_width = (prop->p1 >= 1 || drand() < prop->p1)? prop->W1 : prop->W2; 
  int Ngrid = targp->Ngrid;
  int i = the_state->point[0];
  int j = the_state->point[1];
  int iprop, jprop;
  while(1){ // don't propose same point
    iprop = propose_1dim(i, prop_width, targp->Ngrid);
    jprop = propose_1dim(j, prop_width, targp->Ngrid);
    if(iprop != i || jprop != j) break;
  }
  int x_jumps = 0; // 0 or 1
  int y_jumps = 0; // 0 or 1
  double p_of_proposed_state = f( (double)iprop/(double)targp->Ngrid ) * f( (double)jprop/(double)targp->Ngrid );
  // targp->pi[iprop][jprop];
  double p_of_current_state = the_state->prob; // targp->pi[i][j];
  // prop ratio is 1 for symmetric proposal
  if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept
    // count if jumped from < 0.5 to > 0.5 in x and y
    if( (the_state->point[0] <= Ngrid/2  && iprop > Ngrid/2) || ( the_state->point[0] >= Ngrid/2  && iprop < Ngrid/2 ) ){ 
      x_jumps++;
    }
    if( (the_state->point[1] <= Ngrid/2  && jprop > Ngrid/2) || ( the_state->point[1] >= Ngrid/2  && jprop < Ngrid/2 ) ){ 
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
  hist2d->bin_counts[the_state->point[0]][the_state->point[1]]++;
  hist2d->total_count++;
  //    printf("%8i %8i %10.5g \n", the_state->point[0], the_state->point[1], the_state->prob);
  return x_jumps + y_jumps;
}

double total_variation_distance(Target_distribution* targp, Histogram_2d* hist2d){
  int Ngrid = targp->Ngrid;
  double tvd = 0;
  for(int i=0; i<= Ngrid; i++){
    for(int j=0; j<= Ngrid; j++){
      double mcmc_prob = (double)hist2d->bin_counts[i][j]/(double)hist2d->total_count;
      //          printf("%8.4f %8.4f %8.4f \n", mcmc_prob, targp->pi[i][j], tvd);
      tvd += fabs(mcmc_prob - targp->pi[i][j]);
    }
  }
  tvd *= 0.5;
  return tvd;
}

int propose_1dim(int i, int Width, int Ngrid){ // 1 dimensional proposal
  int delta_i;
  delta_i = (int)((2*Width+1)*drand()) - Width;
  int iprop = i + delta_i; 
  //  printf("width, deltai, iprop: %8i %8i %8i \n", Width, delta_i, iprop);
  if(REFLECT){
    while(iprop < 0 || iprop > Ngrid){
      if(iprop < 0){
	iprop = -iprop;
      }
      if(iprop > Ngrid){
	iprop = Ngrid - (iprop - Ngrid); // so, e.g. Ngrid+1 => Ngrid-1, etc.
      }
    }
  }else{
    // if(iprop < 0 || iprop > Ngrid){ iprop = i; } // if proposal is out of bounds, just propose i (unchanged)
  }
  return iprop;
}

Target_distribution* initialize_pi(int n_dimensions, int Ngrid){
  assert(n_dimensions == 2);
  Target_distribution* target = (Target_distribution*)malloc(1*sizeof(Target_distribution));
  target->n_dimensions = n_dimensions;
  target->Ngrid = Ngrid;
  double** pi = (double**)malloc( (Ngrid+1) * sizeof(double*) );
  for(int i=0; i<=Ngrid; i++){
    pi[i] = (double*)malloc( (Ngrid+1) * sizeof(double) );
  }

  double sum_pi = 0.0;
  for(int i=0; i<=Ngrid; i++){
    double fx = f( (double)i/(double)Ngrid );
    for(int j=0; j<=Ngrid; j++){
      double fy = f( (double)j/(double)Ngrid );
      double point_prob = fx*fy;
      pi[i][j] = point_prob;
      sum_pi += point_prob;
    }
  }
  // normalize pi 
  for(int i=0; i<=Ngrid; i++){
    for(int j=0; j<=Ngrid; j++){
      pi[i][j] /= sum_pi;
    }
  }
  target->pi = pi;
  return target;
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

Histogram_2d* initialize_histogram_2d(int Ngrid){
  Histogram_2d* histogram = (Histogram_2d*)malloc(1*sizeof(Histogram_2d));
  histogram->Ngrid = Ngrid;
  histogram->total_count = 0;
  int** bin_counts = (int**)malloc( (Ngrid+1) * sizeof(double*) );
  for(int i=0; i<=Ngrid; i++){
    bin_counts[i] = (int*)malloc( (Ngrid+1) * sizeof(double) );
  }

  for(int i=0; i<=Ngrid; i++){
    for(int j=0; j<=Ngrid; j++){
      bin_counts[i][j] = 0;
    }
  }
  histogram->bin_counts = bin_counts;
  return histogram;
}

Ndim_array_of_double* initialize_ndim_array_of_double(int Ndim, int Nsize, double init_value){
  Ndim_array_of_double* array_struct;
  array_struct->Ndim = Ndim;
  array_struct->Nsize = Nsize;
  if(Ndim == 1){
    double* the_array = (double*)malloc(Nsize*sizeof(double));
    for(int i=0; i<Nsize; i++){
      the_array[i] = init_value;
    }
     
  }else{
    Ndim_array_of_double** the_array = (Ndim_array_of_double**)malloc(Nsize*sizeof(Ndim_array_of_double*));
    for(int i=0; i<Nsize; i++){
      the_array[i] = initialize_ndim_array_of_double(Ndim-1, Nsize, init_value);
    }
  }
  array_struct
  return array_struct;
}
