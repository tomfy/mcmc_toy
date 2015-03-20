#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "mcmc_toy.h"

int n_accept = 0;
int n_reject = 0;

#define SIGMA (0.06)

// main:
int main(){
  int n_dimensions = 2; 
  int Ngrid_max = 35;
  if(Ngrid_max % 2 == 0){ Ngrid_max++; } // if even, add one, to make it be odd, in case it is not already
  int mcmc_steps  = 40000;
  int burn_in_steps = 20000;
  double sigma_in_steps = SIGMA*Ngrid_max;
  printf("%12.8f \n", SIGMA*Ngrid_max);
  int Nreps = 5;
  int Naccept_sum = 0;
  int Nreject_sum = 0;

  Histogram_ndim* targp = init_target_distribution(n_dimensions, Ngrid_max, 1);

  int k = 1;
  for(int k=1; k<=sigma_in_steps+2; k = (int)(1.1*(k+1)) ) // loop over different proposals
    {    
      double tvd_sum = 0.0;
      int n_jumps = 0;
      for(int j=0; j<Nreps; j++){
    
	Proposal proposal = {k, 2, 1.125};
    
	// ********** do burn-in ***********
	State* the_state = initialize_state(n_dimensions, Ngrid_max);
	for(int n=0; n<=burn_in_steps; n++){
	  int x = mcmc_step_ndim(the_state, &proposal, Ngrid_max, NULL);
	} 

	// ********** do post-burn-in *********
	n_accept = 0;
	n_reject = 0;
	Histogram_ndim* hist_ndim = initialize_histogram_ndim(n_dimensions, Ngrid_max);
     
	for(int n=1; n<=mcmc_steps; n++){   
	  // take a mcmc step
	  n_jumps += mcmc_step_ndim(the_state, &proposal, Ngrid_max, hist_ndim);
	  //  printf("%8i %8i %8i %10.5g \n", n, the_state->point[0], the_state->point[1], the_state->prob);
	}
	//    printf("before normalize mcmc out hist.\n");
	normalize_pdf_ndim(hist_ndim);
	//     printf("after normalize mcmc out hist. total_weight: %g \n", hist_ndim->total_weight);
	tvd_sum += total_variation_distance(targp, hist_ndim);
	Naccept_sum += n_accept;
	Nreject_sum += n_reject;
	free_state(the_state);
	free_histogram_ndim(hist_ndim);
      } // end loop over reps
      printf("%8i  %8i  %10.6f   %8i  %8i\n", k, n_jumps, tvd_sum/Nreps, Naccept_sum, Nreject_sum);
    } // end loop over different proposal widths
  free_histogram_ndim(targp);
} // end of main

// ********** function definitions ******************

// initialize a state to a point chosen u.a.r. among the (Ngrid_max+1)^n_dimensionsoints
State* initialize_state(int n_dimensions, int Ngrid_max){ // , Target_distribution* targp){
  State* the_state = (State*)malloc(sizeof(State));
  the_state->n_dimensions = n_dimensions;
  the_state->point = (int*)malloc(n_dimensions*sizeof(int));
  for(int i=0; i<n_dimensions; i++){
    the_state->point[i] = (int)( drand() * (Ngrid_max+1) );
  }
  the_state->prob = F(n_dimensions, the_state->point, Ngrid_max); // targp->pi[the_state->point[0]][the_state->point[1]];
  return the_state;
}

void free_state(State* s){
  free(s->point);
  free(s);
}

int mcmc_step_ndim(State* the_state, Proposal* prop, int Ngrid_max, Histogram_ndim* hist_ndim){
  int Ndim = the_state->n_dimensions;
  int* index_array = the_state->point;
  int* prop_index_array = propose(Ndim, index_array, prop);
  double p_of_proposed_state = F(Ndim, prop_index_array, Ngrid_max);
  double p_of_current_state = the_state->prob;
  int n_jumps = 0;
  // prop ratio is 1 for symmetric proposal
  if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept
  
    int h0max = Ngrid_max/2;   // count if jumped from < 0.5 to > 0.5 in x and y
    for(int i=0; i<Ndim; i++){
      if( ( the_state->point[i] <= h0max  && prop_index_array[i] > h0max ) 
	  || ( the_state->point[i] > h0max  &&  prop_index_array[i] <= h0max ) ){ 
	n_jumps++;
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
  if(hist_ndim != NULL){
    double* bp = get_bin_pointer(hist_ndim->weights, the_state->point); 
    (*bp) += 1.0;
    hist_ndim->total_weight += 1;
  }
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
    result *= f_1dim(x);
  }
  return result;
}    

double f_1dim(double x){
  if(x<0.0 || x > 1.0){ // out of bounds, return 0
    return 0.0;
  }
  double p1 = 0.5; //0.2;
  double p2 = 0.8;
  double sig1 = SIGMA;
  double sig2 = SIGMA;
  double f = exp(-0.5* ((x-p1)/sig1)*((x-p1)/sig1)) 
    //  + exp(-0.5* ((x-p2)/sig2)*((x-p2)/sig2))
    ;
  return f;
}

double drand(void){
  return (double)rand() / ((double)RAND_MAX + 1);
}

Histogram_ndim* initialize_histogram_ndim(int Ndim, int Ngrid_max){
  Histogram_ndim* histogram = (Histogram_ndim*)malloc(sizeof(Histogram_ndim));
  histogram->Ndim = Ndim;
  histogram->Ngrid_max = Ngrid_max;
  histogram->total_weight = 0;
  //  printf("before initialize_ndim_array_of_int ...\n");
  histogram->weights = initialize_ndim_array_of_double(Ndim, Ngrid_max+1, 0.0);
  // printf("after initialize_ndim_array_of_int ...\n");
  return histogram;
}

void free_histogram_ndim(Histogram_ndim* h){
  free_ndim_array_of_double(h->weights);
  free(h);
}

Ndim_array_of_double* initialize_ndim_array_of_double(int Ndim, int Nsize, double init_value){
  Ndim_array_of_double* array_struct = (Ndim_array_of_double*)malloc(sizeof(Ndim_array_of_double));
  array_struct->Ndim = Ndim;
  array_struct->Nsize = Nsize;
  //  printf("in initialize_ndim_array_of_double, Ndim, Nsize: %i %i \n", Ndim, Nsize);
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

void free_ndim_array_of_double(Ndim_array_of_double* A){
  if(A->Ndim == 1){
    free((double*)A->array);
  }else{
    Ndim_array_of_double** a = ((Ndim_array_of_double**)A->array);
    for(int i=0; i<A->Nsize; i++){      
      free_ndim_array_of_double( a[i] );
    }  
  }
  free(A);
}

double set_ndim_array_of_double_with_function(Ndim_array_of_double* array_struct,
					      double outer_value, double_func_ptr function){
  // for probs. of form Prod(i=1,ndim)(f(x_i)), give it a 1-dim function of x (in [0,1]) as the last argument
  int Ndim = array_struct->Ndim;
  int Nsize = array_struct->Nsize;
  //  printf("ZZZ %i %i \n", Ndim, Nsize);
  double sum = 0.0;
  if(Ndim == 1){
    for(int i=0; i<Nsize; i++){
      double x = (double)i/(double)(Nsize-1);
      //    printf("%8i %8.5g  %8.5g \n", i, outer_value, function(x));
      double innermost_value =  outer_value*function(x);
      ((double*)array_struct->array)[i] = innermost_value;
      sum += innermost_value;
      //   printf("%8.5f ", innermost_value);
    }//printf("\n");
  }else{    
    for(int i=0; i<Nsize; i++){   
      double x = (double)i/(double)(Nsize-1);
      Ndim_array_of_double* inner_array_struct = ((Ndim_array_of_double**)array_struct->array)[i];
      sum +=  set_ndim_array_of_double_with_function(inner_array_struct, outer_value*function(x), function);
    }
  }
  return sum;
}

void normalize_pdf_ndim(Histogram_ndim* pdf){
  // for probs. of form Prod(i=1,ndim)(f(x_i)), give it a 1-dim function of x (in [0,1]) as the last argument
  /* int Ndim = pdf->Ndim; */
  /* int Nsize = pdf->Nsize; */
  /* int ndim_array = pdf->weights; */
  //  printf("in normalize_pdf... total_weight: %12.6g \n", pdf->total_weight);
  multiply_ndim_array_of_double_by_scalar(pdf->weights, 1/(pdf->total_weight));
  pdf->total_weight = 1.0; // normalized now, so sum is 1.0
}

void multiply_ndim_array_of_double_by_scalar(Ndim_array_of_double* array_struct, double multiplier){
  int Ndim = array_struct->Ndim;
  int Nsize = array_struct->Nsize;
  if(Ndim == 1){
    for(int i=0; i<Nsize; i++){
      ((double*)array_struct->array)[i] *= multiplier;
    }
  }else{    
    for(int i=0; i<Nsize; i++){   
      Ndim_array_of_double* inner_array_struct = ((Ndim_array_of_double**)array_struct->array)[i];
      multiply_ndim_array_of_double_by_scalar(inner_array_struct, multiplier);
    }
  }
}


double sum_ndim_array_of_double(Ndim_array_of_double* array_struct){
  int Ndim = array_struct->Ndim;
  int Nsize = array_struct->Nsize;
  double sum = 0.0;
  // printf("top. Ndim, Nsize, sum: %i %i %g \n", Ndim, Nsize, sum);
  if(Ndim == 1){
    //double sum = 0.0;
    double* the_array = (double*)array_struct->array;
    for(int i=0; i<Nsize; i++){
      sum += the_array[i];
    }
    return sum; 
  }else{
    //double sum = 0.0;
    Ndim_array_of_double** the_array = (Ndim_array_of_double**)array_struct->array;
    for(int i=0; i<Nsize; i++){
      sum += sum_ndim_array_of_double(the_array[i]);
    }
  }
  //  printf("bottom. Ndim, Nsize, sum: %i %i %g \n", Ndim, Nsize, sum);
  return sum;
}

double sum_abs_difference_ndim_arrays_of_double(Ndim_array_of_double* a1, Ndim_array_of_double* a2){
  int Ndim = a1->Ndim;
  int Nsize = a1->Nsize;
  //  printf("%i %i %i %i \n", Ndim, Nsize, a2->Ndim, a2->Nsize);
  assert(a2->Ndim == Ndim   &&  a2->Nsize == Nsize);
  double sum_abs_difference = 0.0;
  if(Ndim == 1){
    //   printf("Ndim == 1 branch\n");
    //double sum = 0.0;
    double* A1 = (double*)a1->array;
    double* A2 = (double*)a2->array;
    for(int i=0; i<Nsize; i++){
      //   printf("i, A1[i], A2[2]: %i  %g %g \n", i, A1[i], A2[i]);
      sum_abs_difference += fabs(A1[i] - A2[i]);
    }
    //   printf("row sum: %g\n", sum_abs_difference);
    //   return sum_abs_difference; 
  }else{
    //double sum = 0.0;
    //   printf("Ndim > 1 branch\n");
    Ndim_array_of_double** A1 = (Ndim_array_of_double**)a1->array;
    Ndim_array_of_double** A2 = (Ndim_array_of_double**)a2->array;
    for(int i=0; i<Nsize; i++){
      sum_abs_difference += sum_abs_difference_ndim_arrays_of_double(A1[i], A2[i]);
      //    printf("i, plane sumabsdiff: %i  %g \n", i, sum_abs_difference);
    }
  }
  return sum_abs_difference;
}

Histogram_ndim* init_target_distribution(int Ndim, int Ngrid_max, int normalize){
  Histogram_ndim* targ_pdf = (Histogram_ndim*)malloc(sizeof(Histogram_ndim));
  targ_pdf->Ndim = Ndim;
  targ_pdf->Ngrid_max = Ngrid_max;
  targ_pdf->weights = initialize_ndim_array_of_double(Ndim, Ngrid_max+1, 0.0);
  targ_pdf->total_weight = set_ndim_array_of_double_with_function(targ_pdf->weights, 1.0, f_1dim);
  if(normalize){
    normalize_pdf_ndim(targ_pdf);
  }
  return targ_pdf;
}


// ((int*)hist_ndim->bin_count[the_state->point[0]].array)[the_state->point[1]]++;
double* get_bin_pointer(Ndim_array_of_double* A, int* index_array){
  int Ndim = A->Ndim;
  Ndim_array_of_double* B = A;
  for(int j=0; j<Ndim-1; j++){
    //    printf("j, index_array[j]: %i %i \n", j, index_array[j]);
    Ndim_array_of_double* C = ((Ndim_array_of_double**)B->array)[index_array[j]];
    B = C;
  }
  // at this point, B should point to 1 dim 'Ndim_array_of_int', whose array field is a regular array of double:
  double* array_of_double = (double*)B->array;
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

void print_ndim_array_of_double(Ndim_array_of_double* A){
  int Ndim = A->Ndim;
  int Nsize = A->Nsize;
  //  printf("In print_ndim_array_of_double, Ndim, Nsize: %i %i \n", Ndim, Nsize);
  if(Ndim == 1){
    printf("%4i %4i   ", Ndim, Nsize);
    for(int i=0; i<Nsize; i++){
      printf("%10.8f ", ((double*)A->array)[i]);
    }printf("\n");
    
  }else{
    printf("%4i %4i \n", Ndim, Nsize);
    for(int i=0; i<Nsize; i++){
      print_ndim_array_of_double( ((Ndim_array_of_double**)A->array)[i] );
    }
  }
}

double total_variation_distance(Histogram_ndim* targp, Histogram_ndim* mcmc_out){
  // both distributions must be normalized ahead of time
  Ndim_array_of_double* a1 = targp->weights;
  //printf("targp n elems: %i \n", a1->Nsize);
  Ndim_array_of_double* a2 = mcmc_out->weights;
  /* printf("targp Ndim, mcmc hist Ndim: %i %i \n", targp->weights->Ndim, mcmc_out->weights->Ndim); */
  /* printf("before sum_abs...targp n elems: %i \n", targp->weights->Nsize); */
  /* printf("before sum_abs...mcmc_out hist n elems: %i \n", mcmc_out->weights->Nsize); */
  double tvd = 0.5 * sum_abs_difference_ndim_arrays_of_double(targp->weights, mcmc_out->weights);
  return tvd;
}

// *********** obsolete ? stuff ********************

int propose_1dim(int i, int Width, int Ngrid_max){ // 1 dimensional proposal
  int delta_i;
  delta_i = (int)((2*Width+1)*drand()) - Width;
  int iprop = i + delta_i; 
  //  printf("width, deltai, iprop: %8i %8i %8i \n", Width, delta_i, iprop);
  /* if(REFLECT){ */
  /*   while(iprop < 0 || iprop > Ngrid_max){ */
  /*     if(iprop < 0){ */
  /* 	iprop = -iprop; */
  /*     } */
  /*     if(iprop > Ngrid_max){ */
  /* 	iprop = Ngrid_max - (iprop - Ngrid_max); // so, e.g. Ngrid_max+1 => Ngrid_max-1, etc. */
  /*     } */
  /*   } */
  /* }else{ */
  /*   // just go ahead and propose something out of bounds - it will have pi = 0, so won't go there. */
  /* } */
  return iprop;
}

// int mcmc_step(State* the_state, Proposal* prop, Target_distribution* targp, Histogram_2d* hist2d){
/* int mcmc_step_2d(State* the_state, Proposal* prop, int Ngrid_max, Histogram_ndim* hist_ndim){ */
/*   int prop_width = (prop->p1 >= 1 || drand() < prop->p1)? prop->W1 : prop->W2;  */
/*   //  int Ngrid_max = targp->Ngrid_max; */
/*   int i = the_state->point[0]; */
/*   int j = the_state->point[1]; */
/*   int iprop, jprop; */
/*   while(1){ // don't propose same point */
/*     iprop = propose_1dim(i, prop_width, Ngrid_max); */
/*     jprop = propose_1dim(j, prop_width, Ngrid_max); */
/*     if(iprop != i || jprop != j) break; */
/*   } */
/*   int x_jumps = 0; // 0 or 1 */
/*   int y_jumps = 0; // 0 or 1 */
/*   double p_of_proposed_state = f_1dim( (double)iprop/(double)Ngrid_max ) * f_1dim( (double)jprop/(double)Ngrid_max ); */
/*   double p_of_current_state = the_state->prob; */
/*   // prop ratio is 1 for symmetric proposal */
/*   if((p_of_proposed_state > p_of_current_state) || (p_of_proposed_state > drand() * p_of_current_state)){ // accept */
/*     // count if jumped from < 0.5 to > 0.5 in x and y */
/*     int h0max = Ngrid_max/2; */
/*     if( ( the_state->point[0] <= h0max  && iprop > h0max ) || ( the_state->point[0] > h0max  &&  iprop <= h0max ) ){  */
/*       x_jumps++; */
/*     } */
/*     if( ( the_state->point[1] <= h0max  && jprop > h0max ) || ( the_state->point[1] > h0max  &&  jprop <= h0max ) ){  */
/*       y_jumps++; */
/*     } */
/*     // set state to the proposed state: */
/*     the_state->point[0] = iprop;  */
/*     the_state->point[1] = jprop; */
/*     the_state->prob = p_of_proposed_state; */
/*     n_accept++; */
/*     // printf("%i %i   %i %i   %i %i \n", the_state->point[0], the_state->point[1], iprop, jprop, x_jumps, y_jumps); */
/*   }else{  */
/*     // reject proposed move; do nothing */
/*     n_reject++; */
/*   } */
/*   // ((int*)hist_ndim->bin_count[the_state->point[0]].array)[the_state->point[1]]++; */
/*   int* bp = get_bin_pointer(hist_ndim->bin_count, the_state->point);  */
/*   (*bp)++; */
/*   hist_ndim->total_weight++; */
/*   //    printf("%8i %8i %10.5g \n", the_state->point[0], the_state->point[1], the_state->prob); */
/*   return x_jumps + y_jumps; */
/* } */


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


/* Histogram_2d* initialize_histogram_2d(int Ngrid_max){ */
/*   Histogram_2d* histogram = (Histogram_2d*)malloc(1*sizeof(Histogram_2d)); */
/*   histogram->Ngrid_max = Ngrid_max; */
/*   histogram->total_weight = 0; */
/*   int** bin_count = (int**)malloc( (Ngrid_max+1) * sizeof(double*) ); */
/*   for(int i=0; i<=Ngrid_max; i++){ */
/*     bin_count[i] = (int*)malloc( (Ngrid_max+1) * sizeof(double) ); */
/*   } */

/*   for(int i=0; i<=Ngrid_max; i++){ */
/*     for(int j=0; j<=Ngrid_max; j++){ */
/*       bin_count[i][j] = 0; */
/*     } */
/*   } */
/*   histogram->bin_count = bin_count; */
/*   return histogram; */
/* } */
