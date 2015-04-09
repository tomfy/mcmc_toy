#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>

double avg_deltasquared(int lag, const double* a, int n);
double avg_xy(const double* a1, const double* a2, int n);
double integrated_corr(const double* a, const int n, const int start_lag, const int end_lag, const double mu);
double integrated_corrzz(const double* a, int n, int start_lag, int end_lag, double mu);
double corr_x_var(const double* a, int lag, int n, double mu_sq);
//int first_negative_corr_lag(const double* a, int n, double mu_sq);
int first_negative_corr_lag_old(int n, const double* a, double mu_sq);
double avg_x(int n, const double* a);
int first_negative_corr_lag(int n, const double* a, double mu_sq, double sigma_sq);

int main(int argc, char* argv[]){
  int chosen_col = 3;
  int max_n_input_lines = 100000;
  double target_mean = 0.0;
  double target_variance = 0.0;
  int max_lag = -1;
  if(argc >= 2){
    chosen_col = atoi(argv[1]);
  }
  if(argc >= 3){
    max_n_input_lines = atoi(argv[2]);
  }
  if(argc >= 4){
    target_mean = atof(argv[3]);
  }
  double mu_sq = target_mean*target_mean;
  if(argc >= 5){
    target_variance = atof(argv[4]);
  }
  if(argc >= 6){
    max_lag = atoi(argv[5]);
  }

  int array_size = max_n_input_lines + 2;
  int bytes_read;
  int nbytes = 200;
  char* a_string = (char*)malloc((nbytes+1)*sizeof(char));
  char* token;
  int count_col = 0;
  double* xs = (double*)malloc(array_size*sizeof(double));
  int array_index = 0;
  while(1){
    bytes_read = getline(&a_string, (size_t*)&nbytes, stdin);
    //  printf("String: %s \n", a_string);
    count_col = 0;
    if(bytes_read < 0){ break; } // no bytes read - done.
    token = strtok(a_string, " ");
    if(token[0] == '#' || token[0] == '\n'){ continue; } // skip comment line
    if(count_col == chosen_col){
      //   printf("%i %i %s \n", array_index, count_col, token);
      double x = atof(token);
      xs[array_index] = x;
      array_index++;
      continue;
    }
    while( token != NULL ){
      //   printf("[%s]\n", token);
      count_col++;
      token = strtok(NULL, " ");
      if(token[0] == '\n'){ break; }
      //   printf("%i %i %s \n", array_index, count_col, token);
      if(count_col == chosen_col){
	double x = atof(token);
	xs[array_index] = x;
	array_index++;
	break;
      }
    }
    if(array_index >= array_size){ break; }
  }
  int n_data_points = array_index;
  printf("# array size: %i \n", array_index);
  if(max_lag < 0){
    max_lag = n_data_points/2;
  }

  for(int lag = 1; lag < max_lag; lag += 1000){
    double corr = corr_x_var(xs, lag, n_data_points-max_lag, mu_sq);
    corr /= target_variance;
    printf("lag, corr: %i  %g \n", lag, corr);
  }
 
  int fncl = first_negative_corr_lag(n_data_points, xs, mu_sq, target_variance);
  fncl *= 1;
  printf("fncl: %i \n", fncl);
  // double sum_corr = integrated_corr(xs, n_data_points-3*fncl, 1, fncl, target_mean)/target_variance;
  int n_chunks = 1;
  for(int i=0; i<n_chunks; i++){
    double sum_corri =  integrated_corr(xs, n_data_points-n_chunks*fncl, i*fncl+1, (i+1)*fncl, target_mean)/target_variance;
    // double sum_corr3 =  integrated_corr(xs, n_data_points-3*fncl, 2*fncl+1, 3*fncl, target_mean)/target_variance;
    // double sum_corr4 =   integrated_corr(xs, n_data_points-3*fncl, 1, 3*fncl, target_mean)/target_variance;
    // double efficiency = 1.0/(target_variance * (1.0 + 2.0*sum_corr));  
    //  double eff2 = 1.0/(target_variance * (1.0 + 2.0*scorrz1));
    // printf("#sum_corr, eff: %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n", sum_corr, efficiency, 1.0/efficiency, sum_corr2, sum_corr3, sum_corr4);
    printf("%i %10.5g\n", i, sum_corri);
  }
}

double avg_deltasquared(int lag, const double* a, int n){
  double sum = 0.0;
  for(int i=0; i<n; i++){
    double dx =  a[i] - a[i+lag];
    sum += dx*dx;
  }
  return sum/(double)n;
}

double avg_xy(const double* a1, const double* a2, int n){
  double sum = 0.0;
  for(int i=0; i<n; i++){
    sum += a1[i]*a2[i];
    //   printf("%i %g %g %g  %g \n", i, a1[i], a2[i], a1[i]*a2[i], sum);
   }
  return sum/(double)n;
}

double corr_x_var(const double* a, int lag, int n, double mu_sq){
  return avg_xy(a, a+lag, n) - mu_sq;
}

double avg_x(int n, const double* a){
  double sum = 0.0;
  for(int i=0; i<n; i++){
    sum += a[i];
  }
  return sum/(double)n;
}

int first_negative_corr_lag(int n, const double* a, double mu_sq, double sigma_sq){
  double factor = 1.1;
  int min_increment = 4;
  double corr_threshold = 0.0;
  int max_lag = n/2;
  double sum_corr = 0.0;
  int prev_lag, lag;
  double prev_corr;
  int fncl;
  for(lag=1; lag<=max_lag; lag = (int)((lag+min_increment)*factor) ){
    double corr = corr_x_var(a, lag, n-max_lag, mu_sq)/sigma_sq; 
    if(lag>1){ sum_corr += 0.5*(lag - prev_lag)*(corr + prev_corr); }
    printf("%6i %8.5f %8.5g\n", lag, corr, sum_corr);       
    if(corr < corr_threshold){ 
      sum_corr += 0.5*corr*corr*(lag - prev_lag)/(prev_corr - corr);
      fncl = lag + (int)(corr*(lag - prev_lag)/(prev_corr - corr)) + 1;
      printf("fncl: %i %i %g %g \n", lag, fncl, corr, sum_corr); 
      break; 
    }
    prev_lag = lag;
    prev_corr = corr;
  }
  int ub = lag;
  int lb = prev_lag;
 
  while(ub-lb > 1){
     double fncl_corr = corr_x_var(a, fncl, n-max_lag, mu_sq)/sigma_sq;
    if(fncl_corr < corr_threshold){
      ub = fncl;
    }else{
      lb = fncl;
    }
    fncl = (ub+lb)/2;
  }
  return ub;
}

int first_negative_corr_lag_old(int n, const double* a, double mu_sq){
  // n is <= size of array a
  // will look for smallest lag which gives corr <= 0
  double factor = 1.1;
  int lag_lb = 1;
  int lag_ub;
  int plb = 1;
  int lb = 2;
  int max_lag = n/2;
  double cxv_plb = corr_x_var(a, plb, n-plb, mu_sq);
  if(cxv_plb <= 0){ return plb; }
  double cxv_lb = corr_x_var(a, lb, n-lb, mu_sq);
  if(cxv_lb <= 0){ return lb; }
  int ub = -1;
  double cxv_ub = -1000;
  double next_cxv;
  int do_newton = 1;
  while(1){
    int next_x;  
    if(do_newton){ // no upper bound yet;
      double slope = (cxv_lb - cxv_plb)/(lb-plb);
      next_x = (int)((plb - cxv_plb/slope) + 0.5);  // newton's method, round to int     
      if(next_x > factor*(lb+3) || next_x <= lb){ next_x = (int)(factor*(lb+3)); }
      if(next_x > max_lag){ next_x = max_lag; }
      next_cxv = corr_x_var(a, next_x, n-next_x, mu_sq);
          printf("# a next_x: %i %g\n", next_x, next_cxv);
      if(next_cxv <= 0){
	ub = next_x;
	do_newton = 0;
      }else{
	if(next_x > lb){
	  plb = lb;
	  cxv_plb = cxv_lb;
	  lb = next_x;
	  cxv_lb = next_cxv;	  
	}
      }
    }else{ // upper bound exists
      next_x = (4*lb + ub)/5;
      if(next_x == lb){ next_x++; }
      if(next_x > max_lag){ next_x = max_lag; }
      next_cxv = corr_x_var(a, next_x, n-next_x, mu_sq);
          printf("# b next_x: %i %g\n", next_x, next_cxv);
      if(next_cxv <= 0){
	if(next_x < ub){
	  ub = next_x;
	  cxv_ub = next_cxv;
	}
      }else{
	if(next_x > lb){
	  plb = lb;
	  cxv_plb = cxv_lb;
	  lb = next_x;
	  cxv_lb = next_cxv;	  
	}
      }
    }
    if(ub > 0  &&  ub - lb <= 1){ return ub; }
  } // end of while()
}
 
    

    

double integrated_corr(const double* a, const int n, const int start_lag, const int end_lag, const double mu){
  // n is number of points to include in each avg.
  // k is max lag
  // mu is mean
  double sum = 0.0; // -1.0*k*mu;
  for(int i=start_lag; i<=end_lag; i++){
    sum += a[i] - mu;
  }
  //  printf("j=0 sum over lags: %10.5g\n", sum); 
  double result = 0.0;
  for(int j=0; j<n; j++){
  result += (a[j] - mu)*sum;
  //  printf("j= %i sumoverlags:  %10.5g %10.5g %10.5g \n", j, sum, a[j]*sum, result/(j+1)); 
  sum += a[j+end_lag+1] - a[j+start_lag];
 
}
  result /= (double)n;
  return result;
}

double integrated_corrzz(const double* a, int n, int start_lag, int end_lag, double mu){
  // n is number of points to include in each avg.
  // k is max lag
  // mu is mean
  double sum = 0.0; // -1.0*k*mu;
  for(int i=start_lag; i<=end_lag; i++){
    sum += a[i];
    //   printf("i: %i  %10.5g  ", i, a[i]);
  } //printf("\n");
  //  printf("j=0 sum over lags: %10.5g\n", sum); 
  double result = 0.0;
  for(int j=0; j<n; j++){
    //  printf("j: %i %10.5g  %10.5g\n", j, a[j], (a[j])*sum);
  result += (a[j])*sum;
  //  printf("j= %i sum over lags: %10.5g %10.5g %10.5g \n", j, sum, a[j]*sum, result/(j+1)); 
  sum += a[j+end_lag+1] - a[j+start_lag];
 
}
  result /= (double)n;
  result -= mu*mu;
  return result;
}





 /* int iprev; */
  /* double prev_corr; */
  /* int fncl1; */
  /* for(int i=1; i<=max_lag; i = (int)((i+4)*1.1) ){ */
  /*   double corr = corr_x_var(xs, i, n_data_points-max_lag, mu_sq)/target_variance;  */
  /*   if(i>1){ sum_corr += 0.5*(i - iprev)*(corr + prev_corr); } */
  /*   printf("%6i %8.5f %8.5g\n", i, corr, sum_corr);   */
      
  /*   if(corr < 0.0){  */

  /*     sum_corr += 0.5*corr*corr*(i - iprev)/(prev_corr - corr); */
  /*     fncl1 = i + (int)(corr*(i - iprev)/(prev_corr - corr)) + 1; */
  /*     break; } */
  /*   iprev = i; */
  /*   prev_corr = corr; */
  /* } */
  
  //  int fncl2 = first_negative_corr_lag(n_data_points, xs, mu_sq);
  // printf("#first negative corr lags: %i\n", fncl1);
  //  double scorrz1 = integrated_corrz(xs, n_data_points-fncl1, fncl1, target_mean)/target_variance;
  //int fncl(int n, const double* a, double mu_sq, double sigma_sq){
