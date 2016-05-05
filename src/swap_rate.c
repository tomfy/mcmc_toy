#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

double draw_ndim_gaussian(gsl_rng* rng, int n_dim, double sigma); // , double* rsq);
double avg_swap_Pa_gaussian(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps);
// double draw_ndim_exponential(gsl_rng* rng, int n_dim, double scale, double* xsum);
double draw_ndim_exponential(gsl_rng* rng, int n_dim, double scale);
double avg_swap_Pa_exponential(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps);

int main(int argc, char* argv[]){

  long n_reps = 100;
  if(argc > 1){
  n_reps = atoi(argv[1]);
  }
 int seed = 1234577;
  if(argc > 2){
    seed = atoi(argv[2]);
  }

 
  gsl_rng_env_setup();
  const gsl_rng_type* rng_type = gsl_rng_default;
  gsl_rng* the_rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(the_rng, seed);

  double desired_swap_rate = 0.234; 
  int min_n_dim = 1;
  int max_n_dim = 1;

for(int n_dimensions = min_n_dim; n_dimensions <= max_n_dim; n_dimensions *= 2)
{
    double Tfactor = pow(2.0,0.05);
    double T23_predicted = 1.0 + 2.95/(pow((double)n_dimensions, 0.53));
    for(int k=0; k<32; k++){ // (int k=0; k<=32; k++){
      //  double T_hot = 3.0;
        double T_hot = pow(T23_predicted, 0.1 + k*0.1);
      double avg_swap_rate_gauss = avg_swap_Pa_gaussian(the_rng, n_dimensions, T_hot, n_reps);
      double avg_swap_rate_exp = avg_swap_Pa_exponential(the_rng, n_dimensions, T_hot, n_reps);
      double delta_beta_sq = pow((1.0 - 1.0/T_hot), 2);
      double log_Tratio_sq = pow(log(T_hot), 2);
      printf("%6d  %8g   ", n_dimensions, T_hot);
      printf("%8g  %8g  %8g    ", avg_swap_rate_gauss, avg_swap_rate_gauss*delta_beta_sq,  avg_swap_rate_gauss*log_Tratio_sq);
      printf("%8g  %8g  %8g\n", avg_swap_rate_exp, avg_swap_rate_exp*delta_beta_sq, avg_swap_rate_exp*log_Tratio_sq);
    }
  printf("\n");    
 }
}

double draw_ndim_gaussian(gsl_rng* rng, int n_dim, double sigma){ // , double* rsq){
  double n_dim_density = 1.0;
  double lrsq = 0.0;
  for(int i=0; i<n_dim; i++){
    double x = gsl_ran_gaussian(rng, sigma);
    double X = x/sigma;
    lrsq += x*x;
    double density = exp(-0.5*X*X);
    n_dim_density *= density;
  }
  //  *rsq = lrsq;
  return lrsq; // n_dim_density;
}

double avg_swap_Pa_gaussian(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps){ 
  double sigma_cold = 1.0;
  double sigma_hot = sigma_cold*sqrt(T_hot);
  double var_cold = sigma_cold*sigma_cold;
  double var_hot = sigma_hot*sigma_hot;
  double* cold_rsqs = (double*)malloc(n_reps*sizeof(double));
  double* hot_rsqs = (double*)malloc(n_reps*sizeof(double));
 
  for(int i=0; i<n_reps; i++){
    cold_rsqs[i] = draw_ndim_gaussian(the_rng, n_dimensions, sigma_cold); // , cold_rsqs+i); // density_cold;
  }
  for(int i=0; i<n_reps; i++){
    hot_rsqs[i] =  draw_ndim_gaussian(the_rng, n_dimensions, sigma_hot); // , hot_rsqs+i); // density_hot;
  }

  double xxx = -0.5*(1.0/var_hot - 1.0/var_cold);
  double Pa_sum = 0.0;
  double Pa;
  for(int i=0; i<n_reps; i++){
    //   for(int j=0; j<n_reps; j++){
      double logPa = xxx*(cold_rsqs[i] - hot_rsqs[i]);
      if(logPa > 0.0){ Pa = 1; }else{ Pa = exp(logPa); }
      Pa_sum += Pa;
      // }
  }
  free(cold_rsqs);
  free(hot_rsqs);
  return Pa_sum/(double)(n_reps) ;
}

double draw_ndim_exponential(gsl_rng* rng, int n_dim, double scale){
  double n_dim_density = 1.0;
  double xsum = 0;
  // exponential
  for(int i=0; i<n_dim; i++){
    double x = gsl_ran_exponential(rng, scale);
    xsum += x;
    double density = exp(-1.0*x/scale);
    n_dim_density *= density;
  }
  return xsum;
}

double avg_swap_Pa_exponential(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps){ 
  double T_cold = 1.0;
  double scale_cold = 1.0;
  double scale_hot = scale_cold*T_hot;
  double* cold_xsums = (double*)malloc(n_reps*sizeof(double));
  double* hot_xsums = (double*)malloc(n_reps*sizeof(double));
 
  for(long i=0; i<n_reps; i++){
    cold_xsums[i] = draw_ndim_exponential(the_rng, n_dimensions, scale_cold); //, cold_xsums+i); // density_cold;
  }
  for(long i=0; i<n_reps; i++){
    hot_xsums[i] = draw_ndim_exponential(the_rng, n_dimensions, scale_hot); //, hot_xsums+i); // density_hot;
  }

  double xxx = (1.0/T_cold - 1.0/T_hot); 
  double Pa_sum = 0.0;
  double Pa;
  for(long i=0; i<n_reps; i++){
    // for(int j=0; j<n_reps; j++){
      double logPa = xxx*(cold_xsums[i] - hot_xsums[n_reps - i - 1]);
      if(logPa > 0.0){ Pa = 1.0; }else{ Pa = exp(logPa); }
      Pa_sum += Pa;
      //   printf("%5d  %8g  %8g  %8g\n", (int)(i+1), Pa, Pa_sum, Pa_sum/(double)(i+1)); 
      //  printf("%8g   %8g \n", rsq_cold, density_cold);
      //  }
  }
  // printf("Avg Pa: %8g \n", Pa_sum/(double)n_reps);
  free(cold_xsums);
  free(hot_xsums);
  double Pa_avg = Pa_sum/(double)n_reps;
  return Pa_avg;
}









/* double avg_swap_Pa(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps){  */
/*   // int n_dimensions = 2; */
/*   double sigma_cold = 1.0; */
/*   double sigma_hot = sigma_cold/sqrt(T_hot); */
/*   double var_cold = sigma_cold*sigma_cold; */
/*   double var_hot = sigma_hot*sigma_hot; */
/*   double rsq_cold; */
/*   double rsq_hot; */
/*   double density_cold; */
/*   double density_hot; */

/*   double Pa_sum = 0.0; */
/*   for(int i=0; i<n_reps; i++){ */
/*     density_cold = draw_ndim(the_rng, n_dimensions, sigma_cold, &rsq_cold); */
/*     density_hot = draw_ndim(the_rng, n_dimensions, sigma_hot, &rsq_hot); */

/*     double Pa = exp(-0.5*((rsq_cold - rsq_hot)/var_hot + (rsq_hot - rsq_cold)/var_cold ) ); */
/*     if(Pa > 1){ Pa = 1; } */
/*     Pa_sum += Pa; */
/*     //  printf("%8g   %8g \n", rsq_cold, density_cold); */
/*   } */
/*   // printf("Avg Pa: %8g \n", Pa_sum/(double)n_reps); */
/*   return Pa_sum/(double)n_reps ; */
/* } */

