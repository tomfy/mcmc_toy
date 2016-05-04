#include <stdlib.h>
#include <math.h>
// #include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

double draw_ndim(gsl_rng* rng, int n_dim, double sigma, double* rsq);
double avg_swap_Pa(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps);
double avg_swap_Pa_fast(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps);

int main(int argc, char* argv[]){

  long n_reps = 100;
  if(argc > 1){
  n_reps = atoi(argv[1]);
  }

  int seed = 1234577;
  gsl_rng_env_setup();
  const gsl_rng_type* rng_type = gsl_rng_default;
  gsl_rng* the_rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(the_rng, seed);

  //  long n_reps = 1000;
  double desired_swap_rate = 0.234; 
  int max_n_dim = 2*1025;

for(int n_dimensions = 1; n_dimensions <= max_n_dim; n_dimensions *= 2)
{
    double T_hot_prev = 0;
    double asr = 1.0;
    double asr_prev = 2.0;
    double Tfactor = pow(2.0,0.05);
    double T23_predicted = 1.0 + 2.95/(pow((double)n_dimensions, 0.53));
      for(int k=0; k<=25; k++){ // double T_hot = 1; // + 2.9*pow((double)n_dimensions, -0.55); // 0.8 * exp(1.0)/sqrt((double)n_dimensions); 
        double T_hot = pow(T23_predicted, 0.5 + k*0.05);
        //  T_hot < 100; 
        //  T_hot *= Tfactor){
      double avg_swap_rate = avg_swap_Pa_fast(the_rng, n_dimensions, T_hot, n_reps);
      printf("%6d  %8g   %8g  %8g\n", n_dimensions, T_hot, avg_swap_rate, avg_swap_rate*pow((1.0 - 1.0/T_hot),2));
      asr = avg_swap_rate;
     /*  if(avg_swap_rate < 0.234){  */
/*         double lnT23 = log(T_hot_prev) + (log(T_hot) - log(T_hot_prev))/(log(asr) - log(asr_prev))*(log(desired_swap_rate) - log(asr_prev)); */
/*         double T23 = exp(lnT23); */
/*         avg_swap_rate =  avg_swap_Pa(the_rng, n_dimensions, T23, n_reps); */
/* printf("%6d  %8g   %8g  %8g\n", n_dimensions, T23, avg_swap_rate, avg_swap_rate*(1.0 - 1.0/T_hot)); */
/* //     printf("%6d  %8g   %8g  %6g %6g %6g %6g\n", n_dimensions, T23, avg_swap_rate, T_hot_prev, asr_prev, T_hot, asr); */
/*         break; } */
      T_hot_prev = T_hot;
      asr_prev = asr;
      
    }
  printf("\n");    
 }
// printf("\n");
}


double avg_swap_Pa(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps){ 
  // int n_dimensions = 2;
  double sigma_cold = 1.0;
  double sigma_hot = sigma_cold/sqrt(T_hot);
  double var_cold = sigma_cold*sigma_cold;
  double var_hot = sigma_hot*sigma_hot;
  double rsq_cold;
  double rsq_hot;
  double density_cold;
  double density_hot;

  double Pa_sum = 0.0;
  for(int i=0; i<n_reps; i++){
    density_cold = draw_ndim(the_rng, n_dimensions, sigma_cold, &rsq_cold);
    density_hot = draw_ndim(the_rng, n_dimensions, sigma_hot, &rsq_hot);

    double Pa = exp(-0.5*((rsq_cold - rsq_hot)/var_hot + (rsq_hot - rsq_cold)/var_cold ) );
    if(Pa > 1){ Pa = 1; }
    Pa_sum += Pa;
    //  printf("%8g   %8g \n", rsq_cold, density_cold);
  }
  // printf("Avg Pa: %8g \n", Pa_sum/(double)n_reps);
  return Pa_sum/(double)n_reps ;
}

double avg_swap_Pa_fast(gsl_rng* the_rng, int n_dimensions, double T_hot, long n_reps){ 
  // int n_dimensions = 2;
  double sigma_cold = 1.0;
  double sigma_hot = sigma_cold/sqrt(T_hot);
  double var_cold = sigma_cold*sigma_cold;
  double var_hot = sigma_hot*sigma_hot;
  double* cold_densities = (double*)malloc(n_reps*sizeof(double));
  double* cold_rsqs = (double*)malloc(n_reps*sizeof(double));
  double* hot_densities = (double*)malloc(n_reps*sizeof(double));
  double* hot_rsqs = (double*)malloc(n_reps*sizeof(double));
 
  for(int i=0; i<n_reps; i++){
    cold_densities[i] = draw_ndim(the_rng, n_dimensions, sigma_cold, cold_rsqs+i); // density_cold;
    hot_densities[i] =  draw_ndim(the_rng, n_dimensions, sigma_hot, hot_rsqs+i); // density_hot;
  }

  double xxx = -0.5*(1.0/var_hot - 1.0/var_cold);
  double Pa_sum = 0.0;
  double Pa;
  for(int i=0; i<n_reps; i++){
    for(int j=0; j<n_reps; j++){
      double logPa = xxx*(cold_rsqs[i] - hot_rsqs[j]);
      if(logPa > 0.0){ Pa = 1; }else{ Pa = exp(logPa); }
      Pa_sum += Pa;
      //  printf("%8g   %8g \n", rsq_cold, density_cold);
    }
  }
  // printf("Avg Pa: %8g \n", Pa_sum/(double)n_reps);
  free(cold_densities);
  free(hot_densities);
  free(cold_rsqs);
  free(hot_rsqs);
  return Pa_sum/(double)(n_reps*n_reps) ;
}

double draw_ndim(gsl_rng* rng, int n_dim, double sigma, double* rsq){
  // double* xs = (double*)malloc(n_dim*sizeof(double));
  double n_dim_density = 1.0;
  double lrsq = 0.0;
  for(int i=0; i<n_dim; i++){
    double x = gsl_ran_gaussian(rng, sigma);
    double X = x/sigma;
    lrsq += x*x;
    double density = exp(-0.5*X*X);
    n_dim_density *= density;
  }
  *rsq = lrsq;
  return n_dim_density;
}


