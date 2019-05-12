#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "target.h"
#include "chain_architecture.h"

// *****************************************************************************************************

chain_architecture* set_up_chain_architecture( int n_levels, int n_per_level, int symmetry, double Thot, double min_prop_w, double max_prop_w, double min_kernel_width, double max_kernel_width, double two_body_interpolation_power){
  chain_architecture* chain_arch = (chain_architecture*)malloc(sizeof(chain_architecture));
  chain_arch->symmetry = symmetry;
  chain_arch->n_per_level = n_per_level;
  chain_arch->n_levels = n_levels;
  chain_arch->two_body_interpolation_power = two_body_interpolation_power; //
  chain_arch->inverse_Temperatures = (double*)calloc(n_levels, sizeof(double));
  chain_arch->Lprop_widths = (double*)malloc(n_levels*sizeof(double));
  chain_arch->Rprop_widths = (double*)malloc(n_levels*sizeof(double));
  chain_arch->kernel_widths = (double*)calloc(n_levels, sizeof(double));

  if(n_per_level == 1){ // standard mcmcmc.
    chain_arch->inverse_Temperatures[0] = 1.0;
    double Tratio = 1.0/pow(Thot, 1.0/(n_levels-1)); // < 1
    chain_arch->Lprop_widths[0] = min_prop_w;
    chain_arch->Rprop_widths[0] = min_prop_w;
    //  printf("in set up ch arch: %g  %g \n", min_prop_w, Tratio);
    for(int i=1; i<n_levels; i++){
      chain_arch->inverse_Temperatures[i] = Tratio * chain_arch->inverse_Temperatures[i-1];
      //  printf("%i %g \n", i, chain_arch->inverse_Temperatures[i]);
      chain_arch->Lprop_widths[i] = chain_arch->Lprop_widths[0] / sqrt(chain_arch->inverse_Temperatures[i]);
      chain_arch->Rprop_widths[i] = chain_arch->Rprop_widths[0] / sqrt(chain_arch->inverse_Temperatures[i]);
      //  printf("in suca: %i  %20.12g \n", i, chain_arch->Lprop_widths[i]);
    }
  }else if(n_per_level == 2){ 
    if(chain_arch->symmetry == 1){
      chain_arch->Lprop_widths[0] = min_prop_w;
      chain_arch->Rprop_widths[n_levels-1] = min_prop_w;
      chain_arch->kernel_widths[0] = max_kernel_width;
      double dpropw = (n_levels > 1)? (max_prop_w - min_prop_w)/(n_levels - 1) : 0.0;
      double propw_ratio = (n_levels > 1)? pow(max_prop_w/min_prop_w, 1.0/(n_levels - 1)) : 1.0;
      for(int i=1; i<n_levels; i++){
        chain_arch->Lprop_widths[i] = chain_arch->Lprop_widths[i-1] 
          //  + dpropw;
          * propw_ratio;
        //  printf("Lprop width: %i %g \n", i, chain_arch->Lprop_widths[i]);
        chain_arch->Rprop_widths[n_levels-1 - i] = chain_arch->Lprop_widths[i];
        chain_arch->kernel_widths[i] = max_kernel_width;
      }
    }else{ // asymmetrical

      printf("# Asymmetrical 2-body.\n"); 
      chain_arch->Lprop_widths[0] = min_prop_w;
      chain_arch->Rprop_widths[0] = min_prop_w;
      chain_arch->kernel_widths[0] = min_kernel_width;
      double propw_ratio = pow(max_prop_w/min_prop_w, 1.0/(n_levels - 1));
      double kernelw_ratio = pow(max_kernel_width/min_kernel_width, 1.0/(n_levels-1));
      for(int i=1; i<n_levels; i++){
        chain_arch->Lprop_widths[i] = min_prop_w;
        //    printf("Lprop width: %i %g \n", i, chain_arch->Lprop_widths[i]);
        chain_arch->Rprop_widths[i] = chain_arch->Rprop_widths[i-1] * propw_ratio;
        chain_arch->kernel_widths[i] = chain_arch->kernel_widths[i-1] * kernelw_ratio;
        //  printf("%i  %g\n", i, chain_arch->kernel_widths[i]);
      }
    }
  }else{ // pix*K_i(x,y)
    exit(1);
  }
  return chain_arch;
}


void print_chain_architecture_info(const chain_architecture* const arch){
  printf("# n_levels: %2i  n_per_level: %2i  symmetric?: %2i \n", arch->n_levels, arch->n_per_level, arch->symmetry);
  printf("#     Tinv     Lpropw     Rpropw    Kwidth \n");
  for(int i=0; i<arch->n_levels; i++){
    printf("# %8.6g  %8.5g  %8.6g  %8.5g  \n",  arch->inverse_Temperatures[i], 
           arch->Lprop_widths[i], arch->Rprop_widths[i], arch->kernel_widths[i]);
  }
}

// **************************************************************************************************

double Kernel(double Dsq, double Lsq){ // convolution kernel
  return exp(-0.5*Dsq/Lsq);
  //  return (Dsq >= Lsq)? 0.0 : 1.0 - Dsq/Lsq;
}

double PI(double pix, double piy, int it, int n_Ts, double K, const chain_architecture* const arch){
  if(arch->symmetry == 1){
    if(it == 0){
      return pix*K;
    }else{
      int itop = n_Ts - 1;
      double epsilon = (double)it/itop;
      if(it == itop){
        return piy*K;
      }else{
        if(arch->two_body_interpolation_power == 1){
          double pixpiylc = (1.0 - 1.0*it/itop)*pix + (1.0*it/itop)*piy;  // linear interpolation (equiv. to P = 1)
          return pixpiylc*K;
        }else if(arch->two_body_interpolation_power == 0){
          return pow(pix, (1.0 - 1.0*it/itop)) * pow(piy, 1.0*it/itop) * K;  // geometric interpolation (P -> 0 limit)
        }else{ // p = 0 at ends, arch->two_body_interpolation_power in middle
          if(0){
            double factor = (epsilon <= 0.5)? epsilon : 1.0 - epsilon;
            double power = arch->two_body_interpolation_power * factor;
            double result = pow( 
                                (1.0 - epsilon)*pow(pix,power)  +  (epsilon)*pow(piy, power)
                                , (1.0/power) ) * K; // weighted power mean
            //  printf("%g %g %g  %g\n", pix, piy, result, two_body_interpolation_power);
            return result;
          }else{ //
            double result = pow(pix, 1.0 - epsilon) * pow(piy, 1.0 - epsilon) * K;
            return result;
          }
        }
      }
    }
  }else{ // asymmetrical
    return pix*K;
  }
}
    
