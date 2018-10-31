#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "mcmc2body.h"
#include "target.h"
#include "chain_architecture.h"
#include "mcmc2body_structs.h"

// *****************************************************************************************************

chain_architecture* set_up_chain_architecture(int n_per_level, int symmetry, int n_levels, double Thot, double min_prop_w, double max_prop_w, double min_kernel_width, double max_kernel_width){
  chain_architecture* chain_arch = (chain_architecture*)malloc(sizeof(chain_architecture));
  chain_arch->symmetry = symmetry;
  chain_arch->n_per_level = n_per_level;
  chain_arch->n_levels = n_levels;
  chain_arch->inverse_Temperatures = (double*)malloc(n_levels*sizeof(double));
  chain_arch->Lprop_widths = (double*)malloc(n_levels*sizeof(double));
  chain_arch->Rprop_widths = (double*)malloc(n_levels*sizeof(double));
  chain_arch->kernel_widths = (double*)malloc(n_levels*sizeof(double));

  if(n_per_level == 1){ // standard mcmcmc.
    //  printf("xxxxx\n");
    chain_arch->inverse_Temperatures[0] = 1.0;
    double Tratio = 1.0/pow(Thot, 1.0/(n_levels-1)); // < 1
    chain_arch->Lprop_widths[0] = min_prop_w;
    printf("in set up ch arch: %g  %g \n", min_prop_w, Tratio);
    for(int i=1; i<n_levels; i++){
      chain_arch->inverse_Temperatures[i] = Tratio * chain_arch->inverse_Temperatures[i-1];
      //  printf("%i %g \n", i, chain_arch->inverse_Temperatures[i]);
      chain_arch->Lprop_widths[i] = chain_arch->Lprop_widths[0] / sqrt(chain_arch->inverse_Temperatures[i]);
      //  printf("in suca: %i  %20.12g \n", i, chain_arch->Lprop_widths[i]);
    }
  }else if(n_per_level == 2){ // ((1-e)*pix + e*piy)*K 
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


void print_chain_architecture_info(chain_architecture* arch){
  printf("# symmetric?: %2i  n_levels: %2i  n_per_level: %2i \n", arch->symmetry, arch->n_levels, arch->n_per_level);
  for(int i=0; i<arch->n_levels; i++){
    printf("%8.6g  %8.6g  %8.6g  %8.6g \n", arch->inverse_Temperatures[i], arch->Lprop_widths[i], arch->Rprop_widths[i], arch->kernel_widths[i]);
  }
}

// **************************************************************************************************
