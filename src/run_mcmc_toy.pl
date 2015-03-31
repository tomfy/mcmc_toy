#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# target distribution parameters:
# n_dimensions n_modes x1 sigma1 w1 x1 sigma2 w2 ...

# burn_in_steps mcmc_steps n_replicates
# n_temperatures T_i  shape_i w1_i w2_i p1_i (T and proposals for each T) 

# analysis parameters
# n_bins



# GetOptions(
# 	   'n_dimensions=s'           => \$gg_filename, #
# 	   'nj!'          => \$do_nj, # whether to do NJ tree construction for actual data
# 	   #  exclamation point means can use  -nj  and  -nonj
# 	   'ft!'      => \$do_ft, # whether to do FastTree for actual data
# 	   'phyml!' => \$do_phyml, # whether to do Phyml tree construction for actual data
# 	   'support!'          => \$support, # -nosupport to turn off outputting of local branch support numbers.
# 	   'n_nj_bs=i' => \$n_nj_bs, # number of NJ bootstraps to do
# 	   'n_ft_bs=i' => \$n_ft_bs, # number of FastTree bootstraps to do
# 	   'n_phyml_bs=i' => \$n_phyml_bs, # number of Phyml bootstraps to do
# 	   #    'ml_bs!' => \$ml_bs, # boolean to do ML bs or not (default: 0)
# 	   'species_tree=s' => \$species_tree_newick_file,

my $arg_string = "3   2   -0.5 0.06 0.5   0.5 0.06 0.5   "; # n_dim, n_modes, targ mean, sigma, weight, ...
$arg_string .= "1000 40000 1   "; # burn-in steps, post-burn-in steps, replicates
#$arg_string .= "2   1.0 0.06 1.0 0.9   2.0 0.06 1.0 0.5   "; # n_temperatures, ...
$arg_string .= "1   1.0  ball 0.1 1.5 0.7   "; # n_temperatures, then T_i, proposal_i
$arg_string .= "24";  # n_bins

# my $mcmc_toy_out = `~/mcmc_toy/src/mcmc_toy $arg_string `;
system "~/mcmc_toy/src/mcmc_toy $arg_string";
# print "$mcmc_toy_out \n";



