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

my $n_generations = 100000;
my $n_burn_in = 0;
my $n_replicates = 1; # > 1 not implemented 
my $n_dimensions = 1;
my $peaks_string = '-0.5, 0.05, 0.5; 0.5, 0.05, 0.5';
my $n_peaks;
my $temperatures_string = '1.0';
my $proposals_string = 'gaussian, 0.119, 1.4, 0.7'; # if multiple proposals, separate with ; e.g. 'gaussian 0.07 1.5 0.9; gaussian ... '
my %n_bins_hash = (1 => 4096, 2 => 64, 3 => 16, 4 => 8);
my $n_bins = $n_bins_hash{$n_dimensions};
GetOptions(
	   'n_dimensions=i' => \$n_dimensions,
	   'peaks=s' => \$peaks_string, # e.g. '-0.5,0.05,0.5;0.5,0.05,0.5'  -> 2 peaks
	   'temperatures=s' => \$temperatures_string, # e.g. '1.0, 1.1, 1.2, 2.0';
	   'proposals=s' => \$proposals_string, # e.g. 'ball,0.1,1.5,0.9'
	   'n_bins=i', => \$n_bins,
	   'n_generations=i', => \$n_generations,
	   'n_burnin=i', => \$n_burn_in,
);


my @peak_strings = split(";", $peaks_string);
$n_peaks = scalar @peak_strings;
for(@peak_strings){
  s/,/ /g;
}
# print "# proposals string: $proposals_string \n";

$temperatures_string =~ s/[,;]/ /g;
my @temperatures = split(" ", $temperatures_string);
my $n_temperatures = scalar @temperatures;
my @proposal_strings = split(";", $proposals_string);
for(@proposal_strings){
  s/[,]/ /g;
}
while(scalar @proposal_strings < scalar @temperatures){
  push @proposal_strings, $proposal_strings[scalar @proposal_strings -1]; # duplicate last proposal to get enough
}

my $the_arg_string = "$n_dimensions  $n_peaks  ";
$the_arg_string .= join("  ", @peak_strings) . "  ";
$the_arg_string .= "$n_burn_in $n_generations $n_replicates ";
$the_arg_string .= "$n_temperatures ";
for my $i (0..$n_temperatures-1){
my $T = $temperatures[$i];
  $the_arg_string .= "$T " . $proposal_strings[$i] . "  ";
}
$the_arg_string .= "$n_bins";

print "# $the_arg_string \n";

# my $arg_string = "3   2   -0.5 0.06 0.5   0.5 0.06 0.5   "; # n_dim, n_modes, targ mean, sigma, weight, ...
# $arg_string .= "0 2000000 1   "; # burn-in steps, post-burn-in steps, replicates
# #$arg_string .= "2   1.0 ball 0.06 1.0 0.9   2.0 ball 0.06 1.0 0.5   "; # n_temperatures, ...
# $arg_string .= "1   1.0  ball 0.1 1.5 0.7   "; # n_temperatures, then T_i, proposal_i
# $arg_string .= "25";  # n_bins

# my $mcmc_toy_out = `~/mcmc_toy/src/mcmc_toy $arg_string `;
system "~/mcmc_toy/src/mcmc_toy $the_arg_string";
# print "$mcmc_toy_out \n";



