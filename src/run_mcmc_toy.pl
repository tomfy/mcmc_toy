#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw ( min max sum );
use Statistics::Basic qw ( mean stddev );
my $epsilon = 1e-6;
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
my $n_replicates = 1;		# > 1 not implemented 
my $n_dimensions = 1;
my $peaks_string = '-0.5, 0.1, 0.5; 0.5, 0.1, 0.5';
my $n_peaks;
my $temperatures_string = '1.0';
my $proposals_string = 'gaussian, SIGMA1, 1.4, 1.0'; # if multiple proposals, separate with ; e.g. 'gaussian 0.07 1.5 0.9; gaussian ... '
my %n_bins_hash = (1 => 4096, 2 => 64, 3 => 16, 4 => 8);
my $n_bins = $n_bins_hash{$n_dimensions};
my $n_bins_1d = 216;
my $chain_type = "mcmc";	# or "iid"
GetOptions(
	   'n_dimensions=i' => \$n_dimensions,
	   'peaks=s' => \$peaks_string, # e.g. '-0.5,0.05,0.5;0.5,0.05,0.5'  -> 2 peaks
	   'temperatures=s' => \$temperatures_string, # e.g. '1.0, 1.1, 1.2, 2.0';
	   'proposals=s' => \$proposals_string, # e.g. 'ball,0.1,1.5,0.9'
	   'n_bins=i' => \$n_bins,
	   'n_1d_bins=i' => \$n_bins_1d,
	   'n_generations=i' => \$n_generations,
	   'n_burnin=i' => \$n_burn_in,
	   'type=s' => \$chain_type,
	  );


my @peak_strings = split(";", $peaks_string);
$n_peaks = scalar @peak_strings;
for (@peak_strings) {
  s/,/ /g;
}
# print "# proposals string: $proposals_string \n";

$temperatures_string =~ s/[,;]/ /g;
my @temperatures = split(" ", $temperatures_string);
my $n_temperatures = scalar @temperatures;
my @proposal_strings = split(";", $proposals_string);
for (@proposal_strings) {
  s/[,]/ /g;
}
while (scalar @proposal_strings < scalar @temperatures) {
  push @proposal_strings, $proposal_strings[scalar @proposal_strings -1]; # duplicate last proposal to get enough
}

my $proto_arg_string = "$n_dimensions  $n_peaks  ";
$proto_arg_string .= join("  ", @peak_strings) . "  ";
$proto_arg_string .= "$n_burn_in NGENERATIONS $n_replicates ";
$proto_arg_string .= "$n_temperatures ";
for my $i (0..$n_temperatures-1) {
  my $T = $temperatures[$i];
  $proto_arg_string .= "$T " . $proposal_strings[$i] . "  ";
}
$proto_arg_string .= "$n_bins  $n_bins_1d  ";
#print "$n_generations ;;; $proto_arg_string \n"; exit;
# my $seed = srand();
# $proto_arg_string .= "$seed ";

# print "# $proto_arg_string \n";

# my $arg_string = "3   2   -0.5 0.06 0.5   0.5 0.06 0.5   "; # n_dim, n_modes, targ mean, sigma, weight, ...
# $arg_string .= "0 2000000 1   "; # burn-in steps, post-burn-in steps, replicates
# #$arg_string .= "2   1.0 ball 0.06 1.0 0.9   2.0 ball 0.06 1.0 0.5   "; # n_temperatures, ...
# $arg_string .= "1   1.0  ball 0.1 1.5 0.7   "; # n_temperatures, then T_i, proposal_i
# $arg_string .= "25";  # n_bins

my $tvd_limit = 0.08 ;
my $n_reps = 16;
srand();

my $sig1 = 0.75;
my $the_arg_string;
while (1) {
  $the_arg_string = $proto_arg_string;
  $the_arg_string =~ s/NGENERATIONS/$n_generations/;
  $the_arg_string =~ s/SIGMA1/$sig1/;
  my $seed = int(rand(1000000));
  my $mcmc_toy_out =  `~/mcmc_toy/src/mcmc_toy $the_arg_string $seed mcmc`;
  my $last_line = `tail -1 tvd_vs_gen`;
  #exit;
  my @tvds = split(" ", $last_line);
  shift @tvds;
  my @mcmc_1d_tvds = @tvds[3 .. 3+$n_dimensions-1];
  my $avg_1d_tvd = mean(@mcmc_1d_tvds);
  last if($avg_1d_tvd < $tvd_limit);
  $n_generations = int( 1.2*$n_generations*($avg_1d_tvd/$tvd_limit)**2 );
}
#exit;
my @sig_1s = (0.14, 0.20, 0.32, 0.5, 0.7, 1.0, 1.4, 2.0, 4.0, 7.0);
# @sig_1s = (0.2, 0.45, 0.7, 1.0, 1.4);
# for ($sig1 = ; $sig1 < 5.0; $sig1 *= 1.3)
for my $sig_1 (@sig_1s) {
  trials($proto_arg_string, $n_reps, $n_generations, $sig_1, $n_dimensions);
}

sub trials{
  my $proto_arg_string = shift;
  my $n_reps = shift;
  my $n_generations = shift;
  my $sig1 = shift;
  my $n_dimensions = shift;
  my $the_arg_string = $proto_arg_string;
  $the_arg_string =~ s/NGENERATIONS/$n_generations/;
  $the_arg_string =~ s/SIGMA1/$sig1/;

  my @avg_effs = (0, 0, 0, 0);
  my @ndim_tvds = ();
  my @orthant_tvds = ();
  my @refl_tvds = ();
  my @onedim_tvds = ();
  my @ms_dmus = ();
  my @ksds = ();
  my @mean_qs = ();
  my @mu = (0,0,0,0,0,0,0,0,0,0,0);
  my %dmus = (); # $dmus{$i} is array_ref to array of dmus for component $i
  for (0..$n_dimensions-1) {
    $dmus{$_} = []; $mu[$_] = 0;
  }
  for (1..$n_reps) {
    my $seed = int(rand(1000000));
    system "~/mcmc_toy/src/mcmc_toy  $the_arg_string  $seed  'mcmc' >  mcmc_toy_samples_tmp";
# my $run_params = `cat run_params`;
# my @runparam_lines = split("\n", $run_params);
# my ($T, $shape, $s1, $s2, $p1);
#     for(@runparam_lines){
#       if(/T, Proposal\S+:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
# 	($T, $shape, $s1, $s2, $p1) = ($1, $2, $3, $4, $5);
#       }
#     }

    my $last_line = `tail -1 tvd_vs_gen`;
    print "$sig1 $last_line" ;
    my @cols = split(" ", $last_line);
    #my ($ndim_tvd, $orth_tvd, $refl_tvd) = @cols[1,2,3];
    push @ndim_tvds, $cols[1];
    push @orthant_tvds, $cols[2];
    push @refl_tvds, $cols[3];
    for (my $i=3; $i<$n_dimensions+3; $i++) {
      push @onedim_tvds, $cols[$i];
    }
    my $ms_dmu = 0.0;
    my @muhats = @cols[4+$n_dimensions..3+2*$n_dimensions];
    for (my $i=0; $i < scalar @muhats; $i++) {
      $ms_dmu += ($muhats[$i] - $mu[$i])**2;
      push @{$dmus{$i}}, $muhats[$i];
    }
    my $rms_dmu = sqrt($ms_dmu); ## ??#

    push @ms_dmus, $ms_dmu;
    my $ksd = pop @cols;
    push @ksds, $ksd;
    my $mean_q = pop @cols;
#print "ms_dmu, ksd, mean q: $ms_dmu  $ksd  $mean_q \n";
    push @mean_qs, $mean_q;
  } # loop over reps
#print "ms dmus: ", join("; ", @ms_dmus), "\n";
#print mean(@ms_dmus), "  ", stddev(@ms_dmus), "\n";
  printf("%8i %10.7fg ", $n_generations, $sig1);
  printf("%10.7g +- %10.7g ", mean(@ndim_tvds), stddev(@ndim_tvds)/sqrt(scalar @ndim_tvds));
  printf("%10.7g +- %10.7g ", mean(@orthant_tvds), stddev(@orthant_tvds)/sqrt(scalar @orthant_tvds));
  printf("%10.7g +- %10.7g ", mean(@refl_tvds), stddev(@refl_tvds)/sqrt(scalar @refl_tvds));
  printf("%10.7g +- %10.7g ", mean(@onedim_tvds), stddev(@onedim_tvds)/sqrt(scalar @onedim_tvds));

  for (0..$n_dimensions-1) {
    printf("%10.7g +- %10.7g ", mean(@{$dmus{$_}}), stddev(@{$dmus{$_}})/sqrt(scalar @{$dmus{$_}}) );
  }
  printf("%10.7g +- %10.7g ", mean(@ms_dmus), stddev(@ms_dmus)/sqrt(scalar @ms_dmus));
  printf("%10.7g +- %10.7g ", mean(@mean_qs), stddev(@mean_qs)/sqrt(scalar @mean_qs));
  printf("%10.7g +- %10.7g ", mean(@ksds), stddev(@ksds)/sqrt(scalar @ksds));
  print "\n\n";
}
