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

my $n_generations = 100000;
my $n_burn_in = 0;
#my $n_replicates = 1; # > 1 not implemented 
my $n_dimensions = 1;
my $peaks_string = '-0.5, 0.1, 0.5; 0.5, 0.1, 0.5';
my $n_peaks;
my $temperatures_string = '1.0';
my $rates_string = '1.0';
my $proposals_string = 'gaussian, 1.0, 1.4, 1.0'; # if multiple proposals, separate with ; e.g. 'gaussian 0.07 1.5 0.9; gaussian ... '
#my %n_bins_hash = (1 => 4096, 2 => 64, 3 => 16, 4 => 8);
my $n_bins = 4;			# $n_bins_hash{$n_dimensions};
my $n_bins_1d = 360;
my $chain_type = "mcmc";	# or "iid"
my $n_reps = 8;
GetOptions(
	   'n_dimensions=i' => \$n_dimensions,
	   'peaks=s' => \$peaks_string, # e.g. '-0.5,0.05,0.5;0.5,0.05,0.5'  -> 2 peaks
	   'temperatures=s' => \$temperatures_string, # e.g. '1.0, 1.1, 1.2, 2.0';
	   'rates=s' => \$rates_string,
	   'proposals=s' => \$proposals_string, # e.g. 'ball,0.1,1.5,0.9'
	   'n_bins=i' => \$n_bins,
	   'n_1d_bins=i' => \$n_bins_1d,
	   'n_generations=i' => \$n_generations,
	   'n_burnin=i' => \$n_burn_in,
	   'type=s' => \$chain_type,
	   'n_reps=i' => \$n_reps,
	  );

my %ndim_sig_opt_deltafunction = (1 => 1.000, 2 => 0.8121, 3 => 0.7190, 4 => 0.6541, 5 => 0.5995,
				  6 => 5454, 7 => 0.4811, 8 => 0.4075, 9 => 0.3597, 10 => 0.0.3303,
				  11 => 0.3095, 12 => 0.2934, 13 => 0.2802, 14 => 0.2690, 15 => 0.2593,
				  16 => 0.2507, 17 => 0.2429, 18 => 0.2360, 19 => 0.2296, 20 => 0.2237,
				 );
my @peak_strings = split(";", $peaks_string);
$n_peaks = scalar @peak_strings;
for (@peak_strings) {
  s/,/ /g;
}
my $target_peak_sigma = 0.1;
if ($peak_strings[0] =~ /^\s*(\S+)\s+(\S+)/) {
  $target_peak_sigma = $2;
}
 print "# proposals string: $proposals_string \n";

my $prop_sig_1 = $ndim_sig_opt_deltafunction{$n_dimensions}; # 'opt' peak jumping sigma (peak spacing of 1)
my $prop_sig_2 = 2.38*$target_peak_sigma/sqrt($n_dimensions); # 'opt' withing peak sigma
my $prop_p1 = 0.9;
$temperatures_string =~ s/[,;]/ /g;
$rates_string =~ s/[,;]/ /g;
my @temperatures = split(" ", $temperatures_string);
my $n_temperatures = scalar @temperatures;
my @rates = split(" ", $rates_string);
my @proposal_strings = split(";", $proposals_string);
my @param_values = ();
for (@proposal_strings) {
  if(s/MULTI\(([^)]*)\)/MULTI/){# e.g. -proposals 'gaussian MULTI(0.5,0.7,1.0) 0.15 0.7; gaussian 1.0 0.15 0.9
    @param_values = split(/[, ]+/, $1);
    print "param values: ", join(';', @param_values), "\n";
  }
  s/[,]/ /g;
}
while (scalar @proposal_strings < scalar @temperatures) {
  push @proposal_strings, $proposal_strings[scalar @proposal_strings -1]; # duplicate last proposal to get enough
}

my $proto_arg_string = "$n_dimensions  $n_peaks  ";
$proto_arg_string .= join("  ", @peak_strings) . "  ";
$proto_arg_string .= " $n_generations  ";
$proto_arg_string .= "$n_temperatures ";
for my $i (0..$n_temperatures-1) {
  my $T = $temperatures[$i];
  my $rate = $rates[$i];
  $proto_arg_string .= "$T $rate " . $proposal_strings[$i] . "  ";
}
$proto_arg_string .= "$n_bins  $n_bins_1d  ";
# print "$n_generations ;;; $proto_arg_string \n"; exit;
# my $seed = srand();
# $proto_arg_string .= "$seed ";

# print "# $proto_arg_string \n";

# my $arg_string = "3   2   -0.5 0.06 0.5   0.5 0.06 0.5   "; # n_dim, n_modes, targ mean, sigma, weight, ...
# $arg_string .= "0 2000000 1   "; # burn-in steps, post-burn-in steps, replicates
# #$arg_string .= "2   1.0 ball 0.06 1.0 0.9   2.0 ball 0.06 1.0 0.5   "; # n_temperatures, ...
# $arg_string .= "1   1.0  ball 0.1 1.5 0.7   "; # n_temperatures, then T_i, proposal_i
# $arg_string .= "25";  # n_bins

#my $tvd_limit = 0.05 ;
my $dmu_sq_limit_f = 0.001;
my $dmu_sq_limit;
# my $n_reps = 8;
srand();

my $sig1 = 1.0;
my $the_arg_string;
# while (0) {
#   $the_arg_string = $proto_arg_string;
#   $the_arg_string =~ s/NGENERATIONS/$n_generations/;
#   $the_arg_string =~ s/SIGMA1/$sig1/;
#   $the_arg_string =~ s/SIGMA2/$prop_sig_2/;
#   $the_arg_string =~ s/PROP_P1/$prop_p1/;
#   my $seed = int(rand(1000000));
#   my $mcmc_toy_out =  `~/mcmc_toy/src/mcmc_toy $the_arg_string $seed mcmc`;
#   #  print "$mcmc_toy_out \n";
#   my $last_line = `tail -1 tvd_vs_gen`;
#   #  exit;
#   my @cols = split(" ", $last_line);
#   my $ngens = shift @cols;
#   my $p_accept = shift @cols;
#   my $dmusq = shift @cols;

#   my $run_params = `cat run_params`;
#   my @runparam_lines = split("\n", $run_params);
#   my ($T, $shape, $s1, $s2, $p1);
#   my ($mean, $variance);
 
#   for (@runparam_lines) {
#     if (/T, Proposal\S+:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
#       ($T, $shape, $s1, $s2, $p1) = ($1, $2, $3, $4, $5);
#     } elsif (/target\s+mean,\s*variance:\s+(\S+)\s+(\S+)/) {
#       ($mean, $variance) = ($1, $2);
#     }
#   }
#   $dmu_sq_limit = $dmu_sq_limit_f * $variance;
#   # print "mean, variance: $mean $variance dmusq: $dmusq\n";
#   #exit;

#   # my @mcmc_1d_tvds = @tvds[3 .. 3+$n_dimensions-1];
#   # my $avg_1d_tvd = mean(@mcmc_1d_tvds);
#   last if($dmusq < $dmu_sq_limit);
#   $n_generations = int( 1.2*$n_generations*($dmusq/$dmu_sq_limit) );
# }
#exit;
my @sig_1s = (0.14, 0.20, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 4.0, 7.0);
# @sig_1s = (0.2, 0.45, 0.7, 1.0, 1.4);
# for ($sig1 = ; $sig1 < 5.0; $sig1 *= 1.3)
@sig_1s = (0.15, 0.45, 0.75, 1.5, 5.0);
@sig_1s = (0.15, 0.25, 0.5, 0.8, 1.0, 1.25, 1.6, 2.3, 3.8, 7.0);
@sig_1s = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4);
@sig_1s = (0.55, 0.7, 0.9, 1.1, 1.3, 1.6, 2.0, 5.0);
@sig_1s = (0.4, 1.0, 3.0);

#@sig_1s = (0.128, 0.164, 0.2, 0.256, 0.328, 0.4, 0.512); #, 0.676, 0.8, 1.0, 1.28, 1.64, 2.0, 3.28, 4.0, 5.12, 6.76, 8.0);
#my $sig_1 = $prop_sig_1;
# my @params = (0.3, 0.5, 0.9);
for my $param_val (@param_values) {
  trials($n_reps, $proto_arg_string, $param_val);
}

sub trials{
  my $n_reps = shift;
  my $proto_arg_string = shift;
  my $param_val = shift;

  my $the_arg_string = $proto_arg_string;
#  print "$proto_arg_string \n";
  $the_arg_string =~ s/MULTI/$param_val/;

  my @n_tries = ();
  my @n_accepts = ();
  my @n_pi_evals = ();
  my @ms_dmus = ();
  my @ksds = ();
  my @efdss = ();
  my @sq_of_mean_qs = ();
  my @mean_ps = ();

  my @ndim_tvds = ();
  my @orthant_tvds = ();
  my @refl_tvds = ();
  my @onedim1_tvds = ();
  my @onedimall_tvds = ();

  my %gen_avgdmusq = ();
  my %gen_avgksd = ();
  my %gen_avgefds = ();
  my %gen_avgorthtvd = ();
  my %gen_avg1dim1tvd = ();
  my %gen_avg1dimalltvd = ();
#  print "# the arg string: $the_arg_string \n";
  #  my $mu = 0.0; # get from run_params file - may not always be zero!
  my %dmus = (); # $dmus{$i} is array_ref to array of dmus for component $i
  # for (0..$n_dimensions-1) {
  #   $dmus{$_} = [];		# $mu[$_] = 0;
  # }
my %t_dxhists = ();
  for (1..$n_reps) {
    my $seed = int(rand(1000000));
 #     print "~/mcmc_toy/src/mcmc_toy  $the_arg_string  $seed  'mcmc' \n";
    
    system "~/mcmc_toy/src/mcmc_toy  $the_arg_string  $seed  'mcmc' >  mcmc_toy_samples_tmp";
 my $mt_out_string = `cat mcmc_toy_samples_tmp`;
    dx_histograms($mt_out_string, \%t_dxhists);



    #my ($n_gens, $n_accept, $dmu_sq, $ksd, $ads, $mean_q, $mean_p, $ndim_tvd, $orthant_tvd, $refl_tvd, $oned_tvd);
    my ($n_gens, $n_try, $n_accept, $n_pi_eval, # 0-3
	$ndim_tvd, $orthant_tvd, $refl_tvd, $oned_tvd, $onedall_tvd, # 4-8
	$dmu_sq, $ksd, $efds, $sq_of_mean_q, $mean_p);
    my $tvd_vs_gen_string = `cat tvd_etc_vs_gen`;
    my @tvd_vs_gen_lines = split("\n", $tvd_vs_gen_string);
    my @cols;
    for (@tvd_vs_gen_lines) {
      @cols = split(" ", $_);
      ($n_gens, $n_try, $n_accept, $n_pi_eval, # 0-3
       $ndim_tvd, $orthant_tvd, $refl_tvd, $oned_tvd, $onedall_tvd, # 4-8
       $dmu_sq, $ksd, $efds, $sq_of_mean_q, $mean_p) # 9-13
	= @cols[0..13];

print "n_pi_eval: $n_pi_eval \n";
      $gen_avgdmusq{$n_gens} += $dmu_sq;
      $gen_avgksd{$n_gens} += $ksd;
      $gen_avgefds{$n_gens} += $efds;
      $gen_avgorthtvd{$n_gens} += $orthant_tvd;
      $gen_avg1dim1tvd{$n_gens} += $oned_tvd;
      $gen_avg1dimalltvd{$n_gens} += $onedall_tvd;
    }

 #   print "$sig1 ", join("  ", @cols), "\n"; 
 #   die "n generations inconsistency: $n_generations, $n_gens \n" if($n_gens != $n_generations);
    push @n_tries, $n_try;
    push @n_accepts, $n_accept;
    push @n_pi_evals, $n_pi_eval;
    print "n_pi_evals: ", join("; ", @n_pi_evals), "\n";
    push @ms_dmus, $dmu_sq;
    push @ksds, $ksd;
    push @efdss, $efds;
    push @sq_of_mean_qs, $sq_of_mean_q;
    push @mean_ps, $mean_p;

    push @ndim_tvds, $ndim_tvd;
    push @orthant_tvds, $orthant_tvd;
    push @refl_tvds, $refl_tvd;

    push @onedim1_tvds, $oned_tvd;
    push @onedimall_tvds, $oned_tvd;
  }				# loop over reps
  open my $fh1, ">", "avg_nreps_tvdetc_vs_gen_param_" . $param_val;
  my @sgens = sort {$a <=> $b} keys %gen_avgdmusq;

  for (@sgens) {
    printf $fh1 "%8i %10.7g %10.7g %10.7g %10.7g %10.7g %10.7g \n", $_,
      $gen_avgdmusq{$_}/$n_reps, $gen_avgksd{$_}/$n_reps, $gen_avgefds{$_}/$n_reps, # cols 2-4
	$gen_avgorthtvd{$_}/$n_reps, $gen_avg1dim1tvd{$_}/$n_reps,  $gen_avg1dimalltvd{$_}/$n_reps;
  }
  close $fh1;

  #print "ms dmus: ", join("; ", @ms_dmus), "\n";
  #print mean(@ms_dmus), "  ", stddev(@ms_dmus), "\n";
  printf("%10.7g ", $param_val);
  printf("%8i  ", $n_generations);
  printf("%10.7g +- %10.7g ", mean(@n_accepts)/$n_generations, stddev(@n_accepts)/$n_generations/sqrt(scalar @n_accepts));
  printf("%10.7g +- %10.7g ", mean(@n_pi_evals), stddev(@n_pi_evals)/sqrt(scalar @n_pi_evals));

  printf("%10.7g +- %10.7g ", mean(@ndim_tvds), stddev(@ndim_tvds)/sqrt(scalar @ndim_tvds));
  printf("%10.7g +- %10.7g ", mean(@orthant_tvds), stddev(@orthant_tvds)/sqrt(scalar @orthant_tvds));
  printf("%10.7g +- %10.7g ", mean(@refl_tvds), stddev(@refl_tvds)/sqrt(scalar @refl_tvds));
  printf("%10.7g +- %10.7g ", mean(@onedim1_tvds), stddev(@onedim1_tvds)/sqrt(scalar @onedim1_tvds));
  printf("%10.7g +- %10.7g ", mean(@onedimall_tvds), stddev(@onedimall_tvds)/sqrt(scalar @onedimall_tvds));

  printf("%10.7g +- %10.7g ", mean(@ms_dmus), stddev(@ms_dmus)/sqrt(scalar @ms_dmus));
  printf("%10.7g +- %10.7g ", mean(@ksds), stddev(@ksds)/sqrt(scalar @ksds));
  printf("%10.7g +- %10.7g ", mean(@efdss), stddev(@efdss)/sqrt(scalar @efdss));
  printf("%10.7g +- %10.7g ", mean(@sq_of_mean_qs), stddev(@sq_of_mean_qs)/sqrt(scalar @sq_of_mean_qs));
  printf("%10.7g +- %10.7g ", mean(@mean_ps), stddev(@mean_ps)/sqrt(scalar @mean_ps));
  print "\n\n";

  open my $fhdxh, ">", "dxhist_$param_val";
  my @stemperatures = sort {$a <=> $b} keys %t_dxhists;
#  my $t_string = "temperatures: " . join(", ", @stemperatures);
  # print "XXX: ", $stemperatures[0], "\n";
  # print  "# $t_string \n";
  # printf $fhdxh "# $t_string \n";
  my @dxhist0 = @{$t_dxhists{$stemperatures[0]}};
  print  $fhdxh "# temperature: ", $stemperatures[0], "\n";
#  printf  "XXX: ", $stemperatures[0], "\n";
  while (my ($i, $v) = each @dxhist0) {
    print $fhdxh "$i $v\n";
  }
}




sub dx_histograms{
  my $mt_output_string = shift;
  my $temperature_dxhists = shift; # keys: temperatures, values: arrayref of jumpdistancesquared
  my $inv_bin_width = 16;
  my $n_bins = 2*$inv_bin_width;
  my @lines = split("\n", $mt_output_string);
  my $first_line = shift @lines;
  my @first_line_cols = split(" ", $first_line);
  my @prev_xs = @first_line_cols[3..scalar @first_line_cols - 1];

  for my $a_line (@lines) {
    if($a_line =~ /1.0000/){
    my @cols = split(" ", $a_line);
    my $gen = shift @cols;
    my $temperature = shift @cols;
 #   next if(abs($temperature - 1.0) > 1e-10);
#print "T: $temperature \n";
    my $prob = shift @cols;
    my @xs = @cols;
    my $dxsq = 0.0;
  #  print join(", ", @prev_xs), " :: ", join("; ", @xs), "\n";
    for (my $i=0; $i<scalar @xs; $i++) {

      $dxsq += ($xs[$i] - $prev_xs[$i])**2;
    }
  #  print "dxsq: $dxsq \n";
   
    if (! exists $temperature_dxhists->{$temperature}) {
#print "TT: $temperature \n";
      $temperature_dxhists->{$temperature} = [(0) x $n_bins];
    }
    my $bin = int(sqrt($dxsq)*$inv_bin_width);
 #   print "bin: $bin \n";
    if ($bin >= $n_bins) {
      $bin = $n_bins-1;
    }
    $temperature_dxhists->{$temperature}->[$bin]++;
    @prev_xs = @xs;
  }
}
print "temperatures: ", join("; ", keys %$temperature_dxhists), "\n";
# return $temperature_dxhists;
}
