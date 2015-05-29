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
unlink("avg_nreps_tvdetc_vs_gen");
unlink("avg_accept_info");

my $n_generations = 100000;
my $n_burn_in = 0;
my $n_dimensions = 1;
my $peaks_string = '-0.5, 0.0125, 0.5; 0.5, 0.0125, 0.5';
my $n_peaks;
my $temperatures_string = '1.0';
my $rates_string = '';
my $proposals_string = ''; # 'gaussian, 1.0, 1.4, 1.0'; # if multiple proposals, separate with ; e.g. 'gaussian 0.07 1.5 0.9; gaussian ... '
#my %n_bins_hash = (1 => 4096, 2 => 64, 3 => 16, 4 => 8);
my $n_bins = 4;			# $n_bins_hash{$n_dimensions};
my $n_bins_1d = 360;
my $chain_type = "mcmc";	# or "iid"
my $n_reps = 8;
my $n_factors = 12;
my $f_t_hot = 2;
my $max_t_hot = 10000;
#my $t_factor_exponent = undef;    # t_factor is 2**t_factor_exponent
my $n_temperatures = 1;
#my $max_n_temperatures = 5;
my $cool_rate = 1.0; # rate of within-T moves in all chains but hottest, rel. to hottest.
my $neighbor_swap_only = 1;
#my $max_n_temperatures = 8.01;
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
           'n_temperatures=i' => \$n_temperatures,
           #    'max_n_temperatures=i' => \$max_n_temperatures,
           'n_factors=i' => \$n_factors,
           'f_t_hot=f' => \$f_t_hot,
           'max_t_hot=f' => \$max_t_hot,
           #       't_factor_exponent=f' => \$t_factor_exponent, 
           'cool_rate=f' => \$cool_rate,
           'neighbor_swap_only=i' => \$neighbor_swap_only,
	  );

# 'optimal' (for peak jumping) sigma_prop is A*x(d)  where d is n_dimensions, A is peak spacing;
# here keys are d; values are x(d). (x(d) is approx 1/sqrt(d) for large d)
my %ndim_sig_opt_deltafunction = (1 => 1.000, 2 => 0.8121, 3 => 0.7190, 4 => 0.6541, 5 => 0.5995,
				  6 => 5454, 7 => 0.4811, 8 => 0.4075, 9 => 0.3597, 10 => 0.0.3303,
				  11 => 0.3095, 12 => 0.2934, 13 => 0.2802, 14 => 0.2690, 15 => 0.2593,
				  16 => 0.2507, 17 => 0.2429, 18 => 0.2360, 19 => 0.2296, 20 => 0.2237,
				 );

#$t_factor = 2.0**$t_factor_exponent if(defined $t_factor_exponent);
# my $n_temperatures;
# print "# $min_n_temperatures $max_n_temperatures \n";
#for ($n_temperatures = $min_n_temperatures; $n_temperatures <= $max_n_temperatures; $n_temperatures++) {

# my $tf = $t_factor_factor;
# my @t_factors = ($tf);
# for(2..$n_factors){
#    $tf *= $t_factor_factor;
#    push @t_factors, $tf;
# }
#for my $t_factor (@t_factors){ 
my $Thot = 1.0;
while ($Thot < $max_t_hot*$f_t_hot) {
   print "# $n_temperatures \n";
   my @peak_strings = split(";", $peaks_string);
   $n_peaks = scalar @peak_strings;
   for (@peak_strings) {
      s/,/ /g;
   }
   my $target_peak_sigma = 0.1;
   if ($peak_strings[0] =~ /^\s*(\S+)\s+(\S+)/) {
      $target_peak_sigma = $2;
   }
   # print "# proposals string: $proposals_string \n";

   my $prop_sig_1 = $ndim_sig_opt_deltafunction{$n_dimensions}; # 'opt' peak jumping sigma (peak spacing of 1) for T = 1
   my $prop_sig_2 = 2.38*$target_peak_sigma/sqrt($n_dimensions); # 'opt' within peak sigma for T = 1
   my $prop_p1 = 0.9;
   $temperatures_string =~ s/[,;]/ /g;
   $rates_string =~ s/[,;]/ /g;
   #print "n temperatures: $n_temperatures \n";
   my @temperatures = ((1.0) x ($n_temperatures));
   if ($n_temperatures > 1) {
      while (my ($i, $T) = each @temperatures) { # set up geometric sequence of T's 
         $temperatures[$i] = $Thot**($i/($n_temperatures-1));
      }
   }
   print "temperatures: ", join(", ", @temperatures), "\n";
   my @rates = ();
   if ($rates_string) {
      @rates = split(" ", $rates_string);
   } else {
      if ($n_temperatures == 1) {
         @rates = (1.0);
      } else {
         my $crate = $cool_rate; # /($n_temperatures-1); 
         @rates = (($crate) x ($n_temperatures-1));
         push @rates, 1.0;
      }
   }

   my @proposal_strings = ();
   if ($proposals_string) {
      @proposal_strings = split(";", $proposals_string);
   } else {
      if ($n_temperatures == 1) {
         my $propstr = "gaussian  $prop_sig_1  " . $prop_sig_2*sqrt($temperatures[0]) . " " . 1.0/($cool_rate+1.0);
         push @proposal_strings, $propstr;
      } else {
         while ( my($i, $T) = each @temperatures) {
            my $pr1 = ($i==($n_temperatures-1))? 1.0 : 0.0;
            $pr1 = 0.0; # just all 'local' prop., but growing like sqrt(T)
            my $propstr = "gaussian  $prop_sig_1  " . $prop_sig_2*sqrt($T) . " " . $pr1;
            #  print "XXX $propstr \n";
            push @proposal_strings, $propstr;
         }
      }
   }

   # print "temperatures: ", join(", ", @temperatures), "\n";
   # print "rates: ", join(", ", @rates), "\n";
   # print "proposals: ", join(", ", @proposal_strings), "\n";


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

   #srand();
   #for my $n_temperatures ($min_n_temperatures .. $max_n_temperatures) {


   my $T_hot = $temperatures[-1];
   my $p_string = "$T_hot ";
   my $the_arg_string = $proto_arg_string;
   #  $the_arg_string =~ s/MULTI/$param_val/; # substitute in the particular param value.
   print "# arg string: \n#$the_arg_string \n";
   trials($the_arg_string, $n_reps, $p_string);
   #}
   $Thot *= $f_t_hot;
}

sub trials{                  # do multiple runs with same parameters, 
   my $the_arg_string = shift;
   my $n_reps = shift;
   my $param_val = shift;

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
 #  my %dmus = (); # $dmus{$i} is array_ref to array of dmus for component $i
   # for (0..$n_dimensions-1) {
   #   $dmus{$_} = [];		# $mu[$_] = 0;
   # }
  # my %t_dxhists = ();

   my %type_npieval_avgerr = ('dmusq' => {}, 'ksd' => {}, 'efds' => {}, 'orthtvd' => {}, 'tvd1d' => {}, 'tvd1dall' => {});
   my %accept_info = ('ntry2' => {}, 'nacc2' => {}, 'nswaptry' => {}, 'nswapacc' => {});

   ### 
   #  do the replicate runs ...

   for (1..$n_reps) {
      my $seed = int(rand(1000000)); # get seed for mcmc_toy rng
      system "~/mcmc_toy/src/mcmc_toy  $the_arg_string  $seed  'mcmc' $neighbor_swap_only >  mcmc_toy_samples_tmp";

      ###
      #    read in, store, info from  tvd_etc...
      my ( $nT, $Thot, $rate, $n_pi_eval, $n_gens, $n_try, $n_accept, # 0-3
           $ndim_tvd, $orthant_tvd, $refl_tvd, $oned_tvd, $onedall_tvd, # 4-8
           $dmu_sq, $ksd, $efds, $sq_of_mean_q, $mean_p);
      my $tvd_vs_gen_string = `cat tvd_etc_vs_gen`;
      my @tvd_vs_gen_lines = split("\n", $tvd_vs_gen_string);
      my @cols;
      while (my($i, $tvd_etc_line) = each@tvd_vs_gen_lines) {
         next if($tvd_etc_line =~ /^\s*#/); # skip comments
         @cols = split(" ", $tvd_etc_line);
         ($nT, $Thot, $rate, $n_pi_eval, $n_gens, $n_try, $n_accept, # 0-6
          $ndim_tvd, $orthant_tvd, $refl_tvd, $oned_tvd, $onedall_tvd, # 7-11
          $dmu_sq, $ksd, $efds, $sq_of_mean_q, $mean_p) # 12-16
           = @cols[0..16];
         my $key = $n_pi_eval;
         # $gen_avgdmusq{$key} += $dmu_sq;
         # $gen_avgksd{$key} += $ksd;
         # $gen_avgefds{$key} += $efds;
         # $gen_avgorthtvd{$key} += $orthant_tvd;
         # $gen_avg1dim1tvd{$key} += $oned_tvd;
         # $gen_avg1dimalltvd{$key} += $onedall_tvd;
  $type_npieval_avgerr{dmusq}->{$key} += $dmu_sq;
$type_npieval_avgerr{ksd}->{$key} += $ksd;
$type_npieval_avgerr{efds}->{$key} += $efds;
$type_npieval_avgerr{orthtvd}->{$key} += $orthant_tvd;
$type_npieval_avgerr{tvd1d}->{$key} += $oned_tvd;
$type_npieval_avgerr{tvd1dall}->{$key} += $onedall_tvd;
       
      }

      # store the last tvd etc values (i.e. values at end of each run)
      push @n_tries, $n_try;
      push @n_accepts, $n_accept;
      push @n_pi_evals, $n_pi_eval;
      push @ms_dmus, $dmu_sq;
      push @ksds, $ksd;
      push @efdss, $efds;
      push @sq_of_mean_qs, $sq_of_mean_q;
      push @mean_ps, $mean_p;

      push @ndim_tvds, $ndim_tvd;
      push @orthant_tvds, $orthant_tvd;
      push @refl_tvds, $refl_tvd;

      push @onedim1_tvds, $oned_tvd;
      push @onedimall_tvds, $onedall_tvd;

### store acceptance info:    #############
      store_accept_info('accept_info', \%accept_info);

   }  # end of loop over reps

   my $run_params = `cat run_params`;

   ### output avg over reps of tvdetc_vs_gen ...
   #
   open my $fh1, ">>", "avg_nreps_tvdetc_vs_gen"; # _param_" . $param_val;
   printf $fh1 "$run_params";
   my @sgens = sort {$a <=> $b} keys %gen_avgdmusq;

   for (@sgens) {
      printf $fh1 "%8i %10.7g %10.7g %10.7g %10.7g %10.7g %10.7g \n", $_,
#  $gen_avgdmusq{$_}/$n_reps, $gen_avgksd{$_}/$n_reps, $gen_avgefds{$_}/$n_reps, # cols 2-4
 #         $gen_avgorthtvd{$_}/$n_reps, $gen_avg1dim1tvd{$_}/$n_reps,  $gen_avg1dimalltvd{$_}/$n_reps;
        $type_npieval_avgerr{dmusq}->{$_}/$n_reps, $type_npieval_avgerr{ksd}->{$_}/$n_reps, $type_npieval_avgerr{efds}->{$_}/$n_reps, 
# $gen_avgksd{$_}/$n_reps, $gen_avgefds{$_}/$n_reps, # cols 2-4
#          $gen_avgorthtvd{$_}/$n_reps, $gen_avg1dim1tvd{$_}/$n_reps,  $gen_avg1dimalltvd{$_}/$n_reps;
$type_npieval_avgerr{orthtvd}->{$_}/$n_reps, $type_npieval_avgerr{tvd1d}->{$_}/$n_reps, $type_npieval_avgerr{tvd1dall}->{$_}/$n_reps
   }
   printf $fh1 "\n";
   close $fh1;

   ### output avg over reps of accept info vs T
   #
   open my $fh2, ">>", "avg_accept_info"; # _param_" . $param_val;
   printf $fh2 "$run_params";
   my @sTs = sort {$a <=> $b} keys %{$accept_info{ntry2}};

   for my $T (@sTs) {
      my $acc_rate2 = $accept_info{nacc2}->{$T} / $accept_info{ntry2}->{$T};
      printf $fh2 ("%10.7g %8i %8i %10.7g  ", # %8i %8i %10.7g \n",
      $T, $accept_info{ntry2}->{$T}, $accept_info{nacc2}->{$T}, $acc_rate2);
      if (exists $accept_info{nswaptry}->{$T} and $accept_info{nswaptry}->{$T} > 0) {
         my $swap_acc_rate = $accept_info{nswapacc}->{$T} / $accept_info{nswaptry}->{$T};
         printf $fh2 (" %8i %8i %10.7g \n", $accept_info{nswaptry}->{$T}, $accept_info{nswapacc}->{$T}, $swap_acc_rate );
      }else{
         printf $fh2 "\n";
      }
   }
   printf $fh2 "\n";
   close $fh2;



   #  my $run_params = `cat run_params`;
   printf("$run_params");
   printf("# averages based on $n_reps reps. \n");
   printf("%s   ", $param_val);
   printf("%8i  ", $n_generations);
   printf("%10.6g +- %8.4g  ", sum(@n_accepts)/sum(@n_tries), stddev(@n_accepts)/mean(@n_tries)/sqrt(scalar @n_accepts));
   printf("%10.6g +- %8.4g  ", mean(@n_pi_evals), stddev(@n_pi_evals)/sqrt(scalar @n_pi_evals));

   printf("%10.6g +- %8.4g  ", mean(@ndim_tvds), stddev(@ndim_tvds)/sqrt(scalar @ndim_tvds));
   printf("%10.6g +- %8.4g  ", mean(@orthant_tvds), stddev(@orthant_tvds)/sqrt(scalar @orthant_tvds));
   printf("%10.6g +- %8.4g  ", mean(@refl_tvds), stddev(@refl_tvds)/sqrt(scalar @refl_tvds));
   printf("%10.6g +- %8.4g  ", mean(@onedim1_tvds), stddev(@onedim1_tvds)/sqrt(scalar @onedim1_tvds));
   printf("%10.6g +- %8.4g  ", mean(@onedimall_tvds), stddev(@onedimall_tvds)/sqrt(scalar @onedimall_tvds));
   printf("   ");
   printf("%10.6g +- %8.4g  ", mean(@ms_dmus), stddev(@ms_dmus)/sqrt(scalar @ms_dmus));
   printf("%10.6g +- %8.4g  ", mean(@ksds), stddev(@ksds)/sqrt(scalar @ksds));
   printf("%10.6g +- %8.4g  ", mean(@efdss), stddev(@efdss)/sqrt(scalar @efdss));
   printf("%10.6g +- %8.4g  ", mean(@sq_of_mean_qs), stddev(@sq_of_mean_qs)/sqrt(scalar @sq_of_mean_qs));
   printf("%10.6g +- %8.4g  ", mean(@mean_ps), stddev(@mean_ps)/sqrt(scalar @mean_ps));
   print "\n";
}

sub store_accept_info{
   my $accept_info_filename = shift;
   my $accept_info = shift;     # hashref 
   my $accept_info_string = `cat $accept_info_filename`;
   my @accept_info_lines = split("\n", $accept_info_string);
   for my $a_line (@accept_info_lines) {
      my @cols = split(" ", $a_line);
      my $T = $cols[2];
      $accept_info->{'ntry2'}->{$T} += $cols[5];
      $accept_info->{'nacc2'}->{$T} += $cols[6];
      if (scalar @cols >= 9) {
         $accept_info->{'nswaptry'}->{$T} += $cols[7];
         $accept_info->{'nswapacc'}->{$T} += $cols[8];
      }
   }
}

# sub dx_histograms{
#    my $mt_output_string = shift;
#    my $temperature_dxhists = shift; # keys: temperatures, values: arrayref of jumpdistancesquared
#    my $inv_bin_width = 16;
#    my $n_bins = int(3*$inv_bin_width);
#    my @lines = split("\n", $mt_output_string);
#    my $first_line = shift @lines;
#    while ($first_line =~ /^\s*#/) {
#       $first_line = shift @lines;
#    }
#    my @first_line_cols = split(" ", $first_line);
#    my @prev_xs = @first_line_cols[3..scalar @first_line_cols - 1];

#    for my $a_line (@lines) {
#       next if($a_line =~ /^\s*#/);
#       #   if($a_line =~ /1.0000/){
#       my @cols = split(" ", $a_line);
#       my $gen = shift @cols;
#       my $temperature = shift @cols;
#       #   next if(abs($temperature - 1.0) > 1e-10);
#       #print "T: $temperature \n";
#       my $prob = shift @cols;
#       my @xs = @cols;
#       my $dxsq = 0.0;
#       #   print join(", ", @prev_xs), " :: ", join("; ", @xs), "\n";
#       for (my $i=0; $i<scalar @xs; $i++) {
#          $dxsq += ($xs[$i] - $prev_xs[$i])**2;
#       }
#       #  print "dxsq: $dxsq \n";
   
#       if (! exists $temperature_dxhists->{$temperature}) {
#          # print "TT: $temperature \n";
#          $temperature_dxhists->{$temperature} = [(0) x $n_bins];
#       }
#       my $bin = int(sqrt($dxsq)*$inv_bin_width);
#       #   print "bin: $bin \n";
#       if ($bin >= $n_bins) {
#          $bin = $n_bins-1;
#       }
#       $temperature_dxhists->{$temperature}->[$bin]++;
#       @prev_xs = @xs;
#    }
# }
