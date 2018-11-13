#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Statistics::Basic qw(:all);

my $n_updates = 1000000000; # big, so max_pi_evals will determine stopping point.
my $n_dimensions = 2;

my $n_peaks = 2;
my $peak_separation = 2.0;
my $peak_widths = '0.25,0.25';
my $peak_heights = '1.0,0.25';
my $peak_shape_param = 1.0;

my $levels = 2;
my $Thot_range = '1.2,2.0';
my $kernel_scale = 0.5*$peak_separation;


my $proposal_width_factor = 0.5;
my $rng_seed = 1234567;

my $n_thin = 100;

my $output_order = 'T';         # 1:walker, 0:

my $interpolation_power = 1;
my $n_per_level = 2;
my $symmetry = 1;
my $burn_in = -0.1;
my $max_pi_evals = 100000;
my $n_runs = 10;

my $verbose = 0;
my $Thot_multiplier = 2**0.25;

GetOptions(
           'updates=i' => \$n_updates,
           'dimensions=i' => \$n_dimensions,

           'npeaks=i' => \$n_peaks,
           'peak_separation=f' => \$peak_separation,
           'peak_widths=s' => \$peak_widths,
           'peak_heights=s' => \$peak_heights,
           'peak_shape=f' => \$peak_shape_param,

           # chain architecture params
           'levels=i' => \$levels,
           'n_per_level=i' => \$n_per_level,
           'symmetry=i' => \$symmetry,
           'Thot_range=s' => \$Thot_range,
           'kernel_scale=f' => \$kernel_scale,
           'interpolation_power=s' => \$interpolation_power, # 'geom' is synonym for 0, 'lin' is synonym for 1;
           'proposal_width_factor=f' => \$proposal_width_factor,
   

           'seed=i' => \$rng_seed,
           'burn_in=f' => \$burn_in,
           'thin=i' => \$n_thin,
           'max_pi_evals=i' => \$max_pi_evals,
           'runs=i' => \$n_runs,

           'output_order=s' => \$output_order,
         

       
           'verbose!' => \$verbose,
          );
my $peak_type = '0';
my $x = -0.5*$peak_separation;
my $height = 1.0;
my @widths = split(",", $peak_widths);
my @heights = split(",", $peak_heights);
print "widths: ", join("; ", @widths), "\n";
print "heights: ", join("; ", @heights), "\n";

my $target_string = "$n_dimensions $n_peaks  ";
for (my $i = 0; $i < $n_peaks; $i++) {
   $target_string .= "$peak_type $x," . join(",", ((0) x ($n_dimensions-1))) . "  ";
   $target_string .= $heights[$i] . " " .  $widths[$i] . " $peak_shape_param  ";
   $x += $peak_separation;
 #  $height *= $peak_height_ratio;
}
my @Thots = ();
my @Thot_range = split(',', $Thot_range);
if(scalar @Thot_range == 1){
@Thots = @Thot_range;
}elsif(scalar @Thot_range == 2){
for (my $Thot = $Thot_range[0]; $Thot <= $Thot_range[1] * sqrt($Thot_multiplier); $Thot *= $Thot_multiplier) {
   push @Thots, $Thot;
}}
# print("$target_string \n");
my ($Thot_min, $Thot_max) = split(',', $Thot_range);

#for (my $Thot = $Thot_min; $Thot <= $Thot_max * $Thot_multiplier; $Thot *= $Thot_multiplier) {
for my $Thot (@Thots){ 
  my $arch_string = "  $levels $n_per_level $symmetry $Thot $kernel_scale $proposal_width_factor $interpolation_power  ";

   my $command = "~/mcmc_toy/2body/mcmc2body  $target_string   $arch_string  ";
   # $command .= "$n_peaks $peak_separation $peak_width $peak_height_ratio $peak_shape_param  ";
   $command .= "$rng_seed $burn_in $n_updates $max_pi_evals $n_runs  $n_thin ";
   my $order = "T";
   if ($output_order =~ /walker/) {
      $order = "walker";
   } elsif ($output_order =~ /cold/) {
      $order = "cold";
   }
   $command .= "$order ";

   $command .= "$verbose ";
   # $command .= "$kernel_scale ";

   # if($interpolation_power =~ /geom/){
   #    $interpolation_power = 0;
   # }elsif($interpolation_power =~ /lin/){
   #    $interpolation_power = 1;
   # }
   # $command .= "$interpolation_power  ";

   # $command .= "$n_per_level ";
   # $command .= ($symmetry == 1)? "1 " : "0 ";

   # $command .= "$n_runs ";

   print "#  $command \n";
   $command =~ s/,/ /g;
   print "# $command \n";

   my $output_string = `$command`;

   if ($verbose) {
      print $output_string, "\n\n";
   } else {

      my @lines = split("\n", $output_string);
      my @pooled_abs_errors = ();
      my @abs_errors = (); # elements (1 for each column) are array refs holding errors from the different runs.)
      for (1..$levels*2) {
         push @abs_errors, [];
      }
      for my $a_line (@lines) {
         if ($a_line =~ /^\s*#/) {
            #  print "X $a_line \n";
            if ($a_line =~ /^# Pooled Abs Error:\s*(.*)$/) {
               print "$1 \n";
               push @pooled_abs_errors, $1;
            } elsif ($a_line =~ /^# Abs Errors:\s*(.*)$/) {
               #   print "Z $1\n";
               my @abserrs = split(" ", $1);
               while (my ($i, $err) = each @abserrs) {
                  push @{$abs_errors[$i]}, $err;
               }
            } else {
               print "$a_line \n";
            }
         }
      }
      my $mean = mean(@pooled_abs_errors);
      my $stddev = stddev($mean->query_vector);
      printf("Thot: %8.5f   %8.5f +- %7.5f   ", $Thot, $mean, $stddev);
      while (my ($i, $ar) = each @abs_errors) {
         # my @aerrs = @{$abs_errors[$i]}
         my $mean = mean( $abs_errors[$i] );
         my $stddev = stddev($mean->query_vector);
         printf("%3i  %8.5f +- %7.5f  ", $i, $mean, $stddev/sqrt(scalar @{$abs_errors[$i]}));
      }
      printf("\n");
   }
}
