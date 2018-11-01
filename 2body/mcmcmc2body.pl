#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $n_updates = 10000;
my $n_dimensions = 2;

my $n_peaks = 2;
my $peak_separation = 2.0;
my $peak_width = 0.25;
my $peak_height_ratio = 0.25;
my $peak_shape_param = 1.0;

my $levels = 2;
my $Thot = 2.0;
my $kernel_scale = 0.5*$peak_separation;


my $proposal_width_factor = 0.5;
my $rng_seed = 1234567;

my $n_thin = 100;

my $output_order = 'T'; # 1:walker, 0:

my $interpolation_power = 0;
my $n_per_level = 2;
my $symmetry = 1;

GetOptions(
           'updates=i' => \$n_updates,
           'dimensions=i' => \$n_dimensions,

           'npeaks=i' => \$n_peaks,
           'peak_separation=f' => \$peak_separation,
           'peak_width=f' => \$peak_width,
           'peak_height_ratio=f' => \$peak_height_ratio,
           'peak_shape=f' => \$peak_shape_param,

           'levels=i' => \$levels,
           'Thot=f' => \$Thot,
           'kernel_scale=f' => \$kernel_scale,

           'proposal_width=f' => \$proposal_width_factor,
           'seed=i' => \$rng_seed,

           'thin=i' => \$n_thin,
           'output_order=s' => \$output_order,
           'interpolation_power=s' => \$interpolation_power, # 'geom' is synonym for 0, 'lin' is synonym for 1;

           'n_per_level=i' => \$n_per_level,
           'symmetry=i' => \$symmetry,
);


my $command = "~/mcmc_toy/2body/mcmc2body  $n_updates $n_dimensions $levels $Thot  ";
$command .= "$n_peaks $peak_separation $peak_width $peak_height_ratio $peak_shape_param  ";
$command .= "$proposal_width_factor $rng_seed $n_thin  ";
my $order = "T";
if($output_order =~ /walker/){
$order = "walker";
}elsif($output_order =~ /cold/){
$order = "cold";
}
$command .= "$order ";
$command .= "$kernel_scale ";

if($interpolation_power =~ /geom/){
   $interpolation_power = 0;
}elsif($interpolation_power =~ /lin/){
   $interpolation_power = 1;
}
$command .= "$interpolation_power  ";

$command .= "$n_per_level ";
$command .= ($symmetry == 1)? "1 " : "0 ";

print "#  $command \n";
system "$command";
