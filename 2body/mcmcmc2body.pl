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

my $n_Ts = 2;
my $Thot = 2.0;
my $kernel_scale = 0.5*$peak_separation;


my $proposal_width = 0.5;
my $rng_seed = 1234567;

my $n_thin = 100;

my $output_order = 'T'; # 1:walker, 0: 

GetOptions(
           'updates=i' => \$n_updates,
           'dimensions=i' => \$n_dimensions,

           'npeaks=i' => \$n_peaks,
           'peak_separation=f' => \$peak_separation,
           'peak_width=f' => \$peak_width,
           'peak_height_ratio=f' => \$peak_height_ratio,
           'peak_shape=f' => \$peak_shape_param,

           'nTs=i' => \$n_Ts,
           'Thot=f' => \$Thot,
           'kernel_scale=f' => \$kernel_scale,

           'proposal_width=f' => \$proposal_width,
           'seed=i' => \$rng_seed,

           'thin=i' => \$n_thin,
           'output_order=s' => \$output_order,
);


my $command = "~/mcmc_toy/2body/mcmc2body  $n_updates $n_dimensions $n_Ts $Thot  ";
$command .= "$n_peaks $peak_separation $peak_width $peak_height_ratio $peak_shape_param  ";
$command .= "$proposal_width $rng_seed $n_thin  ";
my $order = "T";
if($output_order =~ /walker/){
$order = "walker";
}elsif($output_order =~ /cold/){
$order = "cold only";
}
$command .= "$order ";
$command .= "$kernel_scale ";

print "#  $command \n";
system "$command";
