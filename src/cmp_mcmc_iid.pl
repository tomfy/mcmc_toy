#!/usr/bin/perl -w
use strict;
use List::Util qw ( min max sum );

my $mu = shift;

my @lines = <>;

my @mcmc_gs = ();
my @iid_gs = ();
my $mcmc_g_sum = 0.0;
my $iid_g_sum = 0.0;
my $mcmc_gsquared_sum = 0.0;
my $iid_gsquared_sum = 0.0;
my $count = 0;
my $n_dimensions;
for my $a_line (@lines) {
  if ($a_line =~ /^\s*#/) {
    if ($a_line =~ /n_dimensions:\s+(\S+)/) {
      $n_dimensions = $1;
    }
    next;
  }
  my @cols = split(" ", $a_line);
  my $gen = shift @cols;
  my $temperature = shift @cols;
  my $prob = shift @cols;
  my $mcmc_g = g_function(@cols[0..$n_dimensions-1]);
  my $iid_g = g_function(@cols[$n_dimensions..2*$n_dimensions-1]);
  $mcmc_g_sum += $mcmc_g;
  $iid_g_sum += $iid_g;
  $mcmc_gsquared_sum += $mcmc_g*$mcmc_g;
  $iid_gsquared_sum += $iid_g*$iid_g;
 # push @mcmc_gs, $mcmc_g;
 # push @iid_gs, $iid_g;
  $count++;
}

#my @mcmc_gsquareds = map {$_ * $_} @mcmc_gs;
#my @iid_gsquareds = map {$_ * $_} @iid_gs;

my $mcmc_mu_hat = $mcmc_g_sum/$count; # sum(@mcmc_gs)/scalar @mcmc_gs;
my $iid_mu_hat = $iid_g_sum/$count; # sum(@iid_gs)/scalar @iid_gs;
my $use_mu_hat = (defined $mu)? 0 : 1;
$mu = ($use_mu_hat)? $mcmc_mu_hat : $mu;
my $mcmc_variance_hat = $mcmc_gsquared_sum/$count - $mu*$mu;    # sum(@mcmc_gsquareds)/scalar(@mcmc_gsquareds);
$mu = ($use_mu_hat)? $iid_mu_hat : $mu;
my $iid_variance_hat = $iid_gsquared_sum/$count - $mu*$mu;    # sum(@iid_gsquareds)/scalar(@iid_gsquareds);

printf("%12.7g %12.7g %12.7g %12.7g \n", $mcmc_mu_hat, $mcmc_variance_hat, $iid_mu_hat, $iid_variance_hat);

sub g_function{
  return sum(@_);
}
