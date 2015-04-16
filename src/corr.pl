#!/usr/bin/perl -w
use strict;
use Math::GSL::Statistics qw( gsl_stats_correlation gsl_stats_lag1_autocorrelation );

my $the_column = shift;

my $max_lag = 100;
my $n_data_points = 400;
my @lines = <>;
# my $the_column = 3;
my @values = ();
my @lagged_values = ();

my $count_points = 0;
for my $the_line (@lines){
  next if($the_line =~ /^\s*#/); # skip comment lines
  my @cols = split(" ", $the_line);
  push @values, $cols[$the_column];
  push @lagged_values, $cols[$the_column];
  $count_points++;
  # last if($count_points >= 100);
}
$n_data_points = scalar @values;
# print join(",", @values), "\n";
print "number of data points read: ", scalar @values, "\n";
 my $lag1_autocorr = gsl_stats_lag1_autocorrelation(\@values, 1, $n_data_points);
print "lag1 autocorr: $lag1_autocorr \n";

# for my $lag (0..100){
#   my $lag1_autocorr = gsl_stats_lag1_autocorrelation(\@values, 1, $n_data_points);
#   my $corr = gsl_stats_correlation(\@values, 1, \@lagged_values, $n_data_points);
#   print "lag: $lag corr: $corr \n";
#   shift @lagged_values;
# }
