#!/usr/bin/perl -w
use strict;
use List::Util qw 'min max sum';
use List::MoreUtils qw 'pairwise';


my ($n_dims, $n_Ts, $peak_separation, $height_ratio) = (undef, undef, undef, undef);
my $n_data_lines = 0;

my $n_print = 100;
my $n_print_factor = 1.2;

while (<>) {
   if (/^\s*#/) {
      $n_dims = $1 if(/n_dims:\s+(\d+)/);
      $n_Ts = $1 if(/n_Ts:\s+(\d+)/);
      $peak_separation = $1 if(/peak separation:\s+([0-9.]+)/);
      $height_ratio = $1 if(/height ratio:\s+([0-9.]+)/);
   }
   last if(defined $n_dims  and  defined $n_Ts  and  defined $peak_separation  and  defined $height_ratio);
}
my @true_x_mean = ( (0) x $n_dims);
$true_x_mean[0] = 0.5*$peak_separation*($height_ratio - 1.0)/(1.0 + $height_ratio);
my @xsum = ( (0) x $n_dims);
my @ysum = ( (0) x $n_dims);
my $i_update;
while (<>) {
   next if(/^\s*#/);
   my @cols = split(" ", $_);
   $i_update = shift @cols;
   my $itx = 0;
   my $ity = 2*$n_Ts - 1;
   $n_data_lines++;

   my $offset = 2;
   my @x =  @cols[$offset+1..$offset+$n_dims];
   accumulate_vector_sum(\@xsum, \@x);

   $offset = 2 + (2*$n_Ts - 1)*($n_dims+3);
   my @y =  @cols[$offset+1..$offset+$n_dims];
   accumulate_vector_sum(\@ysum, \@y);

      if ($i_update >= $n_print) {
         my @estimated_x_mean = map($_/$n_data_lines, @xsum);
         my @estimated_y_mean = map($_/$n_data_lines, @ysum);
         printf("%8i  %10.8g   %10.8g \n", $i_update, 
                abs_error($n_dims, \@true_x_mean, \@estimated_x_mean),
                abs_error($n_dims, \@true_x_mean, \@estimated_y_mean) );
         $n_print *= $n_print_factor;
      }
}

my @estimated_x_mean = map($_/$n_data_lines, @xsum);
my @estimated_y_mean = map($_/$n_data_lines, @ysum);
printf("%8i  %10.8g   ", $i_update, abs_error($n_dims, \@true_x_mean, \@estimated_x_mean));
printf("%10.8g \n", abs_error($n_dims, \@true_x_mean, \@estimated_y_mean));
# my @xavg = map($_/$n_data_lines, @xsum);
# my $outstring = sprintf("$n_data_lines  est mean x: (");
# for(@xavg){
#    $outstring .= sprintf("%5.3f, ", $_);
# } $outstring =~ s/,\s*$/)  /;
# my @dfm = pairwise {$a - $b} @xavg, @true_x_mean;


# my $sqrd_error = sum( map($_*$_, @dfm));

# $outstring .= sprintf("est-true mean x: (");
# for(@dfm){
#    $outstring .= sprintf("%6.4f, ", $_);
# } $outstring =~ s/,\s*$/)  /;
# $outstring .= "abs error: ". sprintf("%10.8g \n", sqrt($sqrd_error));
# printf $outstring;

# print abs_error($n_dims, \@true_x_mean, \@xavg), "\n";


sub abs_error{
   my $n_dims = shift;
   my $true_mean = shift;
   my $estimated_mean = shift;
   my @ds = ();
   for(my $i = 0; $i < $n_dims; $i++){
      push @ds, $estimated_mean->[$i] - $true_mean->[$i];
   }
   return sqrt(sum(map($_*$_, @ds)));
}

sub accumulate_vector_sum{
   my $vsum = shift;
   my $vec = shift;
  while (my ($i, $c) = each @{$vec}) {
         $vsum->[$i] += $c;
      }
return $vsum;
}
