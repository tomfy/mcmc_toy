#!/usr/bin/perl -w
use strict;
use lib '/usr/local/share/perl/5.14.2';
use lib '/home/tomfy/myGnuplotIF';
use Graphics::GnuplotIF qw(GnuplotIF);
use Capture::Tiny qw(:all);
use List::Util qw( min max sum );
#my $enhanced      = 1;
my $persist       = 0;
#my $dmax = 8;

# read in a file (avg_nreps...) 
# and for each block (each T_hot), fir
# dmusqavg to 1/(1/a + x/b)

my $paths_string = shift;
my $xaxis = shift || 'Thot';
my $factor = shift || 0.1;

print "$paths_string   $xaxis   $factor \n";
my @paths = split(/\s*,\s*|\s+/, $paths_string); # comma or space separated
my %path_opt = ();

my @arglist = ();
for my $path (@paths) {
 #  print "path $path \n";
   #   $filename .= '/avg_nreps_tvdetc_vs_gen';
   my ($Thots, $sigma_prop_hots, $effs, $eff_errs, $rchisqrs) = get_effs($path . '/avg_nreps_tvdetc_vs_gen');
   my ($Thot_opt, $max_eff) = get_max_eff($Thots, $effs);
   print "path, opt Thot, eff: $path $Thot_opt, $max_eff \n";
   $path_opt{"$path $Thot_opt"} = $max_eff;
   if ($xaxis eq 'Thot') {
      push @arglist, ($Thots, $effs, $eff_errs, $path);
   } else {
      push @arglist, ($sigma_prop_hots, $effs, $eff_errs, $path);
   }

}

my @sorted_paths = sort { $path_opt{$b} <=> $path_opt{$a} } keys %path_opt;
for(@sorted_paths){
   print "$_  ", $path_opt{$_}, "\n";
}

my $plot1 = Graphics::GnuplotIF->new( persist => $persist, style => 'points ps 0.5 pt 7');
my ($ymin, $ymax) = (5e-7, '*'); #(5e-7, 1e-3);
$plot1->gnuplot_cmd("set log");
$plot1->gnuplot_cmd("set xrange [*:*]"); # [1:40000]");
$plot1->gnuplot_cmd("set yrange [$ymin:$ymax]");
$plot1->gnuplot_cmd("set key left");
# my @cs = map {1e-1*$_} @bs;
# my @c_errs = map {1e-1*$_} @b_errs;
# my @xs = map {0.1*$_} @Thots;
#$plot1->gnuplot_plot_xye(\@Thots, \@bs, \@b_errs);
$plot1->gnuplot_plot_many_xyet(@arglist); # $Thots, $effs, $eff_errs); #, \@xs, \@cs, \@c_errs);
#$plot1->gnuplot_plot_xye(\@Thots, [\@bs, \@cs], [\@b_errs, \@c_errs] );
#my @cs = map {1e-6*$_} @bs;
#$plot1->gnuplot_plot_xy(\@Thots, \@as, \@cs); # {\@b_errs);
$plot1->gnuplot_pause();



sub get_effs{
   my $filename = shift;
   open my $fh, "<", "$filename"; # or "die couldnt open $filename for reading.\n";
   $/ = undef;
   my $data_string = <$fh>;     # read whole file.
   $/ = "\n";

   my @Thots = ();
   my @sigma_prop_hots = ();
   my @effs = ();
   my @eff_errs = ();
   my @rchisqrs = ();

   my @blocks = split(/\n\s*\n/, $data_string);

   my ($std_out, $std_err) = ('AAA', 'BBB');
   for (@blocks) {
      my $reps_summary = read_avgdmusq_vs_npieval($_);
      my @npievs = sort {$a <=> $b} keys %{$reps_summary->{npieval_dmusqavg}};
      my @dmusqavgs = ();
      my @dmusqavg_errs = ();
      while (my ($i, $npiev) = each @npievs) {
         $dmusqavgs[$i] = $reps_summary->{npieval_dmusqavg}->{$npiev};
         $dmusqavg_errs[$i] = $factor*$dmusqavgs[$i];
      }
      my $aguess = $dmusqavgs[0];
      my $bguess = $npievs[-1]*$dmusqavgs[-1];
      my ($a, $a_err, $b, $b_err, $red_chisqr) = 
        fit_ab(\@npievs, \@dmusqavgs, \@dmusqavg_errs, $aguess, $bguess);
      #  print $reps_summary->{Ts}->[-1], "  $a  $a_err     $b  $b_err     $red_chisqr \n";
      my $eff = $reps_summary->{n_dim}*$reps_summary->{variance}/$b;
      my $eff_err = $reps_summary->{n_dim}*$reps_summary->{variance}/$b**2*$b_err;
      print $reps_summary->{Ts}->[-1], " ", $reps_summary->{sigma_props}->[-1],
        " $eff $eff_err  $a $a_err  $b $b_err $red_chisqr \n";
      push @Thots, $reps_summary->{Ts}->[-1];
      push @sigma_prop_hots,  $reps_summary->{sigma_props}->[-1];
      push @effs, $reps_summary->{n_dim}*$reps_summary->{variance}/$b;
      push @eff_errs, $reps_summary->{n_dim}*$reps_summary->{variance}/$b**2*$b_err;
      push @rchisqrs, $red_chisqr;
   }
   print "\n";
   return (\@Thots, \@sigma_prop_hots, \@effs, \@eff_errs, \@rchisqrs);
}                               # end sub get_effs

sub read_avgdmusq_vs_npieval{
   my $block = shift;
   my @lines = split("\n", $block);
   my $obj = {n_dim => 'undef', variance => 'undef', Ts => [], npieval_dmusqavg => {}, sigma_props => []};
   for my $line (@lines) {
      #    print $line, "\n";
      if ($line =~ /^#\s*n_dimensions:\s*(\d+)/) {
         $obj->{n_dim} = $1;
      } elsif ($line =~ /^#.*variance:\s+\S+\s+(\S+)/) {
         $obj->{variance} = $1;
      } elsif ($line =~ /^#.*T, rate, Proposal.*:\s+(\S+)/) {
         push @{$obj->{Ts}}, $1;
         $line =~ /gaussian\s+\S+\s+(\S+)/;
         push @{$obj->{sigma_props}}, $1;
      } elsif ($line =~ /^\s*(\d+)\s+(\S+)/) {
         #   print "$1 $2 \n";
         $obj->{npieval_dmusqavg}->{$1} = $2;
      }
   }
   return $obj;
}

# fit to f(a,b,x) = 1.0/(1/a + x/b)
#return a, e_a, b, e_b, reduced_chisqr
sub fit_ab{
   my $xs = shift;              # array ref
   my $ys = shift;
   my $es = shift; # std errs (so weights will be 1/e**2), or if not an array ref, then
   # a factor f, with  e_i = y_i * f
   my $aguess = shift;
   my $bguess = shift;
   my $fit_limit = shift || 1.0e-10;

   my $errstr = (ref $es eq 'ARRAY')? '3' : '($2*$es)';
   my $tmp_file = $$ . "_xye.data";
   open my $fhtmp, ">", $tmp_file;
   while (my($i, $x) = each @$xs) {
      last if($i >= scalar @$ys  or (ref $es eq 'ARRAY') and $i >= scalar @$es);
      my $err = (ref $es eq 'ARRAY')? $es->[$i] : $ys->[$i]*$es;
      print $fhtmp "$x ", $ys->[$i], " ", $err, "\n";
   }
   close $fhtmp;

   my ($std_out, $std_err) = capture {
      open my $fhgp, '| gnuplot';
      print $fhgp "a=$aguess \n";
      print $fhgp "b=$bguess \n";
      print $fhgp "f(a,b,x) = 1.0/(1.0/a + x/b) \n";
      print $fhgp "FIT_LIMIT = $fit_limit \n";
      print $fhgp "fit f(a,b,x) '$tmp_file'  using 1:2:3 via a,b \n";
   };
   unlink $tmp_file;
   # extract fit parameters, errors, and chisqr from output:
   my @out_lines = split("\n", $std_err); # the fit output goes to stderr, captured to $std_err
   my ($a, $a_err, $b, $b_err, $red_chisqr);
   my $ready = 0;
   for (@out_lines) {
      if (/reduced chisquare.*WSSR\/ndf\s+:\s+(\S+)/) {
         $red_chisqr = $1;
         $ready = 1;
      }
      if ($ready) {
         if (/^\s*a\s+=\s*(\S+)\s+[+]\/[-]\s*(\S+)/) {
            $a = $1; $a_err = $2;
         } elsif (/^\s*b\s+=\s*(\S+)\s+[+]\/[-]\s*(\S+)/) {
            $b = $1; $b_err = $2;
            last;
         }
      }
   }
   return ($a, $a_err, $b, $b_err, $red_chisqr);
}

sub get_max_eff{
my $Thots = shift;
my $effs = shift;

my $max_eff = -1;
my $Thot_opt;
while(my($i, $eff) = each @$effs){
   my $thot = $Thots->[$i];
   if($eff > $max_eff){
      $max_eff = $eff;
      $Thot_opt = $thot;
   }
}
return ($Thot_opt, $max_eff);
}
