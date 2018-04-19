package dscf;
use strict;
use warnings;

use base qw(Exporter);
our @ISA = qw(Exporter);
our @EXPORT = qw(chdft_func add_cube chscfconv so_params no_rel_params);

sub chdft_func {
	my ($file, $functional) = @_;
	open my $in, "<", $file or die "Can't open $file\n";
	open my $out, ">", 'control.tmp' or die "Can't open control.tmp\n";
	while( my $str = <$in> ) {
		$str = "\$dft-functional $functional\n" if $str =~ m/^\$dft-functional/;
		print $out $str;
	}
	close $in;
	close $out;
	system ('mv control.tmp '.$file);
}
sub add_cube {
	my ($points, $step) = @_;
	my $file = 'control';
	my $init = -$points*$step/2;
	open my $in_coord, "<", 'coord' or die "Can't open coord\n";
	my @coords = <$in_coord>;
	close $in_coord;
	my $atom = $coords[1];
	$atom =~ s/^\s+//g;
	my ($x, $y, $z, undef) = split /\s+/, $atom;
	my $x0 = sprintf("%.2f", $init + $x);
	my $y0 = sprintf("%.2f", $init + $y);
	my $z0 = sprintf("%.2f", $init + $z);
	my $cube_params = <<"EOF";
\$density-cube
    x0= $x0 y0= $y0 z0= $z0
    nx= $points ny= $points nz= $points
    dx= $step dy= $step dz= $step
EOF
	open my $in, "<", $file or die "Can't open $file\n";
	open my $out, ">", 'control.tmp' or die "Can't open control.tmp\n";
	while( my $str = <$in> ) {
		$str.=$cube_params  if $str =~ m/^\$dft-functional/;
		print $out $str;
	}
	close $in;
	close $out;
	system ('mv control.tmp '.$file);
}
sub chscfconv {
	my ($file, $limit) = @_;
	open my $in, "<", $file or die "Can't open $file\n";
	open my $out, ">", 'control.tmp' or die "Can't open control.tmp\n";
	while( my $str = <$in> ) {
		$str = "\$scfconv $limit\n" if $str =~ m/^\$scfconv/;
		print $out $str;
	}
	close $in;
	close $out;
	system ('mv control.tmp '.$file);
}
sub so_params {
	my %so_params = ('test' =>   'ENERGY CONVERGED !', #  line to be found in the listing to ensure 
			  #successful termination of a QC code
			 'line' =>   '! Total energy', #   line containing the answer (can coincide with <<test>>)
			 'word' =>   4, #  position of the answer in the abovementioned line (0 - first, 1-second etc)
			 'maxi' =>   100, # maximal number of iterations
			 'rmin' =>   1.2, # for lr mode, minimum ratio
			 'step' =>   0.1, # numerical differentiation step 
			 'grad' =>   'y', # use wullen gradient
			 # (())   1  12.0 # this option abesent now
			 'step' =>   0.1, #  numerical differentiation step
			 'mode' =>   'xc',
			 # MODE:
			 #     lr - для первого ряда
			 #     xc - для второго и расстояния между атомами
			 #    first symbol: logarithmic sale if <<l>>, linear otherwise
			 #    second symbol: if <<r>> and the first one is <<l>>, minimal ratio of 
			 #                   two neighbor parameter is given by <<rmin>> line
			 #                   if <<c>>, the step for fd gradient calcns is independent 
			 #                   on the parameter value
			 #    lx (or lr) are recommended for exponent optimization,
			 #    xc for weight optimization
			 'without_energy' => 'n',
			 'direct' => 'y',
			);
	return  %so_params;
}
sub no_rel_params {
	my %no_rel_params = ('test' =>   'ENERGY CONVERGED !', #  line to be found in the listing to ensure 
			  #successful termination of a QC code
			 'line' =>   'Total energy', #   line containing the answer (can coincide with <<test>>)
			 'word' =>   4, #  position of the answer in the abovementioned line (0 - first, 1-second etc)
			 'maxi' =>   100, # maximal number of iterations
			 'rmin' =>   1.2, # for lr mode, minimum ratio
			 'step' =>   0.1, # numerical differentiation step 
			 'grad' =>   'y', # use wullen gradient
			 # (())   1  12.0 # this option abesent now
			 'step' =>   0.1, #  numerical differentiation step
			 'mode' =>   'xc',
			 # MODE:
			 #     lr - для первого ряда
			 #     xc - для второго и расстояния между атомами
			 #    first symbol: logarithmic sale if <<l>>, linear otherwise
			 #    second symbol: if <<r>> and the first one is <<l>>, minimal ratio of 
			 #                   two neighbor parameter is given by <<rmin>> line
			 #                   if <<c>>, the step for fd gradient calcns is independent 
			 #                   on the parameter value
			 #    lx (or lr) are recommended for exponent optimization,
			 #    xc for weight optimization
			 'without_energy' => 'n',
			 'direct' => 'y',
			);
	return  %no_rel_params;
}
1;
