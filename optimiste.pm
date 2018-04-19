package optimiste;
use strict;
use warnings;

use base qw(Exporter);
our @ISA = qw(Exporter);
our @EXPORT = qw(all_grad wul_grad conv_fort_num);

sub all_grad {
	my ($control, $grad_file) = @_;
	print "WULEN GRAD START PARAMETERS\n";
	print $control, "\n";
	print $grad_file, "\n";
	print "----------------------------\n";
	open my $res, "<", $control or die "Can't open control: ", $control, "\n";
	my @all_control = <$res>;
	close $res;
	my $end = scalar @all_control;
	my $natom;
	foreach my $str (@all_control) {
		 $natom= $1 if $str =~ m/natoms=(\d+)/g
	}
	my $grad_str = $all_control[$end - 2 - $natom*2];
	my @tmp =  split /\s+/, $grad_str;
	my $grad = $tmp[10];
	open my $in_grad, ">", $grad_file or die "Can't open file for grad";
	print "grad: ",  -$grad/$natom, "\n";
	print $in_grad -$grad/$natom, "\n";
	close $in_grad;
	print "calc grad ended\n";
}
sub wul_grad {
	my ($in_coord, $control, $grad_file) = @_;
	print "WULEN GRAD START PARAMETERS\n";
	print $in_coord, "\n";
	print $control, "\n";
	print $grad_file, "\n";
	print "----------------------------\n";
	open my $res, "<", $control or die "Can't open control: ", $control, "\n";
	my @all_control = <$res>;
	close $res;
	my $end = scalar @all_control;
	my $natom;
	foreach my $str (@all_control) {
		 $natom= $1 if $str =~ m/natoms=(\d+)/g
	}
	my @work = @all_control[$end-$natom-1 .. $end-2];
	foreach my $str (@work) {
		print $str;
	}
	open my $in, "<", $in_coord or die "Can't open coord";
	my @coord = ();
	while (my $str = <$in>) {
		next if $str =~ m/^[\$\#]/g;
		push @coord, $str;
	}
	foreach my $str (@coord) {
		print $str;
	}
	close $in;
	open my $grad, ">", $grad_file or die "Can't open file for grad";
	print scalar @work, "\t", scalar @coord, "\n";
	my @grad;
	foreach my $coord_str (@coord) {
		my $grad_str = shift @work;
		$coord_str =~ s/^\s+//g;
		my @coord_arr = split /\s+/, $coord_str;
		my @grad_arr = split /\s+/, $grad_str;
		shift @grad_arr;
		pop @coord_arr;
		foreach my $cur_coord (@coord_arr) {
			my $cur_grad = shift @grad_arr;
			$cur_grad = conv_fort_num($cur_grad);
			$grad[$1-1] = $cur_grad, if $cur_coord =~ m/^_x(\d+)_/g;
			$grad[$1-1] = -$cur_grad, if $cur_coord =~ m/^-_x(\d+)_/g;
		}
	}
	my $grad_str = join ' ', @grad;
	print $grad $grad_str, "\n";
	close $grad;
	print "calc grad ended\n";
}
sub conv_fort_num {
	my $number = shift;
	$number =~ s/[dD]/e/g;
	$number = sprintf("%.6f", $number);
	return $number;
}
1;
