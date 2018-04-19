#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use lib '/s/ls2/home/users/yuriyd/bin/perl_lib';
use optimiste;
use dscf;
use utils;

my $DEBUG = 'yes';
my $root_dir = getcwd;
our %params; # FOR PARAMETERS
our @target_files;

help() if !@ARGV;
my ($param_file) = @ARGV;
print "dir with optimiste: ", $root_dir, "\n" if $DEBUG eq 'yes';
parse_param_file($param_file, \%params, $DEBUG); # PARSING INPUT FILE, INSERT PARAMETERS INTO %params
my ($template_file, @inits) = make_start_from_template($DEBUG); # MAKE START FILE FROM TEMPLATE
################################################################################################
my $grad_str ='';
$grad_str ="\$MPIRUN  grad_mpi | tee grad.log" if $params{'grad'} eq 'y' || $params{'grad'} eq 'all';
my $energy_str ="\$MPIRUN  dscf_mpi | tee xx.log";
#$energy_str ='' if $params{'without_energy'} eq 'y';
my $time_str ='';
if ($params{'type'} eq 'cluster') {
	$time_str ="#SBATCH --time=".$params{'time'} if exists $params{'time'};
}
$params{'ncpu'} ='' if !exists $params{'ncpu'};
$params{'tail'} ='' if !exists $params{'tail'};
$params{'scale_opt'} = '' if !exists $params{'scale_opt'};

my $task_message = <<"EOF";
#!/bin/sh
#SBATCH -n $params{'ncpu'} 
#SBATCH -o \%j.res
#SBATCH -e \%j.err
#SBATCH -p $params{'tail'}
$time_str

module load intel-compilers mpi

echo "CPU list: \$SLURM_JOB_NODELIST"

$energy_str
#rm realmos imagmos
$grad_str
EOF
################################################################################################
my @results = ();
foreach my $target_file (@target_files) {
	$target_file = Cwd::abs_path($target_file);
	print "Full path to target file: ", $target_file, "\n" if $DEBUG eq 'yes';
	my $target_dir = $target_file; 
	$target_dir =~ s/[^\/]+$//; 
	print "\n\n", $target_dir, "\n\n";
	if ($params{'scale_opt'} eq 'old') {
		change_scale($template_file, $target_file, \@inits, 'old');
	}
	elsif($params{'grad'} eq 'y') {
		change_scale($template_file, $target_file, \@inits, 'new');
	}
	else {
		make_in_file($template_file, $target_file, \@inits);
	}
	my $task_name = 'task';
	my $task_file = $target_dir.$task_name;
	print "Full path to task file: ", $task_file, "\n" if $DEBUG eq 'yes';
	make_task($task_file, $task_message) if $params{'type'} eq 'cluster';
	start_task($target_dir,$task_name) if $params{'type'} eq 'cluster';
	chdir $target_dir;
	my  $my_dir = getcwd;
	if ($params{'without_energy'} ne 'y') {
		system "dscf > xx.log" if $params{'type'} eq 'cpu';
	}
	system "grad > grad.log" if $params{'type'} eq 'cpu' ;
	chdir $my_dir;
	my $result_file = $target_dir.'xx.log';
	test_result($result_file,$params{'test'}) if $params{'without_energy'} ne 'y' && $params{'direct'} eq 'y';
	my $result = get_result($result_file,$params{'line'}, $params{'word'}) if $params{'without_energy'} ne 'y';
	push @results, $result  if $params{'without_energy'} ne 'y';
	my $control_file = $target_dir.'control';
	my $grad_file = $root_dir.'/xx.grad';
	wul_grad($template_file,$control_file,$grad_file) if $params{'grad'} eq 'old';
	all_grad($control_file,$grad_file) if $params{'grad'} eq 'y';
	if ($params{'scale_opt'} eq 'y') {
		my $gradient = '';
		open my $grad_read, "<", $grad_file or die "Can't open";
		while( my $st = <$grad_read>) {
                	$st =~s/[\r\n]+/ /g;
                	$st =~s/^\s+//g;
                	$gradient.=$st;
        	}
		close $grad_read;
		my @tmp_grad_ar = split /\s+/, $gradient;
		my $tmp_res = 0;
		foreach my $cur_grad (@tmp_grad_ar) {
			$cur_grad =~ s/D/e/g;
			$cur_grad = sprintf("%.10g", $cur_grad);
			$tmp_res += abs $cur_grad;
		}
		@results = ($tmp_res);
	}

}
chdir $root_dir;
open my $res_out, ">", 'xx.yy' or die "Can't open xx.yy";
print $res_out  scalar @results, "\n";
foreach my $res (@results) {
	print "qqqq\n";
	print $res_out $res, "\n";
}
close $res_out;

sub make_start_from_template {
	my ($DEBUG) = @_;
	my $template_file;
	($template_file, @target_files) = split /\s+/, $params{'file'};
	$template_file = Cwd::abs_path($template_file);
	print "Full path to template file: ", $template_file, "\n" if $DEBUG eq 'yes';
	open my $coord_file, "<", 'xx.xx' or die "Can't open  xx.xx\n";
	my $str_coord ="";
	while( my $st = <$coord_file>) {
		$st =~s/[\r\n]+/ /g;
		$st =~s/^\s+//g;
		$str_coord.=$st;
	}
	close $coord_file;
	print "NEW: ", $str_coord, "\n";
	$str_coord =~ s/^\s+//;
	my @inits = split /\s+/, $str_coord;
	if ($DEBUG eq 'yes') {
		print "Start parameters optimization: ";
	 	print $_ foreach @inits;
	 	print "\n";
	}
	return ($template_file, @inits);
}

sub change_scale {
	my ($template_file, $target_file, $ref2inits, $flag) = @_;
	print $flag, "\n";
	my @for_coord = ($params{'scale_init'}) if $flag eq 'old';
	make_in_file($template_file, 'xx.tmp', \@for_coord) if $flag eq 'old';
	my $in;
	open $in, "<", 'xx.tmp' or die "Can't open", 'xx.tmp', "\n" if $flag eq 'old';

	open $in, "<", $template_file or die "Can't open", 'xx.tmp', "\n" if $flag eq 'new';

	open my $out, ">", $target_file or die "Can't open", $target_file, "\n";

	while (my $str = <$in>) {
		if ($str =~ m/(\s+)([\d.eE-]+)(\s+)([\d.eE-]+)(\s+)([\d.eE-]+)(.*?)$/g){
			my $x = $2*$$ref2inits[0];
			my $y = $4*$$ref2inits[0];
			my $z = $6*$$ref2inits[0];
			$str = $1.$x.$3.$y.$5.$z.$7."\n";
		}
		print $out $str;
	}
	close $in;
	close $out;
	unlink 'xx.tmp' if $flag eq 'old';
}

sub parse_param_file {
	my ($in_file, $ref_params, $DEBUG) = @_;
	open my $in, "<", $in_file or die "Can't open file with parametres\n";
	print "PARSING INPUT FILE RESULTS\n" if $DEBUG eq 'yes';
	while (my $str = <$in>) {
		next if $str =~ m/^#/g ||$str =~ m/^[\r\n]+/g;
		my ($key, $value) =($1, $2) if $str =~ m/^(.*?)\s+(.*?)[\r\n]+/g;
		$ref_params -> {$key} = $value;
	}
	close $in;
	my %additional = ();
	%additional = so_params() if $params{'program'} eq 'dscf_so' ;
	%additional = no_rel_params() if $params{'program'} eq 'dscf_no_rel' ;
	# объединяем хеши
	foreach my $key (keys %additional) {
		$ref_params -> {$key} = $additional{$key} unless exists $ref_params -> {$key};
	}
	if($DEBUG eq 'yes') {
		foreach my $key (keys %$ref_params) {
			print $key, "\t", $ref_params -> {$key}, "\n";
		}
	}
	print "END PARSING INPUT FILE RESULTS\n" if $DEBUG eq 'yes';
	return 0; 
}

sub make_in_file { 
	# Делаю входной файл для задачи
	my ($template_file , $target_file, $ref2inits) = @_;
	open my $in, "<", $template_file or die "Can't open", $template_file, "\n";
	open my $out, ">", $target_file or die "Can't open", $target_file, "\n";
	while (my $str = <$in>) {
		$str =~s/_x(\d)_/$$ref2inits[$1-1]/g;
		print  $out  $str;
	}
	close $in;
	close $out;
	return 0; 
}
sub make_task {
# делаю батник с задачей
	my ($task_file, $task_message) = @_;
	open my $out, ">", $task_file or die "Can't open", $task_file, "\n";
	print $out $task_message;
	close $out;
	my $mode = 0700;   chmod $mode, $task_file;
}
sub start_task {
	my ($target_dir,$task_name) = @_;
	my  $my_dir = getcwd;
	chdir $target_dir;
	my $pid_str =`sbatch $task_name`;
	print  $pid_str, "\n";
	my $pid = $1 if $pid_str =~ m/Submitted batch job (\d+)$/;
	# Жду когда эта задача кончится
	while (1) {
		my $res ='';
        	$res = `squeue | grep  $pid`;
       		sleep(10);
        	return 0 if $res eq '';
	}
 	chdir $my_dir; 

}
sub test_result {
	# проверяю сходимость
	my ($result_file, $test_str) = @_;
	open my $in, "<", $result_file or die "Can't open", $result_file, "\n";
	while (my $str = <$in>) {
		return if $str =~ m/$test_str/g; 
	}
	close $in;
	print "test line not found \n";
	exit(2);

}
sub get_result {
        # ищу результат
        my ($result_file, $res_str, $word_num) = @_;
        open my $in, "<", $result_file or die "Can't open", $result_file, "\n";
        while (my $str = <$in>) {
		$str = strip($str);
                if ($str =~ m/$res_str/g) {
			my @tmp = split ' ', $str;
			print 'WORD: ', $tmp[$word_num], "\n" if $DEBUG eq 'yes';
			return $tmp[$word_num];
		}
        }
        close $in;
        print "get_result faled\n";
	exit 2;

}
sub help {
	print "Needed file with parameters for optimiste\n";
	print "FAILED\n";
	exit 1;
}
