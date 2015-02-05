#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use FindBin;
use lib "$FindBin::Bin";
use modules::SystemCall;
use modules::Exception;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "vep_bin=s",
	   "vep_extra=s",
	   "samtools=s",
	   "ref_fasta=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});



=pod

=head1 SYNOPSIS

test_vasp.pl -vep_bin test_with_local_vep_installation -vep_extra pass_extra_arguments_to_vep -samtools test_with_bam_list -fasta test_with_fasta [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

test_vasp.pl -> Test the installation of vasp using sample data

=head1 DESCRIPTION

Feb 5th, 2015

Wrapper script for testing VASP installatino

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

#Test simple case
>test_vasp.pl

#Test local install of vep
>test_vasp.pl -vep_bin /path/to/local/vep -vep_extra '--fork 4'

#Test using bam files (preferred)
>test_vasp.pl -samtools /path/to/local/samtools -ref_fasta /path/to/local/fasta

=cut

my $out = "test_vasp.tsv";
my $default_command = "./vasp.pl -vcf sample/sample.vcf -vep sample/sample.vep -ped sample/sample.ped -out $out";
my $sys_call = modules::SystemCall->new();

print "Testing VASP default....\n";
&test_vasp($out,$default_command,8743);

#Test with locally installed vep
if ($OPT{vep_bin}) {
	my $out_vep = "test_vasp_vep.tsv";
	my $vep_bin = $OPT{vep_bin};
	if (!-e $vep_bin) {
		modules::Exception->throw("ERROR: vep_bin $vep_bin binary doesn't exist");
	}
	if (!-x $vep_bin) {
		modules::Exception->throw("ERROR: vep_bin $vep_bin binary is executable");
	}
	my $vep_command = "./vasp.pl -vcf sample/sample.vcf -ped sample/sample.ped -vep_bin $vep_bin -out $out_vep";
	
	
	if ($OPT{vep_extra}) {
		$vep_command .= " -vep_extra \'$OPT{vep_extra}\'";
	} 
	print "\nTesting VASP with local VEP installation....\n";
	&test_vasp($out_vep,$vep_command,8743);
}

#Test with bam files
if ($OPT{samtools} && $OPT{ref_fasta}) {
	my $out_bam = "test_vasp_bam.tsv";
	my $samtools = $OPT{samtools};
	if (!-e $samtools) {
		modules::Exception->throw("ERROR: samtools $samtools binary doesn't exist");
	}
	if (!-x $samtools) {
		modules::Exception->throw("ERROR: vep_bin $samtools binary is executable");
	}
	my $fasta = $OPT{ref_fasta};
	if (!-e $fasta) {
		modules::Exception->throw("ERROR: fasta $fasta doesn't exist");
	}
	
	
	my $bam_command = "./vasp.pl -vcf sample/sample.vcf -ped sample/sample.ped -vep sample/sample.vep -ref_fasta $fasta -bam_list sample/bam_list -samtools $samtools -out $out_bam";
	print "\nTesting VASP with bam files....\n";
	&test_vasp($out_bam,$bam_command,8809);
}


print "\nInstallation Success!\n";

#subroutine to check output
sub test_vasp {
	my ($out,$command,$expected_lines) = @_;
	
	$sys_call->run($command);
	
	if ( !-e $out ) {
		modules::Exception->throw("ERROR: Installation failed; Output file $out wasn't generated");	
	} 

	if ( !-s $out) {
		modules::Exception->throw("ERROR: Installation failed; Output file $out is empty");
	}
	
	open(FILE,"$out") || modules::Exception->throw("Can't open file $out\n");
	my @lines = <FILE>;

	if (@lines != $expected_lines) {
		modules::Exception->throw("ERROR: Installation failed; Output file $out doesn't have correct number of lines");
	}
	close FILE;
}


