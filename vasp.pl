#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use FindBin qw($Bin);
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);
use lib catdir(dirname($Bin), 'perl');
use modules::VASP;
use modules::Exception;
use modules::PED;
use modules::SystemCall;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "ped=s",
	   "vcf=s",
	   "vep=s",
	   "bam_list=s",
	   "ref=s",
	   "samtools=s",
	   "out=s",
	   "vcf_cutoff=i",
	   "gene_list=s",
	   "coord=s",
	   "chrom=s",
	   "denovo", 
	   "min_phase_block=i", 
	   "phase_var_num=i", 
	   "inheritance=s", 
	   "comhet", 
	   "max_allele_freq=s",
	   "sift=s",
	   "polyphen=s",
	   "min_num_aff=i",
	   "debug"
	   ) || modules::Exception->throw("Invalid command-line option for script\n");

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{ped} || !$OPT{vcf} || !$OPT{vep});



=pod

=head1 SYNOPSIS

vasp.pl -ped ped_file -vcf input_vcf -vep vep_annotation_file

[OPTIONAL ARGS]

-bam_list list_containing_id_and_bamfile_locations (preferred method)

-ref ref_fasta_location (required with -bam_list)

-samtools samtools_path (required with -bam_list; default=/usr/bin/samtools)

-vcf_cutoff quality_cutoff(default=20, ignored if no quality) 

-min_num_aff min_number_of_affected_samples_variant

-polyphen filter_on_polyphen_category(default=probably_damaging)

-sift filter_on_sift_category(default=deleterious) 

-comhet only_report_comhet 

-max_allele_freq filter_on_allele_freq 

-phase_var_num min_number_variants_to_make_block(default=2) 

-min_phase_block min_phase_block_size(default=10kb)

-denovo only_report_denovo

-inheritance filter_on_inheritance_type(options=ad,ar,xd,xr) 

-coord genomic_coordinates(format=chr:start-end)

-gene_list only_report_on_gene_in_file 

-chrom only_report_chr   

-out outfile(default=./vasp.tsv) 

Required flags: -ped -vep -vcf

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

vasp.pl -> wrapper script for generating pedigree summaries

=head1 DESCRIPTION

Oct 31/2014

A script that summarises seqeuenced pedigrees. Takes in three required files (ped, vcf, and vep file) and many optional filters.

=head1 AUTHOR

Matthew Field

See https://github.com/mattmattmattmatt/VASP

=cut



#Read required arguments and check files exist
my $vcf_file = $OPT{vcf};

if ( !-e $vcf_file ) {
	modules::Exception->throw("File $vcf_file doesn't exist");	
}

#Flag to determine if filtering is required
my $filter_flag = 0;

my $ped_file = $OPT{ped};

if ( !-e $ped_file ) {
	modules::Exception->throw("File $ped_file doesn't exist");	
}

my $vep_file = $OPT{vep};

if ( !-e $vep_file ) {
	modules::Exception->throw("File $vep_file doesn't exist");	
}
my $out = defined $OPT{out}?$OPT{out}:'vasp.tsv';
$out .= '.tsv' unless $out =~ /tsv$/;

my $cwd = `pwd`;
chomp $cwd;

#default args for object
my %args = (
			-vcf_file=>$vcf_file,
			-ped_file=>$ped_file,
			-vep_file=>$vep_file,
			-out=>$out,
			);


if ($OPT{debug}) {
	$args{-debug} = 1;
}


my %bam_data = ();
my $allele_data;

if (defined $OPT{bam_list}) {
	if ( !-e $OPT{bam_list} ) {
		modules::Exception->throw("File $OPT{bam_list} doesn't exist");	
	}
	
	if (!defined $OPT{ref}) {
		modules::Exception->throw("ERROR: Need to pass in reference fasta with bam argument");
	}
	
	if ( !-e $OPT{ref} ) {
		modules::Exception->throw("File $OPT{ref} doesn't exist");	
	}
	
	$args{-ref} = $OPT{ref};
	
	
	open(BAMLIST,"$OPT{bam_list}") || modules::Exception->throw("Can't open file $OPT{bam_list}\n");
	
	while (<BAMLIST>) {
		chomp;
		my @fields = split;
		if (! -e $fields[1]) {
			modules::Exception->throw("ERROR: File $fields[1] doesn't exist; please use full path");
		}
		
	}
	$allele_data = 'bam';
	$args{-bam_list} = $OPT{bam_list};
	
	my $samtools;
	if ($OPT{samtools}) {
		$samtools = $OPT{samtools};
	} else {
		$samtools = '/usr/bin/samtools';
	}
	
	if (!-e $samtools) {
		modules::Exception->throw("ERROR: Samtools $samtools doesn't exists");
	}
	
	
	$args{-samtools} = $samtools;
	
} else {
	$allele_data = 'vcf';
}

#Either get allele info from bam files via pileups or from vcf files
$args{-allele_data} = $allele_data;

if (defined $OPT{vcf_cutoff}) {
	$args{-vcf_cutoff} = $OPT{vcf_cutoff};
	$filter_flag = 1;
}

#Now deal with all the optional filtering 

#Filtering by gene
if (defined $OPT{gene_list}) {
	my %genes = ();
	open(GENES,$OPT{gene_list}) || modules::Exception->throw("Can't open file $OPT{gene_list}\n");
	while (<GENES>) {
		next unless /\S/;
		my ($gene_name) = $_ =~ /(\S+)/;
		$genes{$gene_name} = 1;
	}
	$args{-gene_list} = \%genes;
	$filter_flag = 1;
}

if (defined $OPT{chrom} && defined $OPT{coord}) {
	modules::Exception->throw("ERROR: Only filter by -chrom or -coord, not both");
}

#Filter by chrom
if (defined $OPT{chrom}) {
	if ($OPT{chrom} !~ /^c?h?r?[0-9XYMT]+$/) {
		modules::Exception->throw("ERROR: -chrom argument not formatted properly");
	}
	$args{-chrom} = $OPT{chrom};
	$filter_flag = 1;
}

#Filter by chrom
if (defined $OPT{coord}) {
	if ($OPT{coord} =~ /[0-9XY]+:\d+\-\d+/) {
		$args{-coord} = $OPT{coord};
	} elsif ($OPT{coord} =~ /[0-9XY]+:(\d+)/) {
		my $coord = $OPT{coord} .'-'.$1;
		$args{-coord} = $coord;		
	} else {
		modules::Exception->throw("ERROR: -coord argument must be of the form chr:start-end OR chr:coord");
	}
	$filter_flag = 1;
}

if (defined $OPT{min_phase_block}) {
	$args{-phase_block_size} = $OPT{min_phase_block}
} 

if (defined $OPT{phase_var_num}) {
	$args{-phase_var_num} = $OPT{phase_var_num};
} 


if (defined $OPT{denovo} && $OPT{inheritance} ) {
	modules::Exception->throw("ERROR: Only filter by -denovo or -inheritance. not both");
}

if (defined $OPT{denovo}) {
	$args{-denovo} = 1;
	$filter_flag = 1;
}

my @dis_inh_options = qw(ad ar xd xr);
if (defined $OPT{inheritance}) {
	my $input = lc($OPT{inheritance});
	
	my $match = 0;
	for my $dis_type ( @dis_inh_options ) {
	    if ($dis_type eq $input) {
	    	$match = 1;
	    }
	}
	
	
		
	if (!$match) {
		modules::Exception->throw("ERROR: problem with inheritance; must be ad,ar,xd,or xr");
	} else {
		my %inheritance_map = (
								ad => 'auto-dominant',
								ar => 'auto-recessive',
								xd => 'x-dominant',
								xr => 'x-recessive'
							  );
		$args{-inheritance} = $inheritance_map{$input};
	}
	$filter_flag = 1;
}

if ($OPT{comhet}) {
	$args{-comhet} = 1;
	$filter_flag = 1;
}



if ($OPT{min_num_aff}) {
	$args{-min_num_aff} = $OPT{min_num_aff};
	$filter_flag = 1;
}


#First check vep was run with canonical flag
if (! `grep CANONICAL $vep_file | head -1`) {
	modules::Exception->throw("ERROR: Can't use vep file as canonical genes not annotated; must rerun with --canonical flag");
}


#These filters are in vep annotation files so check they exist
if ($OPT{max_allele_freq}) {
	if ($OPT{max_allele_freq} >= 0 && $OPT{max_allele_freq} <= 1) {
		$args{-max_allele_freq} = $OPT{max_allele_freq};
	} else {
		modules::Exception->throw("ERROR: max_allele_freq must be between 0 and 1");
	}
	if (! `grep GMAF $vep_file | head -1`) {
		modules::Exception->throw("ERROR: Can't filter on max_allele_freq as vep file doesn't contain gmaf info; must rerun with --gmaf flag");
	}
	$filter_flag = 1;
}

if ($OPT{polyphen}) {
	my $match = 0;
	my @polyphen_options = qw(probably_damaging possibly_damaging benign);
	for my $polyopt ( @polyphen_options ) {
	    if ($polyopt eq $OPT{polyphen}) {
	    	$match = 1;
	    	$args{-polyphen} = $OPT{polyphen};
	    }
	}
	if (!$match) {
		modules::Exception->throw("ERROR: polyphen input must be probably_damaging,possibly_damaging, or benign");
	}
	if (! `grep PolyPhen $vep_file | head -1`) {
		modules::Exception->throw("ERROR: Can't filter on polyphen as vep file doesn't contain polyphen info; must rerun with '--polyphen b' flag");
	}
	$filter_flag = 1;
}

if ($OPT{sift}) {
	my @sift_options = qw(deleterious tolerated);
	my $match = 0;
	for my $siftopt ( @sift_options ) {
	    if ($siftopt eq $OPT{sift}) {
	    	$match = 1;
	    	$args{-sift} = $OPT{sift};
	    }
	}
	if (!$match) {
		modules::Exception->throw("ERROR: sift input must be tolerated or deleterious");
	}
	
	if (! `grep SIFT $vep_file | head -1`) {
		modules::Exception->throw("ERROR: Can't filter on SIFT as vep file doesn't contain SIFT info; must rerun with '--sift b' flag");
	}
	$filter_flag = 1;
}

my $vasp = modules::VASP->new(%args);

if ($allele_data eq 'bam') {
	print "Generate pileups\n";
	$vasp->generate_pileups();
} 

#Calculates inheritance, etc; different behavior if bams available
print "Get inheritance info\n";
$vasp->get_inheritance_info();
	
print "Get Com het\n";
$vasp->get_compound_het();
	
print "Get Phasing\n";
$vasp->get_phasing_blocks();

if ($filter_flag) {
	print "Filtering results\n";
	$vasp->filter_results();
}

print "Get line data\n";
$vasp->generate_line_data();

print "Write files\n";
$vasp->write_to_files();



