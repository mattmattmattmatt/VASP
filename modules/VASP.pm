package modules::VASP;

use strict;
use Data::Dumper;
use modules::Exception;
use modules::SystemCall;
use modules::Utils;
use File::Basename;

sub new {
	my ($class, @args) = @_;
	
	my $self = bless {}, $class;

    my %args = @args;


    my @required_args = (
    					-vcf_file,
    					-ped_file,
    					-vep_file,
    					-out,
    					-allele_data
    					);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set $args{$required_arg}");
		} 
    }
   
   	#Load the bam file info
   	if (exists $args{-bam_list}) {
   		open(BAM,$args{-bam_list}) || modules::Exception->throw("Can't open file $args{-bam_list}\n");
   		while (<BAM>) {
   			my ($sample_id,$bam_file) = split;
   			$self->{sample_data}{$sample_id}{bam} = $bam_file;
   		}
   		$self->{samtools} = $args{-samtools};
   		$self->{ref} = $args{-ref};
   	}
   
   	
   	#Determines whether to use bam files to generate pileups or the vcf file when obtaining allele information
   	$self->{allele_data} = $args{-allele_data};
   	
   	if ($args{-out}  =~ /\//) {
   		#If it's full directory
	   	$self->{out} = $args{-out};
   	} else {
   		$self->{out} = './'.$args{-out};
   	}
   	
   	
   	
   	#Now parse the ped file
   	my $ped = modules::PED->new();
   	if (my $error = $ped->validate_ped(-ped_file=>$args{-ped_file})) {
   		modules::Exception->throw("ERROR: Problem with ped file $error");
   	}
	
	
	#Load the sample info into the VASP object   	
   	my $ped_data = $ped->parse_ped(-ped_file=>$args{-ped_file});
	$self->{ped_data} = $ped_data;

	my %parent_mapping = ();

	
	my %parent_count = ();

   	for my $ped_identifier (sort {my $a_count = split(":",$a); my $b_count = split(":",$b); $a_count<=>$b_count} keys %{$ped_data} ) {
		my ($ped_count,$ped_id) = split(":",$ped_identifier);
		$self->{sample_data}{$ped_id}{count} = $ped_count;


		#Generate total affected and unaffected number as well as parent info
		for my $ped_key (keys %{$ped_data->{$ped_identifier}}) {
			my $mother = 0;
			my $father = 0;
			if ($ped_key eq 'mother' || $ped_key eq 'father') {
				$parent_mapping{$ped_id}{$ped_key} = $ped_data->{$ped_identifier}{$ped_key};
				$parent_count{$ped_id}++;
			} 
			if ($ped_key eq 'affected') {
				$self->{total_affected}++;
			}
			if ($ped_key eq 'unaffected') {
				$self->{total_unaffected}++;
			}
			
			$self->{sample_data}{$ped_id}{$ped_key} = $ped_data->{$ped_identifier}{$ped_key};
		}
   	}
   	
	#Check we have two parents for one pedigree member if filtering by denovo
	if (exists $args{-denovo}) {
		my $two_parents = 0;
		for my $sample ( keys %parent_count ) {
		    $two_parents = 1 if $parent_count{$sample} == 2;
		}
		if (!$two_parents) {
			modules::Exception->throw("ERROR: Can't filter for denovo mutations without two parents");
		}
	}
   	
   	if (keys %parent_mapping) {
   		$self->{parent} = 1;
   		$self->{parent_mapping} = \%parent_mapping;
   	}
   	
   	$self->{total_samples} = keys %{$ped_data};

   	
   	my $vcf_cutoff;
   	if (exists $args{-vcf_cutoff}) {
   		$vcf_cutoff = $args{-vcf_cutoff};
   	} else {
   		$vcf_cutoff = 20;
   	}
   	
	
   	
   	#Parse the optional output filters
   	my $filter_output = 0;
   	
   	if (exists $args{-chrom}) {
   		$filter_output = 1;
   		$self->{filters}{chrom} = $args{-chrom};
   	}
   	
   	if (exists $args{-coord}) {
   		$filter_output = 1;
   		my ($chr,$start,$end) = $args{-coord} =~ /([0-9XYMT]+):(\d+)\-(\d+)/;
   		$self->{filters}{chrom} = $chr; #overwrite $args{-chr} if also passed in
   		$self->{filters}{start} = $start;
   		$self->{filters}{end} = $end;
   	}
   	
   	if (exists $args{-gene_list}) {
   		$filter_output = 1;
   		$self->{filters}{gene_list} = $args{-gene_list};
   	}
   	
   	if (exists $args{-denovo}) {
   		$filter_output = 1;
   		$self->{filters}{denovo} = 1;
   	} 
   	 	
   	if (exists $args{-phase_block_size}) {
   		$filter_output = 1;
   		$self->{filters}{phase_block_size} = $args{-phase_block_size};
   	} 
   	
   	if (exists $args{-phase_var_num}) {
   		$filter_output = 1;
   		$self->{filters}{phase_var_num} = $args{-phase_var_num};
   	} 
   	
   	if (exists $args{-inheritance}) {
   		$filter_output = 1;
   		$self->{filters}{inheritance} = $args{-inheritance};
   	} 
   	
   	if (exists $args{-max_allele_freq}) {
   		$filter_output = 1;
   		$self->{filters}{max_allele_freq} = $args{-max_allele_freq};
   	} 
   	
   	if (exists $args{-min_num_aff}) {
   		$filter_output = 1;
   		$self->{filters}{min_num_aff} = $args{-min_num_aff};
   	} 
   	
   	if (exists $args{-comhet}) {
   		$filter_output = 1;
   		$self->{filters}{comhet} = 1;
   	}
   	
   	if (exists $args{-sift}) {
   		$filter_output = 1;
   		$self->{filters}{sift} = $args{-sift};
   	} 
   	
   	if (exists $args{-polyphen}) {
   		$filter_output = 1;
   		$self->{filters}{polyphen} = $args{-polyphen};
   	} 
   	
   	if ($filter_output) {
   		$self->{filter_output} = 1;
   	}
   	
   	if (exists $args{-debug}) {
   		$self->{debug} = 1;
   	}
   	
   	if ($self->{debug}) {
	   	print Dumper $self;
   	}
   	
   	print "Parsing vep...\n";
   	$self->_parse_vep(-vep_file=>$args{-vep_file});
   	print "Parsing vcf...\n";   	
   	$self->_parse_vcf(-vcf_file=>$args{-vcf_file},-vcf_cutoff=>$vcf_cutoff);
   	
    return $self;
}

sub _parse_vcf {
    my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file,
			             -vcf_cutoff
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    my $vcf_cutoff = $args{-vcf_cutoff};
    
    open(VCF,$args{-vcf_file}) || modules::Exception->throw("Can't open file $args{-vcf_file}\n");
    
    while (<VCF>) {
    	next if /^#/;
    	next unless /\w/;
    	chomp;
    	
    	
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest,$gt_fields,@ped_alleles) = split;
    	
    	if ($chr eq 'MT') {
			$chr = 'M';
		}
    	
    	if (@ped_alleles != $self->{total_samples}) {
    		my $ped_data_count = @ped_alleles;
    		my $ped_total = $self->{total_samples};
    		modules::Exception->throw("ERROR: Expecting $ped_total sample and vcf only has data for $ped_data_count samples");	
    	}
    	

    	if ($qual !~ /\d/ && $qual ne '.') {
    		modules::Exception->throw("ERROR: Error with qual $qual format at line $_");
    	}

		if ($ref =~ /N/ || $var_str =~ /N/) {
			next;
		}

		if ($var_str eq '.') {
			next;
		}

		my @vars = split(",",$var_str);
		
		for my $var ( @vars ) {
			my ($var_key,$var_type) = _get_variant_key(-type=>'vcf',-chrom=>$chr,-first=>$first_coord,-ref_seq=>$ref,-var_seq=>$var);


			my ($start,$end) = $var_key =~ /(\d+)-(\d+)/;

			#Skip by coordinates if filter set
			if ($self->{filter_output}) {
				if (exists $self->{filters}{chrom}) {
					next unless $self->{filters}{chrom} eq $chr;
				}
				
			}
			
			if (!exists $self->{var_data}{$var_key}) {
				#Skip these intergenic cases as won't have canonical vep annotations
				next;
			}

			#If fails quality test; only apply if quality score available
			if ($qual ne '.' && $qual <= $vcf_cutoff) {
				delete $self->{var_data}{$var_key};
				next;
			}


			my $length_ref = length($ref);
			my $length_var = length($var);
			
			if ($length_ref == $length_var && $length_ref != 1) {
				modules::Exception->warning("Skip vcf entry $ref -> $var");
			}
			
			if ($qual eq '.') {
				$self->{var_data}{$var_key}{qual} = "N/A";				
			} else {
				$self->{var_data}{$var_key}{qual} = $qual;
			}
			$self->{var_data}{$var_key}{var_type} = $var_type;

			
			$self->{var_data}{$var_key}{aff_count} = 0;
			$self->{var_data}{$var_key}{unaff_count} = 0;
			
			if ($self->{allele_data} eq 'vcf') {
				#Need to use the sorted ped information here
				for my $ped_key (sort {my $a_count = split(":",$a); my $b_count = split(":",$b); $a_count<=>$b_count} keys %{$self->{ped_data}}) {
					my $zyg;
					my ($ped_count,$ped_id) = split(":",$ped_key);
					my $affected = $self->{ped_data}{$ped_key}{affected}?1:0;
					my $ped_lookup = $ped_count-1;
					if ($ped_alleles[$ped_lookup] eq '.') {
						$zyg = 'ref';
					} else {
						my ($alleles) = split(':',$ped_alleles[$ped_lookup]);
						my ($allele1,$allele2) = split('/',$alleles);
						if ($alleles eq '0/0' || $alleles eq '.') {
							$zyg = 'ref';
						} elsif ($allele1 eq $allele2) {
							$zyg = 'hom';
							if ($affected) {
								$self->{var_data}{$var_key}{aff_count}++;
							} else {
								$self->{var_data}{$var_key}{unaff_count}++;
							}
						} elsif ($allele1 ne $allele2) {
							$zyg = 'het';
							if ($affected) {
								$self->{var_data}{$var_key}{aff_count}++;
							} else {
								$self->{var_data}{$var_key}{unaff_count}++;
							}
						}  else {
							modules::Exception->throw("ERROR: Can't parse zygosity from $alleles for line $_");
						}
					}
					$self->{var_data}{$var_key}{zyg}{$ped_id} = $zyg;
				}
			}
		}
    }
    
    $self->_check_vars_remain();
}

sub _parse_vep {
    my ($self,@args) = @_;
    
    my @required_args = (
			             -vep_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    open(my $OUTPUT, $args{-vep_file}) || modules::Exception->throw("ERROR: Can't open output file $args{-vep_file}");

	my %grouping_data = ();

	my %aa_found = ();

	while (<$OUTPUT>){
		chomp;
		#Skip headers
		next if $_ =~ /^#/;
		next unless $_ =~ /\w/;
		
		my ($identifier, $coord_str, $var_base, $ens_gene, $ens_transcript, $classifier, $aa_type, undef, undef, undef, $aa_change, $codon_change, $rs, $attribute_str ) = split /\t/;
		
		#Only record results from canonical transcripts
		next unless $attribute_str =~ /CANONICAL/; 
		
		if ($identifier =~ /N/) {
			next;
		}
		
		my $chr = my $start = my $end;
		if ($coord_str =~ /([0-9XYMT]+):(\d+)\-(\d+)/) {
			$chr = $1;
			$start = $2;
			$end = $3;
		} else {
			($chr,$start) = $coord_str =~ /([0-9XYMT]+):(\d+)/;
			$end = $start;
		}
		
		if ($chr eq 'MT') {
			$chr = 'M';
		}
		
		if ($self->{filter_output} && exists $self->{filters}{chrom}) {
			#Save filtering for the end to allow for comhet and phase blocks stats to be accurate
			#Ok to filter here as above measures will occur within single chromosome
			next unless $self->{filters}{chrom} eq $chr;
		}
		
		
		my $ref_base;
		my $var_identifier;
		
		my $first_coord = my $ref;
		
		if ($identifier =~ /(\S+)_(\S+)_(\S+)/) {
			$first_coord = $2;
			($ref) = split('/',$3);
		} else {
			modules::Exception->throw("ERROR: Can't parse vep identifier $identifier");
		}
			
		my ($var_key,$var_type) = _get_variant_key(-type=>'vep',-chrom=>$chr,-first=>$first_coord,-ref_seq=>$ref,-var_seq=>$var_base);		
		
		if (length($var_base) == length($ref) && length($ref) != 1) {
			modules::Exception->warning("Skip entry vep $identifier");
		}
		
		#update ens gene data
		$self->{gene_data}{$ens_gene}{uniq_count}{$var_key}++;
		
		
		#Check general variant fields first
		my $allele_freq = "N/A";  
    	$self->{var_data}{$var_key}{allele_freq} = $allele_freq;
    	
    	
    	my @attribute_pairs = split /;/, $attribute_str;
		my $exon = my $intron = 0;
		
		
		
		
		if ($rs =~ /rs/) {
			$self->{var_data}{$var_key}{dbsnp} = $rs;
		}
		
		my $fail_allele_freq = 0;
		
		for my $attribute_pair ( @attribute_pairs ) {
		    
		   if ($attribute_pair =~ /GMAF=(.*)/) {
		   		my @entries = split(",",$attribute_pair);
		   	
		   		for my $entry ( @entries ) {
			    	my ($base,$freq) = split(':',$entry);
			    	$base =~ s/GMAF=//;
			    	if ($base eq $var_base) {
			    		#Here the variant allele is more common so invert the frequency to approximate; needed if filtering
			    		$allele_freq = $freq;
			    	} else {
			    		$allele_freq = 1-$freq;
			    	}
		   		}
		   	
			    $self->{var_data}{$var_key}{allele_freq} = $allele_freq;
		    	
		    	
		   }
		    
		   if ($attribute_pair =~ /HGNC=(.*)/) {
		   		#update hgnc gene info
		   		$self->{var_data}{$var_key}{hgnc} = $1;
		   		$self->{gene_data}{$ens_gene}{hgnc} = $1;
		   } 
		    
		}
		
		
		#Check exon specific fields 	

		if ($aa_change =~ /([A-Z\*])\/([A-Z\*])/){
			my $aa_ref = $1;
			my $aa_var = $2;
			my $polyphen_pred = "N/A";
			my $polyphen_score = "N/A";
			my $sift_pred = "N/A";	
			my $sift_score = "N/A";
		
			my @attribute_pairs = split /;/, $attribute_str;
			#Sample attribute line: 
			#PolyPhen=possibly_damaging(0.593);CANONICAL=YES;SIFT=deleterious(0);EXON=5/11
			for my $attribute_pair ( @attribute_pairs ) {			    
			    if ($attribute_pair =~ /PolyPhen/) {
			    	($polyphen_pred, $polyphen_score) = $attribute_pair =~ /PolyPhen=(.+)\(([0-9\.]+)\)/;
					$self->{var_data}{$var_key}{poly_pred} = $polyphen_pred;
					$self->{var_data}{$var_key}{poly_score} = $polyphen_score;
					
					
			    } elsif ($attribute_pair =~ /SIFT/) {
			    	($sift_pred, $sift_score) = $attribute_pair =~ /SIFT=(.+)\(([0-9\.]+)\)/;
			    	$self->{var_data}{$var_key}{sift_pred} = $sift_pred;
			    	$self->{var_data}{$var_key}{sift_score} = $sift_score;
			    	
			    } 
			}
		
			
		
		
			$aa_ref = 'Stop' if ($aa_ref eq '*');
			$aa_var = 'Stop' if ($aa_var eq '*');
			
			my $aa_string = "$aa_ref->$aa_var";
			#This case trumps any previous annotations
			$self->{var_data}{$var_key}{aa_change} = $aa_change;
			$self->{var_data}{$var_key}{transcript} = $ens_transcript;
			$self->{var_data}{$var_key}{gene} = $ens_gene;
			$self->{var_data}{$var_key}{vep} = $aa_type;
			$aa_found{$var_key}++;
			
			
		} 
		
		#Only add gene info if not already annotated with non-syn change
		if (!exists $aa_found{$var_key}) {
			if ($ens_gene =~ /ENS/) {
				$self->{var_data}{$var_key}{gene} = $ens_gene;
			}
			if ($ens_transcript =~ /ENS/) {
					$self->{var_data}{$var_key}{transcript} = $ens_transcript;
			}
			if (!exists $aa_found{$var_key}) {
				$self->{var_data}{$var_key}{vep} = $aa_type;
			}
		}	
		
		
	}
	
	
    close($OUTPUT);
}

#subroutine creates the pileup files for getting allele strings later
sub generate_pileups {
	my ($self) = @_;
	
	my $pileup_coord_file = 'pileup.coord';
	open(COORD,">$pileup_coord_file") || modules::Exception->throw("Can't open file to write pileup.coord\n");
	
	my %pileup_lookup = ();
	#First generate the pileup coords file
	for my $var_key (sort {my ($a_chr,$a_start) = $a =~ /([0-9XYMT]+):(\d+)/; my ($b_chr,$b_start) = $b =~ /([0-9XYMT]+):(\d+)/; $a_chr gt $b_chr || $a_start <=> $b_start } keys %{$self->{var_data}}) {
		my ($chr,$start) = $var_key =~ /([0-9XYMT]+):(\d+)/;
		my $pileup_coord = $start;
		if ($self->{var_data}{$var_key}{var_type} eq 'DEL') {
			$pileup_coord--;		
		} 
		print COORD "$chr\t$pileup_coord\n";
			
		
		$pileup_lookup{"$chr:$pileup_coord"} = $var_key;
	}
	my $sys_call = modules::SystemCall->new();
	my $ref = $self->{ref};
	
	
	for my $sample (keys %{$self->{sample_data}}) {
		my $bam_file = $self->{sample_data}{$sample}{bam};
		my $bamdir = dirname($bam_file);
		
		my $pileup_file = $bamdir . '/' . $sample . '.pileup';
		my $samtools_bin = $self->{samtools};
		
		my $mpileup_snv_command = "$samtools_bin mpileup -A -E  -l $pileup_coord_file -f $ref $bam_file  > $pileup_file";
		$sys_call->run($mpileup_snv_command) unless -e $pileup_file; #Don't run longish command if already run
		open(PILEUP,"$pileup_file") || modules::Exception->throw("Can't open pileup file $pileup_file");
		while (<PILEUP>) {
			my @fields = split("\t");
			my ($pileup_string,$zyg) = modules::Utils->pileup_string($fields[4]);
			
			my $var_key_lookup = $pileup_lookup{$fields[0].':'.$fields[1]};	
				
			if (!exists $pileup_lookup{$fields[0].':'.$fields[1]}){
				if (exists $self->{filter_output} && exists $self->{filters}{chrom}) {
					#Can happen when use previously generated pileup file and rerun with genome coordinate filters
					next;					
				} else {
					#Here the is a problem as there is no information from the var_data field
					modules::Exception->throw("ERROR: Problem with pileup line for coord $fields[0] $fields[1]");
				}
				next; 
			}
				
			$self->{var_data}{$var_key_lookup}{pileup_str}{$sample} = $pileup_string;
			if ($zyg ne 'ref') {
				if ($self->{sample_data}{$sample}{affected}) {
					$self->{var_data}{$var_key_lookup}{aff_count}++;
				} elsif ($self->{sample_data}{$sample}{unaffected}) {
					$self->{var_data}{$var_key_lookup}{unaff_count}++;
				}
			}
			$self->{var_data}{$var_key_lookup}{zyg}{$sample} = $zyg;
		}
		
		#Set the no data entries
		for my $var_key (keys %{$self->{var_data}}) {
			if (!exists $self->{var_data}{$var_key}{pileup_str}{$sample}) {
				$self->{var_data}{$var_key}{pileup_str}{$sample} = "No data";
			}	
			
			if (!exists $self->{var_data}{$var_key}{zyg}{$sample}) {
				$self->{var_data}{$var_key}{zyg}{$sample} = "No data";
			}
		}
		
	}
	
   	
}	


#subroutine gets the raw parent alleles and sets the mendelian flag	
sub get_inheritance_info {
	
	my ($self) = @_;
	
	my %parent_mapping = ();
	if ($self->{parent}) {
		%parent_mapping = %{$self->{parent_mapping}};
	}

	#Iterate through variants and populate parent_allele info
	for my $var_key (keys %{$self->{var_data}}) {
		
		my $affected_alleles_from_one_parent = 'yes'; #default is yes; set to zero if misproven or ? if alleles unknown
		my %affected_allele_data = ();
		
		#Get the comhet key for later
		my $het_key = $var_key . ':GMAF='. $self->{var_data}{$var_key}{allele_freq};
		
		#Flag for deleting comhet entries if conditions are broken
		my $del_com_het = 0;
		
		my $gene = $self->{var_data}{$var_key}{gene};
		
		
		for my $local_sample (keys %{$self->{var_data}{$var_key}{zyg}}) {
			my $cz = $self->{var_data}{$var_key}{zyg}{$local_sample}; #child zygosity
			if (exists $parent_mapping{$local_sample}) {
				if (exists $parent_mapping{$local_sample}{mother} && exists $parent_mapping{$local_sample}{father}) {
					#Both parents available; easy case					
					my $mother = $parent_mapping{$local_sample}{mother};
					my $mz = $self->{var_data}{$var_key}{zyg}{$mother};
					my $father = $parent_mapping{$local_sample}{father};
					my $fz = $self->{var_data}{$var_key}{zyg}{$father};
					my $mother_first_allele = my $mother_second_allele = my $father_first_allele = my $father_second_allele;
					
					if ($self->{allele_data} eq 'bam') {
						#Here we use the pileup_str which is preferred
						($mother_first_allele,$mother_second_allele) = $self->_get_alleles_from_pileup($mz,$self->{var_data}{$var_key}{pileup_str}{$mother});
						($father_first_allele,$father_second_allele) = $self->_get_alleles_from_pileup($fz,$self->{var_data}{$var_key}{pileup_str}{$father});						
					} else {
						($mother_first_allele,$mother_second_allele) = $self->_get_alleles_from_varkey($cz,$mz,$var_key);
						($father_first_allele,$father_second_allele) = $self->_get_alleles_from_varkey($cz,$fz,$var_key);						
					}
					my $inherited_mother_allele = '?';
					my $inherited_father_allele = '?';
					
					if (($mz eq 'ref' && $fz eq 'ref' && $cz eq 'ref') || ($mz eq 'het' && $fz eq 'ref' && $cz eq 'ref') || ($mz eq 'ref' && $fz eq 'het' && $cz eq 'ref') || ($mz eq 'het' && $fz eq 'het' && $cz eq 'ref')) {
						#ref/ref cases
						$inherited_mother_allele = $mother_first_allele;
						$inherited_father_allele = $father_first_allele;
					} elsif (($mz eq 'hom' && $fz eq 'het' && $cz eq 'hom') || ($mz eq 'het' && $fz eq 'hom' && $cz eq 'hom') || ($mz eq 'hom' && $fz eq 'hom' && $cz eq 'hom') || ($mz eq 'het' && $fz eq 'het' && $cz eq 'hom')) {
						#hom/hom cases
						$inherited_mother_allele = $mother_second_allele;
						$inherited_father_allele = $father_second_allele;
					} elsif (($mz eq 'ref' && $fz eq 'het' && $cz eq 'het') || ($mz eq 'ref' && $fz eq 'hom' && $cz eq 'het') || ($mz eq 'het' && $fz eq 'hom' && $cz eq 'het')) {
						#father variant cases
						$inherited_mother_allele = $mother_first_allele;
						$inherited_father_allele = $father_second_allele;
		
						if ($self->_affected($local_sample)) {
							$affected_allele_data{father}++;
						} else {
							$affected_alleles_from_one_parent = 'No';
						}
					} elsif (($mz eq 'het' && $fz eq 'ref' && $cz eq 'het') || ($mz eq 'hom' && $fz eq 'ref' && $cz eq 'het') || ($mz eq 'hom' && $fz eq 'het' && $cz eq 'het')) {
						#mother variant cases
						$inherited_mother_allele = $mother_second_allele;
						$inherited_father_allele = $father_first_allele;

						if ($self->_affected($local_sample)) {
							$affected_allele_data{mother}++;
						} else {
							$affected_alleles_from_one_parent = 'No';
						}
					} 
					
					
					if ($cz eq 'het') {
						
						#Flag if it's definite compound het 
						if (($mz eq 'het' && $fz eq 'ref') || ($mz eq 'hom' && $fz eq 'ref') || ($mz eq 'hom' && $fz eq 'het')) {
							#mother variant allele
							$self->{gene_data}{$gene}{comhet}{$local_sample}{mother}{$het_key}++;
						} 
						
						if (($fz eq 'het' && $mz eq 'ref') || ($fz eq 'hom' && $mz eq 'ref') || ($fz eq 'hom' && $mz eq 'het')) {
							#father variant allele
							$self->{gene_data}{$gene}{comhet}{$local_sample}{father}{$het_key}++;
						}
		
						if ($fz eq 'het' && $mz eq 'het') {
							#ambiguous case
							$self->{gene_data}{$gene}{comhet}{$local_sample}{either}{$het_key}++;
						}
					} 
					
					

					if ($inherited_father_allele eq '?' || $inherited_mother_allele eq '?') {
						$affected_alleles_from_one_parent = 'No';
					}
					
					#If unaffected and not reference or unknown
					if ($self->_unaffected($local_sample)) {
						if ($inherited_father_allele ne 'ref' || $inherited_mother_allele ne 'ref') {
							$affected_alleles_from_one_parent = 'No';
						}
					}
					
					#If affected and is reference or unknown
					if ($self->_affected($local_sample)) {
						if ($inherited_father_allele eq 'ref' && $inherited_mother_allele eq 'ref') {
							$affected_alleles_from_one_parent = 'No';
						}
					}
					
					$self->{var_data}{$var_key}{alleles}{$local_sample}{mother} = $inherited_mother_allele;
					$self->{var_data}{$var_key}{alleles}{$local_sample}{father} = $inherited_father_allele;
					
					
				} elsif (exists $parent_mapping{$local_sample}{mother} || exists $parent_mapping{$local_sample}{father}) {
					#Single parent case
					my $parent_name = my $parent;
					if (exists $parent_mapping{$local_sample}{mother}) {
						$parent_name = 'mother';
						$parent = $parent_mapping{$local_sample}{mother};
					} else {
						$parent_name = 'father';
						$parent = $parent_mapping{$local_sample}{father}
					}
					
					
					my $parent_first_allele = my $parent_second_allele;
					my $pz = $self->{var_data}{$var_key}{zyg}{$parent};

					if ($self->{allele_data} eq 'bam') {
						($parent_first_allele,$parent_second_allele) = $self->_get_alleles_from_pileup($pz,$self->{var_data}{$var_key}{pileup_str}{$parent});
					} else {
						($parent_first_allele,$parent_second_allele) = $self->_get_alleles_from_varkey($cz,$pz,$var_key);
					}
					my $inherited_allele = '?';
						
					#Data for mother available; can still deduce some cases
					if (($pz eq 'ref' && $cz eq 'ref') || ($pz eq 'het' && $cz eq 'ref') || ($pz eq 'ref' && $cz eq 'het')) {
						$inherited_allele = $parent_first_allele; #inherits ref allele
					} elsif (($pz eq 'het' && $cz eq 'hom') || ($pz eq 'hom' && $cz eq 'het') || ($pz eq 'hom' && $cz eq 'hom')) {
						$inherited_allele = $parent_second_allele; #inherits variant allele
					}
					
					
					if ($cz eq 'het') {
						#Flag if it's definite compound het 
						if ($pz eq 'hom') {
							#single parent variant allele						
							$self->{gene_data}{$gene}{comhet}{$local_sample}{$parent_name}{$het_key}++;
							if ($self->_affected($local_sample)) {
								$affected_allele_data{$parent_name}++;
							} else {
								$affected_alleles_from_one_parent = 'No';
							}
						} 
						
						#Here the mutant allele must have come from the other parent
						if ($pz eq 'ref') {
							my $other_parent = $parent_name eq 'mother'?'father':'mother';
							$self->{gene_data}{$gene}{comhet}{$local_sample}{$other_parent}{$het_key}++;
							if ($self->_affected($local_sample)) {
								$affected_allele_data{$other_parent}++;
							} else {
								$affected_alleles_from_one_parent = 'No';
							}
						}
						
						if ($pz eq 'het') {
							#ambiguous case
							$self->{gene_data}{$gene}{comhet}{$local_sample}{either}{$het_key}++;
						}
						
						
					}
					
					if ($self->_unaffected($local_sample) && $cz ne 'ref' && $cz ne '?') {
						$affected_alleles_from_one_parent = 'No';
					}		
					
					$self->{var_data}{$var_key}{alleles}{$local_sample}{$parent_name} = $inherited_allele;
					
								
				} 
			}
			#delete het entries if any pre-conditions are broken
			if ($self->_unaffected($local_sample)) {
				#unaffected can't be homozygous
				if ($cz eq 'hom') {
					$del_com_het = 1;
				}
			} 
			if ($self->_affected($local_sample)) {
				#all affected must be hets
				if ($cz ne 'het') {
					$del_com_het = 1;
				}
			}
		}
		
		#If we broke comhet conditions for variant we need to remove existing entries
		if (exists $self->{gene_data}{$gene}{comhet}) {
			if ($del_com_het) {
				for my $sample (keys %{$self->{gene_data}{$gene}{comhet}}) {
			
					for my $het_parent (keys %{$self->{gene_data}{$gene}{comhet}{$sample}}) {
						if (exists $self->{gene_data}{$gene}{comhet}{$sample}{$het_parent}{$het_key}) {
							if (keys %{$self->{gene_data}{$gene}{comhet}{$sample}{$het_parent}} > 1) {
								#If other entries exist only remove coordinate
								delete $self->{gene_data}{$gene}{comhet}{$sample}{$het_parent}{$het_key};																								
							} else {
								#If last entry remove the het_parent field as well
								delete $self->{gene_data}{$gene}{comhet}{$sample}{$het_parent};
							}
						}
					}
					if (keys %{$self->{gene_data}{$gene}{comhet}{$sample}} == 0) {
						#Nothing left for sample entry
						delete $self->{gene_data}{$gene}{comhet}{$sample};
						if (keys %{$self->{gene_data}{$gene}{comhet}} == 0) {
							#If nothing left for gene remove entire entry
							delete $self->{gene_data}{$gene}{comhet};
						}
					}
				} 
			}
		}
		
		#Check conditions weren't broken and that one parent only gave alleles
		if ($affected_alleles_from_one_parent eq 'yes' && keys %affected_allele_data == 1) {
			if (exists $affected_allele_data{mother}) {
				$self->{var_data}{$var_key}{parent_allele_affected} = "Mother gives affected allele";
			} elsif (exists $affected_allele_data{father}) {
				$self->{var_data}{$var_key}{parent_allele_affected} = "Father gives affected allele";
			} else {
				print Dumper \%affected_allele_data;
				modules::Exception->throw("ERROR: Incorrect key for affected single allele");
			}
		} elsif ($affected_alleles_from_one_parent eq '?') {
			$self->{var_data}{$var_key}{parent_allele_affected} = '?';
		} else {
			$self->{var_data}{$var_key}{parent_allele_affected} = 'No';
		}
	
		#Get mendelian inheritance info
		my $mendel_inheritance = 'mendelian_rules_followed';
		if ($self->{parent}) {
			#Check if it's mendelian
			my $mendel_str = $self->_mendel_inheritance($var_key);
			if ($mendel_str ne 'yes') {
				$mendel_inheritance = $mendel_str;
			} 
		} else {
			$mendel_inheritance = 'No Parental info';
		}
	
		$self->{var_data}{$var_key}{mendel} = $mendel_inheritance;
		
		#Get disease inheritance info
		my $disease_inheritance = $self->_disease_inheritance($var_key);
		$self->{var_data}{$var_key}{dis_inh} = $disease_inheritance;
	
		
	
		#Generate parent allele string
		if ($self->{parent}) {
			my %parent_data = ();
			my $mother_allele_string = my $father_allele_string;
			for my $sample (keys %{$self->{var_data}{$var_key}{alleles}}) {
				if (exists $parent_mapping{$sample}{mother}) {
					$mother_allele_string .= $sample . ' ('.$self->{var_data}{$var_key}{alleles}{$sample}{mother}.'),';
				}
				if (exists $parent_mapping{$sample}{father}) {
					$father_allele_string .= $sample . ' ('.$self->{var_data}{$var_key}{alleles}{$sample}{father}.'),';
				}
			}
			
			if ($mother_allele_string && $mother_allele_string =~ /\w/) {
				$mother_allele_string =~ s/,$//;
				$self->{var_data}{$var_key}{mother_allele} = $mother_allele_string;
			}
			if ($father_allele_string && $father_allele_string =~ /\w/) {
				$father_allele_string =~ s/,$//;
				$self->{var_data}{$var_key}{father_allele} = $father_allele_string;
			}
		}	
	}	
}


#Get the compound het info and update disease_inheritance field
sub get_compound_het {
	my ($self) = @_;
	
	
	return unless $self->{parent};
	my %compound_hets = ();
	
	#First figure out if gene is comhet gene
	for my $gene (keys %{$self->{gene_data}}) {
		if (exists $self->{gene_data}{$gene}{comhet}) {
			#Candidate gene here
			for my $sample (keys %{$self->{gene_data}{$gene}{comhet}}) {
				if (exists $self->{gene_data}{$gene}{comhet}{$sample}{mother} && exists $self->{gene_data}{$gene}{comhet}{$sample}{father}) {
					#def com het
					my $def_compound_het = '(M='.join(',',keys %{$self->{gene_data}{$gene}{comhet}{$sample}{mother}}).' : F='.join(',',keys %{$self->{gene_data}{$gene}{comhet}{$sample}{father}}).')';
					push @{$compound_hets{$gene}{def}{$def_compound_het}}, $sample;
				} elsif (exists $self->{gene_data}{$gene}{comhet}{$sample}{mother} && exists $self->{gene_data}{$gene}{comhet}{$sample}{either}) {
					my $pos_compound_het = '(M='.join(',',keys %{$self->{gene_data}{$gene}{comhet}{$sample}{mother}}).' : AMB='.join(',',keys %{$self->{gene_data}{$gene}{comhet}{$sample}{either}}).')';
					push @{$compound_hets{$gene}{pos}{$pos_compound_het}}, $sample;
				} elsif (exists $self->{gene_data}{$gene}{comhet}{$sample}{father} && exists $self->{gene_data}{$gene}{comhet}{$sample}{either}) {
					my $pos_compound_het = '(F='.join(',',keys %{$self->{gene_data}{$gene}{comhet}{$sample}{father}}).' : AMB='.join(',',keys %{$self->{gene_data}{$gene}{comhet}{$sample}{either}}).')';
					push @{$compound_hets{$gene}{pos}{$pos_compound_het}}, $sample;
				} elsif (exists $self->{gene_data}{$gene}{comhet}{$sample}{either} && keys %{$self->{gene_data}{$gene}{comhet}{$sample}{either}} > 1) {
					my $pos_compound_het = '(AMB='.join(',',keys %{$self->{gene_data}{$gene}{comhet}{$sample}{either}}).')';
					push @{$compound_hets{$gene}{pos}{$pos_compound_het}}, $sample;
				}
			}
		}
	}
	
	
	#Then update the var_data fields with this info
	for my $var_key (keys $self->{var_data}) {				
		my $def_com_het = 'No';	
		my $pos_com_het = 'No';
		my $inheritance = $self->{var_data}{$var_key}{dis_inh};
		my $gene = $self->{var_data}{$var_key}{gene};
		if (exists $compound_hets{$gene}) {
			if (exists $compound_hets{$gene}{def}) {
				my $def_ok = 1;
				#Check that unaffected and affected don't have same alleles inherited
				for my $def_str (keys %{$compound_hets{$gene}{def}}) {
					#Confirm that this variant is included
					if ($def_str !~ /$var_key/) {
						$def_ok = 0
					}
					if (@{$compound_hets{$gene}{def}{$def_str}} > 1) {
						my $aff = my $unaff = 0;
						for my $sample (@{$compound_hets{$gene}{def}{$def_str}}) {
							if ($self->_affected($sample)) {
								$aff = 1;
							} else {
								$unaff = 1;
							}
						}
						if ($aff && $unaff) {
							$def_ok = 0;
						}
					}
				}
				if ($def_ok) {
					for my $def_str (keys %{$compound_hets{$gene}{def}}) {
						$def_com_het .= join(",",@{$compound_hets{$gene}{def}{$def_str}}) . $def_str . ',';
					}
					$def_com_het =~ s/^No//;	
					$def_com_het =~ s/,$//;
					if ($inheritance eq 'none') {
						$inheritance = "def_compound_het";
					} else {
						$inheritance .= ",def_compound_het";
					}
				} else {
					$def_com_het = 'No';
				}
			}
			
			if (exists $compound_hets{$gene}{pos}) {
				my $pos_ok = 1;
				#Check that unaffected and affected don't have same alleles inherited
				for my $pos_str (keys %{$compound_hets{$gene}{pos}}) {
					#Confirm that this variant is included
					if ($pos_str !~ /$var_key/) {
						$pos_ok = 0
					}
					if (@{$compound_hets{$gene}{pos}{$pos_str}} > 1) {
						my $aff = my $unaff = 0;
						for my $sample (@{$compound_hets{$gene}{pos}{$pos_str}}) {
							if ($self->_affected($sample)) {
								$aff = 1;
							} else {
								$unaff = 1;
							}
						}
						if ($aff && $unaff) {
							$pos_ok = 0;
						}
					}
				}
				if ($pos_ok) {
					for my $pos_str (keys %{$compound_hets{$gene}{pos}}) {
						$pos_com_het .= join(",",@{$compound_hets{$gene}{pos}{$pos_str}}) . $pos_str . ',';
					}
						
					$pos_com_het =~ s/,$//;
					$pos_com_het =~ s/^No//;
					if ($inheritance eq 'none') {
						$inheritance = "pos_compound_het";
					} else {
						$inheritance .= ",pos_compound_het";
					}
				} else {
					$pos_com_het = 'No';
				}
			}
		}
		
		
		#Finally update fields
		$self->{var_data}{$var_key}{dis_inh} = $inheritance;
		$self->{var_data}{$var_key}{def_com_het} = $def_com_het;
		$self->{var_data}{$var_key}{pos_com_het} = $pos_com_het;
	}
	
   	
}






#Get phasing blocks
sub get_phasing_blocks {
	my ($self) = @_;
	
	return unless $self->{parent};

	#flags for detecting blocks
	my $mother_block_start = 0;
	my $father_block_start = 0;
	my $mother_block_end = 0;
	my $father_block_end = 0;
	my $current_var_count = 0;
	my %block_coords = ();

	#params for comparison
	my $min_variant_num = 2;
	my $min_block = 10000;
	
	#ok to use these filters here as they don't remove results, they change the definition of block
	if ($self->{filter_output}) {
		if ($self->{filters}{phase_block_size}) {
			#overwrite if necessary
			$min_block = $self->{filters}{phase_block_size};			
		}
		if ($self->{filters}{phase_var_num}) {
			#overwrite if necessary
			$min_variant_num = $self->{filters}{phase_var_num};			
		}
	} 

	my $prev_chr = 'N';
	for my $var_key (sort {my ($a_chr,$a_coord) = $a =~ /([0-9XYMT]+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /([0-9XYMT]+):(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys %{$self->{var_data}}) {
		$self->{var_data}{$var_key}{phase_data} = 'No'; #Set the default
		my ($chr,$coord) = $var_key =~ /^([0-9XYMT]+):(\d+)/;
		if ($prev_chr ne $chr) {
			if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
				if ($mother_block_start != 0) {
					my $block_count = keys %block_coords;
					my $size = $mother_block_end - $mother_block_start;
					for my $block_key (keys %block_coords) {
						$self->{var_data}{$block_key}{phase_data} = $block_count . ' variants; Mother; ' . $chr .':'.$mother_block_start.'-'.$mother_block_end.'; '.$size .'bp';
					}
				} elsif ($father_block_start != 0) {
					my $block_count = keys %block_coords;
					my $size = $father_block_end - $father_block_start;
					for my $block_key (keys %block_coords) {
						$self->{var_data}{$block_key}{phase_data} = $block_count . ' variants; Father; ' . $chr .':'.$father_block_start.'-'.$father_block_end.'; '.$size .'bp';
					}							
				} else {
					modules::Exception->throw("ERROR: Can't find block");
				}
			}
			$prev_chr = $chr;
			$current_var_count = 0;
			$father_block_start = 0;
			$father_block_end = 0;
			$mother_block_start = 0;
			$mother_block_end = 0;
			%block_coords = ();
			
		}
		
		if ($self->{var_data}{$var_key}{parent_allele_affected} =~ /Mother/) {
			#Mother allele
			if ($father_block_start != 0 ) {
				#Check if we're coming out of father block
				if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
					my $block_count = keys %block_coords;
					my $size = $father_block_end - $father_block_start;
					for my $block_key (keys %block_coords) {
						$self->{var_data}{$block_key}{phase_data} = $block_count . ' variants; Father; ' . $chr .':'.$father_block_start.'-'.$father_block_end.'; '.$size .'bp';
					}	
				}
				
				$current_var_count = 0;
				$father_block_start = 0;
				$father_block_end = 0;
				%block_coords = ();
			}
			
			if ($mother_block_start == 0) {
				$mother_block_start = $coord;
			}
			$mother_block_end = $coord;
			$block_coords{$var_key}++;
			$current_var_count++;
			
		} elsif ($self->{var_data}{$var_key}{parent_allele_affected} =~ /Father/) {
			#Father_allele found
			if ($mother_block_start != 0) {
				#Check if we're coming out of mother block
				if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
					my $block_count = keys %block_coords;
					my $size = $mother_block_end - $mother_block_start;
					for my $block_key (keys %block_coords) {
						$self->{var_data}{$block_key}{phase_data} = $block_count . ' variants; Mother; ' . $chr .':'.$mother_block_start.'-'.$mother_block_end.'; '.$size .'bp';
					}
				}
				
				$current_var_count = 0;
				$mother_block_start = 0;
				$mother_block_end = 0;
				%block_coords = ();
			}

			if ($father_block_start == 0) {
				$father_block_start = $coord;
			}			
			$father_block_end = $coord;
			$block_coords{$var_key}++;
			$current_var_count++;
			
		} elsif ($self->{var_data}{$var_key}{parent_allele_affected} eq 'No') {
			#Trigger to check for end of a block
			if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
				if ($mother_block_start != 0) {
					my $block_count = keys %block_coords;
					my $size = $mother_block_end - $mother_block_start;
					for my $block_key (keys %block_coords) {
						$self->{var_data}{$block_key}{phase_data} = $block_count . ' variants; Mother; ' . $chr .':'.$mother_block_start.'-'.$mother_block_end.'; '.$size .'bp';
					}
				} elsif ($father_block_start != 0) {
					my $block_count = keys %block_coords;
					my $size = $father_block_end - $father_block_start;
					for my $block_key (keys %block_coords) {
						$self->{var_data}{$block_key}{phase_data} = $block_count . ' variants; Father; ' . $chr .':'.$father_block_start.'-'.$father_block_end.'; '.$size .'bp';
					}	
				} else {
					modules::Exception->throw("ERROR: Can't find block");
				}
			}					
			%block_coords = ();
			$current_var_count = 0;
			$father_block_start = 0;
			$father_block_end = 0;
			$mother_block_start = 0;
			$mother_block_end = 0;
		}
	}

	#Check final case if it's at end of last chromosome	
	if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
		if ($mother_block_start != 0) {
			my $block_count = keys %block_coords;
			my $size = $mother_block_end - $mother_block_start;
			for my $block_key (keys %block_coords) {
				$self->{var_data}{$block_key}{phase_data} = $block_count . ' variants; Mother; ' . $prev_chr .':'.$mother_block_start.'-'.$mother_block_end.'; '.$size .'bp';
			}
		} elsif ($father_block_start != 0) {
			my $block_count = keys %block_coords;
			my $size = $father_block_end - $father_block_start;
			for my $block_key (keys %block_coords) {
				$self->{var_data}{$block_key}{phase_data} = $block_count . ' variants; Father; ' . $prev_chr .':'.$father_block_start.'-'.$father_block_end.'; '.$size .'bp';
			}	
		} else {
			modules::Exception->throw("ERROR: Can't find block");
		}			
	}
	
	#print Dumper $self;
	#exit;					
}

#Filter out results before reporting
sub filter_results {
	my ($self) = @_;
	
	if (!$self->{filter_output}) {
		return;
	} else {
		if ($self->{debug}) {
			print "Apply the following filters....\n";
			for my $filter (keys %{$self->{filters}}) {
				if ($filter eq 'gene_list') {
					print "Filter on gene_list\n";				
				} else {
					print "Filter $filter $self->{filters}{$filter}\n";
					
				}
			}
		}
		
		for my $var_key (keys %{$self->{var_data}}) {
			my ($chr,$start,$end) = $var_key =~ /([0-9XYMT]+):(\d+)\-(\d+)/;
			
			#Filter by chr; should be handled in parse_vep and parse_vcf
			if (exists $self->{filters}{chrom} && $self->{filters}{chrom} ne $chr) {
				delete $self->{var_data}{$var_key};
				next;
			}
			
			#Filter by genomic coordinate
			if (exists $self->{filters}{start} && exists $self->{filters}{end} && ($self->{filters}{start} > $end || $self->{filters}{end} < $start)) {
				delete $self->{var_data}{$var_key};
				next;
			}
			
			#Check the min number of affected samples are variant for the base
			if (exists $self->{filters}{min_num_aff} && $self->{var_data}{$var_key}{aff_count} < $self->{filters}{min_num_aff}) {
				delete $self->{var_data}{$var_key};
				next;
			}
			
			#Fail allele freq if has freq and too high;  leave novel or unknown cases in
			if (exists $self->{filters}{max_allele_freq} && $self->{var_data}{$var_key}{allele_freq} ne 'N/A' && $self->{filters}{max_allele_freq} <= $self->{var_data}{$var_key}{allele_freq}) {
				delete $self->{var_data}{$var_key};
				next;
			}
			
			#Filter out cases without polyphen when polyphen filter applied
			if (exists $self->{filters}{polyphen} && !exists $self->{var_data}{$var_key}{poly_pred}) {
				delete $self->{var_data}{$var_key};
				next;
			} elsif (exists $self->{filters}{polyphen} && $self->{filters}{polyphen} ne $self->{var_data}{$var_key}{poly_pred}) {
				delete $self->{var_data}{$var_key};
				next;
			}
			
			#Same for SIFT
			if (exists $self->{filters}{sift} && !exists $self->{var_data}{$var_key}{sift_pred}) {
				delete $self->{var_data}{$var_key};
				next;
			} elsif (exists $self->{filters}{sift} && $self->{filters}{sift} ne $self->{var_data}{$var_key}{sift_pred}) {
				delete $self->{var_data}{$var_key};
				next;
			}
			
			#Filter by gene name; check for ensembl gene and transcripts as well as hgnc
			if (exists $self->{filters}{gene_list}) {
				my $gene_found = 0;
				if (!exists  $self->{var_data}{$var_key} ) {
					
				}
				my $ens_gene = $self->{var_data}{$var_key}{gene};
				my $ens_trans = $self->{var_data}{$var_key}{transcript};
				
				for my $include_gene ( keys %{$self->{filters}{gene_list}} ) {
					#Check ensembl gene, transcripts, and hgnc to allow filtering on all three
				    if ($ens_gene eq $include_gene) {
				    	$gene_found = 1;
				    }
				    if ($ens_trans eq $include_gene) {
				    	$gene_found = 1;
				    }
				    
				    #This one's optional so check exists first 
				    if (exists $self->{gene_data}{$ens_gene}{hgnc} && $self->{gene_data}{$ens_gene}{hgnc} eq $include_gene) {
				    	$gene_found = 1;
				    }
				}
				if (!$gene_found) {
					delete $self->{var_data}{$var_key};
					next;
			   	}
			}
			
			#Filter on inheritance; use regex as can contain comhet string as well
			if ($self->{filters}{inheritance} && $self->{var_data}{$var_key}{dis_inh} !~  /$self->{filters}{inheritance}/) {
				delete $self->{var_data}{$var_key};
				next;
			}
			
			#Filter for denovo
			if ($self->{filters}{denovo} && $self->{var_data}{$var_key}{mendel} !~ /de novo/) {
				delete $self->{var_data}{$var_key};
				next;
			}
			
			#Filter on comhet; allow definite or possible cases through
			if ($self->{filters}{comhet} &&  $self->{var_data}{$var_key}{def_com_het} eq 'No' && $self->{var_data}{$var_key}{pos_com_het} eq 'No') {
				delete $self->{var_data}{$var_key};
				next;
			}
		
				
		}
		
		
	}
	
	$self->_check_vars_remain();
	
}




#Generate the full summary lines
sub generate_line_data {
	my ($self) = @_;
	
	$self->_check_vars_remain();
	#Add up total variants for each gene
	my %ens_gene_count = ();
	for my $var_key ( keys %{$self->{var_data}} ) {
	    $ens_gene_count{$self->{var_data}{$var_key}{gene}} +=  $self->{var_data}{$var_key}{aff_count} + $self->{var_data}{$var_key}{unaff_count};
	}
	for my $gene (keys %ens_gene_count) {
		$self->{gene_data}{$gene}{total_count} = $ens_gene_count{$gene};
	}
	
	#print Dumper $self;
	
	my @headers = qw(variant_samples);
	
	if ($self->{total_affected}) {
		push @headers, 'number_affected_variant/total_affected';
	}
	if ($self->{total_unaffected}) {
		push @headers, 'number_unaffected_variant/total_unaffected';
	}
	
	my $mother = my $father = 0;

	if ($self->{parent}) {
		for my $sample (keys %{$self->{parent_mapping}}) {
			if (exists $self->{parent_mapping}{$sample}{mother}) {
				$mother = 1;
			}
			if (exists $self->{parent_mapping}{$sample}{father}) {
				$father = 1;
			}
		}
		
		push @headers, 'mendelian_inheritance','disease_inheritance','gene','total_gene_variants','unique_coord_gene_variants';

		if ($self->{allele_data} eq 'bam') {
			for my $sample (sort keys %{$self->{sample_data}}) {
				push @headers, $sample.'_pileup(zyg)';
			}
		} else {
			for my $sample (sort keys %{$self->{sample_data}}) {
				push @headers, $sample.'_zyg';
			}
		}	
		
		
		if ($mother && $father) {
			push @headers,'mother_allele','father_allele','parent_allele_common_to_affected','affected_allele_block','definite_compound_het','possible_compound_het';
		} elsif ($mother) {
			push @headers,'mother_allele','parent_allele_common_to_affected','affected_allele_block','definite_compound_het','possible_compound_het';
		} elsif ($father) {
			push @headers,'father_allele','parent_allele_common_to_affected','affected_allele_block','definite_compound_het','possible_compound_het';		
		}
	} else {
		push @headers, 'disease_inheritance','gene','total_gene_variants','unique_coord_gene_variants';
		if ($self->{allele_data} eq 'bam') {
			for my $sample (sort keys %{$self->{sample_data}}) {
				push @headers, $sample.'_pileup(zyg)';
			}
		} else {
			for my $sample (sort keys %{$self->{sample_data}}) {
				push @headers, $sample.'_zyg';
			}
		}	
	}	


	push @headers, 'chr','start_coord','end_coord','var_type','variant','var_quality','allele_freq','dbsnp','ens_transcript','vep_classifier','aa_change','polyphen','sift';
			
	$self->{headers} = \@headers;
	#print Dumper \@headers;
	for my $var_key (keys $self->{var_data}) {
		my @line_data = (); #build up the lines
		my %variant_samples = ();
		my $aff_count = 0;
		my $unaff_count = 0;
		
		for my $sample (keys %{$self->{var_data}{$var_key}{zyg}}) {
			if ($self->{var_data}{$var_key}{zyg}{$sample} eq 'hom' || $self->{var_data}{$var_key}{zyg}{$sample} eq 'het') {
				$variant_samples{$sample}++;
				if ($self->_affected($sample)) {
					$aff_count++;
				} else {
					$unaff_count++;
				}
			}
		}
		
		push @line_data, join(",",sort keys %variant_samples);
		if ($self->{total_affected}) {
			push @line_data, $aff_count . ' of '.$self->{total_affected};
		}
		if ($self->{total_unaffected}) {
			push @line_data, $unaff_count . ' of '.$self->{total_unaffected};
		}
		
		if ($self->{parent}) {
			push @line_data, $self->{var_data}{$var_key}{mendel};
		}
		my $gene = exists $self->{var_data}{$var_key}{hgnc}?$self->{var_data}{$var_key}{hgnc}:$self->{var_data}{$var_key}{gene};
		
		my $ens_gene = $self->{var_data}{$var_key}{gene};
		my $uniq_count = 0;
		
		#Don't count entries that have been removed during vcf parsing; can't do it there because of overlapping canonical transcript cases
		for my $uniq_coord (keys %{$self->{gene_data}{$ens_gene}{uniq_count}}) {
			if (exists $self->{var_data}{$uniq_coord} && $self->{var_data}{$uniq_coord}{gene} eq $ens_gene) {
				$uniq_count++;
			}
		}
		my $total_count = $self->{gene_data}{$ens_gene}{total_count};
		
		
		push @line_data, $self->{var_data}{$var_key}{dis_inh}, $gene, $total_count,$uniq_count;
		
		
		if ($self->{allele_data} eq 'bam') {
			for my $sample (sort keys %{$self->{sample_data}}) {
				push @line_data, $self->{var_data}{$var_key}{pileup_str}{$sample} .'('.$self->{var_data}{$var_key}{zyg}{$sample}.')';			
			}
		} else {
			for my $sample (sort keys %{$self->{sample_data}}) {
				push @line_data, $self->{var_data}{$var_key}{zyg}{$sample};			
			}
		}
		
		
		if ($mother) {
			push @line_data, $self->{var_data}{$var_key}{mother_allele};
		}
		if ($father) {
			push @line_data, $self->{var_data}{$var_key}{father_allele};
		}
		if ($self->{parent}) {
			push @line_data, $self->{var_data}{$var_key}{parent_allele_affected}, $self->{var_data}{$var_key}{phase_data}, $self->{var_data}{$var_key}{def_com_het}, $self->{var_data}{$var_key}{pos_com_het};
		}
		
		
		
		
		my ($chr,$range,$event) = split(':',$var_key);
		my ($start,$end) = split('-',$range);
		
		my $allele_freq = $self->{var_data}{$var_key}{allele_freq};
		my $dbsnp = exists $self->{var_data}{$var_key}{dbsnp}?$self->{var_data}{$var_key}{dbsnp}:'Novel';

		push @line_data, $chr,$start,$end,$self->{var_data}{$var_key}{var_type},$event,$self->{var_data}{$var_key}{qual}, $allele_freq,$dbsnp,$self->{var_data}{$var_key}{transcript},$self->{var_data}{$var_key}{vep};	
		
		if (exists $self->{var_data}{$var_key}{aa_change}) {
			push @line_data, $self->{var_data}{$var_key}{aa_change};
		} else {
			push @line_data, "N/A";
		}
		
		if (exists $self->{var_data}{$var_key}{poly_pred}) {
			push @line_data, $self->{var_data}{$var_key}{poly_pred} .':'.$self->{var_data}{$var_key}{poly_score};
		} else {
			push @line_data, "N/A";
		}
		
		if (exists $self->{var_data}{$var_key}{sift_pred}) {
			push @line_data, $self->{var_data}{$var_key}{sift_pred} .':'.$self->{var_data}{$var_key}{sift_score};
		} else {
			push @line_data, "N/A";
		}
		
		my $pass = 0;
		
		if ($allele_freq eq 'N/A' || $allele_freq <= 0.02) {
			$pass = 1;
		}
		if ($allele_freq eq 'N/A') {
			$allele_freq = -1; #For sorting
		}
		
	
		
		my $var_count = keys $self->{var_data};
	   	
	   	if ($var_count == 0) {
	   		modules::Exception->warning("All variants are filtered out; please rerun with different filters");
	   		exit;
	   	}
		
		
		#print Dumper \@line_data;
		#print Dumper $self->{var_data}{$var_key};
		my $line = join("\t",@line_data);
		$self->{lines}{$pass}{$aff_count}{$unaff_count}{$allele_freq}{$chr}{$start} = $line;
	}		
}



	
#sort lines and write to file
sub write_to_files {
	my ($self) = @_;
		
	my $out = $self->{out};	
		
	open(OUT,">$out") || modules::Exception->throw("Can't open file to write $out\n");
	my $header_line = join("\t",@{$self->{headers}});
	print OUT $header_line . "\n\n";
	
	
	#Sort by non-mendelian first; then highest likely affected variant and lowest unaffected non-variant 
	for my $pass (sort {$b<=>$a} keys %{$self->{lines}}) {
		if ($pass == 1) {
			print OUT "\nNOVEL/RARE VARIANTS\n\n";
		} else {
			print OUT "\n\nCOMMON VARIANTS\n\n";
		}
		for my $aff_count (sort {$b<=>$a} keys %{$self->{lines}{$pass}}) {
			for my $unaff_count (sort {$a<=>$b} keys %{$self->{lines}{$pass}{$aff_count}}) {
				for my $dbsnp_freq (sort {$a<=>$b} keys %{$self->{lines}{$pass}{$aff_count}{$unaff_count}}) {
					for my $chr (sort {$a cmp $b} keys %{$self->{lines}{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}}) {
						for my $coord (sort {$a<=>$b} keys %{$self->{lines}{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}{$chr}}) {
							#print "Pass $pass Aff $aff_count Unaff $unaff_count dbsnp $dbsnp_freq Chr $chr Coord $coord\n";
							print OUT $self->{lines}{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}{$chr}{$coord}. "\n";
						}
					}
				}
			}
		}
	}
	
	close OUT;
	
}
	

	
#Check whether conditions qualify as a block
sub _check_block {
	my ($self,$mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block) = @_;
	
	if ($current_var_count >= $min_variant_num) {
		#print "if $current_var_count >= $min_variant_num && $mother_block_start != 0 && ($mother_block_end-$mother_block_start) > $min_block\n";
		if ($mother_block_start != 0 && ($mother_block_end-$mother_block_start) > $min_block) {
			return 1;
		}
		if ($father_block_start != 0 && ($father_block_end-$father_block_start) > $min_block) {
			return 1;
		}
	} 
		
	return 0;
	
}	




#Checks whether Mendelian inheritance is upheld with one or more parent sequenced
sub _mendel_inheritance {
	my ($self,$var_key) = @_;
	
	my %parent_mapping = %{$self->{parent_mapping}};
	
	for my $sample (keys %parent_mapping) {
		for my $parent (keys %{$parent_mapping{$sample}}) {
			my $sample_zyg = $self->{var_data}{$var_key}{zyg}{$sample};
			my $parent_zyg = $self->{var_data}{$var_key}{zyg}{$parent_mapping{$sample}{$parent}};
			if ($sample_zyg eq 'No data' || $parent_zyg eq 'No data') {
				return 'missing allele info';
			}
			
			if ($sample_zyg eq 'ref' && $parent_zyg eq 'hom') {
					return "hom_mother ref_child";
			} elsif ($sample_zyg eq 'het') {
				if (keys %{$parent_mapping{$sample}} == 2) {
					my $mother_zyg = $self->{var_data}{$var_key}{zyg}{$parent_mapping{$sample}{mother}};
					my $father_zyg = $self->{var_data}{$var_key}{zyg}{$parent_mapping{$sample}{father}};
					if ($mother_zyg eq 'ref' && $father_zyg eq 'ref') { 
						return "de novo (ref_parents het_child)"; #Special de novo case
					} 
					if ($mother_zyg eq 'hom' && $father_zyg eq 'hom') {
						return "hom_parents het_child";
					}
				} 
			} else {
				if ($sample_zyg eq 'ref' && $parent_zyg eq 'hom') {
					return 'hom_parent ref_child';
				} elsif ($sample_zyg eq 'hom' && $parent_zyg eq 'ref') {
					return 'ref_parent hom_child';
				} elsif ($sample_zyg eq 'hom' && $parent_zyg eq 'ref') {
					return "ref_mother hom_child";
				} 
			}
		}
	}

	return 'yes';
}						

#Get the disease inheritance patterns
sub _disease_inheritance {
	my ($self,$var_key) = @_;
	
	my ($chr) = split(":",$var_key);
	
	my $x_or_auto = $chr eq 'X'?'x':'auto';
	my $aff_var_count = 0;
	my $recessive = 1;
	my $dominant = 1;
	my $disease_inheritance = 'none';
	
	for my $sample (keys %{$self->{sample_data}}) {
		my $sample_zyg = $self->{var_data}{$var_key}{zyg}{$sample};
		if ($sample_zyg eq 'No data') {
			return 'missing allele info';
		}
		if ($self->_affected($sample)) {
			
			if ($sample_zyg ne 'hom') {
				#Affected must be homs for auto-recessive
				$recessive = 0;
			} 
			if ($sample_zyg ne 'het') {
				#Affected must be het for auto-dominant
				$dominant = 0;
			}
			
			$aff_var_count++ unless $sample_zyg eq 'ref';
		} elsif ($self->_unaffected($sample)) {
			if ($sample_zyg eq 'hom') {
				#Unaffected can't be homs for auto-recessive
				$recessive = 0;
			}
			if ($sample_zyg ne 'ref') {
				#Unaffected can't be variant for auto-dominant
				$dominant = 0;
			}
			
		}
		
	}
	
	if ($self->{total_affected} == $aff_var_count && $recessive) {
		$disease_inheritance = $x_or_auto . '-recessive';
	}
	
	if ($self->{total_affected} == $aff_var_count && $dominant) {
		$disease_inheritance = $x_or_auto . '-dominant';
	}
	
	return $disease_inheritance;
}


#Gets variant key from either vcf or vep; standardises naming for loading into data structure
sub _get_variant_key {
	 my @args = @_;
	
	 my %args = @args;


    my @required_args = (
    					-chrom,
    					-first,
    					-ref_seq,
    					-var_seq,
    					-type
    					);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set $args{$required_arg}");
		} 
    }
    
    my $ref = $args{-ref_seq};
    my $var = $args{-var_seq};
    my $first_coord = $args{-first};
    my $chr = $args{-chrom};
    my $type = $args{-type};
    
    my $start_coord = my $end_coord = my $bases;
    my $length_ref = length($ref);
    my $length_var = length($var);
    my $var_type;
    
    if ($type eq 'vcf') {
		if ($length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord + 1;
			$end_coord = $start_coord + $length_ref - $length_var - 1;				
			my $del_length = $length_ref - $length_var;
			#print "VCF R $ref L $del_length\n";
			
			$bases = '-'. substr($ref,1,$del_length);
		} elsif ($length_ref < $length_var) {
			#Add the ref length and var length difference to the coord 
			#$start_coord = $end_coord = $first_coord + 1;
			$var_type = 'INS';
			$start_coord = $end_coord = $first_coord;
			my $ins_length = $length_var - $length_ref;
			$bases = '+'.substr($var,1,$ins_length);
		} else {
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			$bases = $ref . '->' .$var;
		}
    	
    } else {
    	if ($ref eq '-' || $length_ref < $length_var) {
	    	$start_coord = $end_coord = $first_coord - 1;
			$var_type = 'INS';	
			my $ins_length = $length_var - $length_ref;
			$ins_length++ if $ref eq '-';
			$bases = '+'.substr($var,0,$ins_length);
		}  elsif ($var eq '-' || $length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord;
			my $del_length = $length_ref - $length_var;
			$end_coord = $start_coord + $length_ref - $length_var - 1;
			$del_length++ if $var eq '-';
			$end_coord++ if $var eq '-';
			$bases = '-'. substr($ref,0,$del_length);	
		} elsif ($length_ref == $length_var && $length_ref == 1) {
			#single snvs
			$var_type = 'SNV'; 
			$bases = $ref .'->'.$var;
			$start_coord = $end_coord = $first_coord;
		} else {
			modules::Exception->warning("ERROR: Can't identify var_type doesn't match any var type\n");
			#next;
		}
    }
    
	my $var_key = $chr . ':'.$start_coord .'-'.$end_coord .':' .$bases;
	return($var_key,$var_type);
	
}

#Get alleles from variant key and zygosity (when no bam available)
sub _get_alleles_from_varkey {
	my ($self,$child_zyg,$parent_zyg,$var_key) = @_;
	if ($child_zyg eq 'No data' || $parent_zyg eq 'No data') {
		return ('?','?');
	} elsif ($parent_zyg eq 'ref') {
		return ('ref','ref');
	} elsif ($parent_zyg eq 'hom') {
		my $var_allele;
		if ($var_key =~ /->([ATCG])/) {
			$var_allele = $1;
		} else {
			my @fields = split(':',$var_key);
			$var_allele = $fields[-1];
		}
		return ($var_allele,$var_allele);
	} else {
		my $var_allele;
		if ($var_key =~ /->([ATCG])/) {
			$var_allele = $1;
		} else {
			my @fields = split(':',$var_key);
			$var_allele = $fields[-1];
		}
		return ('ref',$var_allele);
	}
}


#Get alleles from pileup string and zygosity (when bam available)
sub _get_alleles_from_pileup {
	my ($self,$zyg,$pileup_str) = @_;
	if ($zyg eq 'No data') {
		return ('?','?');
	} elsif ($zyg eq 'ref') {
		return ('ref','ref');
	} elsif ($zyg eq 'hom') {
		my ($var_allele) = split(':',$pileup_str);
		return ($var_allele,$var_allele);
	} else {
		#het case
		my ($block1,$block2) = split('_',$pileup_str);
		
		if (!defined $block2) {
			modules::Exception->warning("ERROR: Het zyg and pileup str $pileup_str");
			return('?','?');
		}
		
		my ($first_allele) = split(':',$block1);
		my ($second_allele) = split(':',$block2);
		#Return ref as first allele always
		if ($block2 =~ /ref/i){
			return ($second_allele,$first_allele)
		} else {
			return ($first_allele,$second_allele);
		}
	}
}

sub _affected {
	my ($self,$sample) = @_;
	if (exists $self->{sample_data}{$sample}{affected} && $self->{sample_data}{$sample}{affected} == 1) {
		return 1;
	} else {
		return 0;
	}
	
}

sub _unaffected {
	my ($self,$sample) = @_;
	if (exists $self->{sample_data}{$sample}{unaffected} && $self->{sample_data}{$sample}{unaffected} == 1) {
		return 1;
	} else {
		return 0;
	}
	
}

sub _check_vars_remain {
	my ($self) = @_;
	my $var_count = keys %{$self->{var_data}};
	if ($var_count == 0) {
		modules::Exception->throw("ERROR: No variants remain to report; change filtering options");
	} else {
		print "Step run: $var_count variants analysed\n" if $self->{debug};
	}
}

return 1;