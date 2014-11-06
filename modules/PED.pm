package modules::PED;

use strict;
use modules::Exception;
use Data::Dumper;

sub new {
    my ($class) = @_;

    my $self = bless {}, $class;

    return $self;
}

#Validate
sub validate_ped {
	my ($self, @args) = @_;

    my @required_args = (
			             -ped_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-ped_file} ) {
    	modules::Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    open(PED,$args{-ped_file}) || modules::Exception->throw("Can't open file $args{-ped_file}\n");
    
    my $error = 0;
    
    my @ped_lines  = <PED>;
    
    #First get a list of all ped identifiers
    my %ped_ids = ();
    
    for my $line (@ped_lines) {
    	next unless $line =~ /\S/;
    	chomp $line;
    	my @fields = split(' ',$line);
    	if (exists $ped_ids{$fields[1]}) {
    		modules::Exception->throw("ERROR: Multiple entries for $fields[1] in ped file");
    	}
    	$ped_ids{$fields[1]}++;
    }
    
    for my $line (@ped_lines) {
    	next unless $line =~ /\S/;
    	chomp $line;
    	my @fields = split(' ',$line);
    	if (@fields < 6) {
    		$error = "ERROR: Expecting Ped file with at least six fields (family_id, id, father, mother, sex, and affected)";
    	}
    	my ($family,$id,$father,$mother,$sex,$affected) = @fields;
    	
    	if ($affected !~ /[12]/) {
    		$error = "ERROR: Affected value ($affected) must be 1 or 2 for line $line";
    	}
    	if ($sex !~ /[12]/) {
    		$error = "ERROR: Sex value ($sex) must be 1 or 2 for line $line";
    	}
    	if (!exists $ped_ids{$mother} && $mother != 0) {
    		$error = "ERROR: Mother value ($mother) must match previous ped name for line $line";
    	}
    	if (!exists $ped_ids{$father} && $father != 0) {
    		$error = "ERROR: Father value ($father) must match previous ped name for line $line";
    	}
    }
    return $error;
    
}

#Parse a ped file
sub parse_ped {
    my ($self, @args) = @_;

    my @required_args = (
			             -ped_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-ped_file} ) {
    	modules::Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    open(PED,$args{-ped_file}) || modules::Exception->throw("Can't open file $args{-ped_file}\n");
    
    my %ped_data;
    my $count = 1; #Keep order consistent
    
    while (<PED>) {
    	chomp;
    	next unless /\S/; #ignore blank lines
    	my @fields = split;
    	if (@fields < 6) {
    		modules::Exception->throw("ERROR: Expecting Ped file with at least six fields (family_id, id, father, mother, sex, and affected");
    	}
    	my ($family,$id,$father,$mother,$sex,$affected) = @fields;
    	
    	$ped_data{"$count:$id"}{father} = $father unless $father eq '0';
    	$ped_data{"$count:$id"}{mother} = $mother unless $mother eq '0';
    	$ped_data{"$count:$id"}{family} = $family;
    	
    	
    	if ($affected == 2) {
    		$ped_data{"$count:$id"}{affected} = 1;
    	} else {
    		$ped_data{"$count:$id"}{unaffected} = 1;
    	}
    	
    	if ($sex == 1) {
    		$ped_data{"$count:$id"}{sex} = 'male';
    	} else {
    		$ped_data{"$count:$id"}{sex} = 'female';
    	}
    	$count++;
    	
    }
    
    return \%ped_data;
}

1;