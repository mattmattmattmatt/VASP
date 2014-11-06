package modules::Utils;

use strict;
use modules::Exception;
use Data::Dumper;

sub pileup_string {
	my ( $self, $pileup_str ) = @_;
	my @bases = split("",uc($pileup_str));
    my $count = 0;
    my $allele_count = 0;
    my $base_count = @bases;
    my %final_bases = ();
    
    while ($count < $base_count) {
     	if (defined $bases[$count+2] && ($bases[$count+1] eq '+' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
    		my $rest_indel = join("",@bases[$count+2..$#bases]);
			my ($length_indel) = $rest_indel =~ /(\d+)/;
    		my $indel_end = $length_indel+1+length($length_indel)+$count;
    		my $indel_bases = join("",@bases[$count+1..$indel_end]);
    		$count +=  $length_indel+2+length($length_indel);
    		$final_bases{$indel_bases}++;
    		$allele_count++;
    		next;
    	} elsif (defined $bases[$count+2] && ($bases[$count+1] eq '-' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
    		my $rest_indel = join("",@bases[$count+2..$#bases]);
			my ($length_indel) = $rest_indel =~ /(\d+)/;
    		my $indel_end = $length_indel+1+length($length_indel)+$count;
	    	my $indel_bases = join("",@bases[$count+1..$indel_end]);
    		$count +=  $length_indel+2+length($length_indel);
    		$final_bases{$indel_bases}++;
    		$allele_count++;
    		next;
    		
   		} elsif ($bases[$count] =~ /[ATCGN]/) {
	    	$final_bases{$bases[$count]}++;
	    	$allele_count++;
	    } elsif ($bases[$count] eq '.' || $bases[$count] eq ',') {
	    	$final_bases{'ref'}++;
	    	$allele_count++;
	    } elsif ($bases[$count] eq '*') {
	    	$final_bases{'other_del'}++;
	    	$allele_count++;
	    }
	    $count++;
	}
    
    my $zyg = 'het'; #default is het
    my @pileup_strings;
    for my $allele (sort {$final_bases{$b}<=>$final_bases{$a}} keys %final_bases) {
    	#Simply hom test (90% non reference allele) and REF test
    	if ($final_bases{$allele}/$allele_count >= 0.9 && $allele ne 'ref') {
    		$zyg = 'hom';
    	} elsif ($final_bases{$allele}/$allele_count >= 0.9 && $allele eq 'ref') {
    		$zyg = 'ref';
    	}
    	
    	push @pileup_strings, $allele.':'.$final_bases{$allele};
    }
    my $pileup_string = join('_',@pileup_strings);
    
    return ($pileup_string,$zyg);    
	
	
	
}

return 1;
