package  modules::SystemCall;

use strict;


sub new {
    my ($class) = @_;

    my $self = bless {}, $class;

    return $self;
}

sub run {
    my ($self, $command) = @_;
   
	if (defined $command){
		#split up any commands consisting of >1 command so we can treat each one separately
		my @commands = split(';', $command);
		
		#special handling for cd as returns non zero by default but need to run as batch command for archive step
		my $no_check = 0;
		
		if ($command =~ /cd/) {
			my ($dir) = $command =~ /cd\s+(\S+);?/;
			$dir =~ s/;$//;
			if (!-d $dir) {
				die "Directory $dir doesn't exist\n";
			}
			@commands = ();
			push @commands, $command;
		}
		
		
	
		for my $local_command (@commands) {
			if ($local_command =~ /cd/) {
				$local_command .= ';';
			}

			if ($local_command =~ /rm.*null/) {
				#skip calls we're trying to ignore stderr (rm and don't care if file is there or not)
				system($local_command);
				next;
			}
			
			print STDERR "Running command:\n$local_command\n";
			if (my $return = system($local_command)){
			    die "Command exited with non-zero status $return";
			} 
	    }
	    
	}  else {
			warn "No command specified to system call";
			return 0;
	    }

    return 1;
}

return 1;
