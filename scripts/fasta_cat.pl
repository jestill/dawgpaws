#!/usr/bin/perl -w
# STARTED: 08/20/2006
# UPDATED: 04/30/2007

=head1 NAME

FastaCat.pl - Concatenate FASTA files

=head1 SYNOPSIS

   FastaCat.pl -i InputDir -o OutFile.fasta

=head1 DESCRIPTION

Given a directory, concatenate all of the fasta files in
the directory into a single file.

=head1 ARGUMENTS

=over 2

=item -i InDir

Path to the input directory containing the fasta files.

=item -o OutFile.fasta

Path to the fasta file that will be produced.

=back

=head1 AUTHOR

James C. Estill<lt>JamesEstill at gmail.comE<gt>

=cut

use Getopt::Std;

#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
getopts('i:o:', \%Options);

my $InDir = $Options{i};
my $OutFile = $Options{o};

&CatFasta ($InDir, $OutFile);
exit;

sub CatFasta
    
{ #  Begin of CatBlast Subfunction


    my $FileNum = 0;
    my $Directory = $_[0];
    my $FastaFileOut = $_[1];
    # First warn the user if the output file already exists
    my $ArrayRefI = 0;                # Start the arrray ref at 0
    
    # Set the current working directory to the directory that was selected
    chdir $Directory;
    
    # -------------------------------------------
    # Load an array with the files to concateate
    # -------------------------------------------
    # Open the directory containing the BLAST output files
    opendir( DIR, $Directory ) ||
	die "Can not open directory\n$Directory\n";
    # Read the directory ignoring the . and ..
    my @FastaList = grep !/^\.\.?$/, readdir DIR ;
    closedir( DIR );
    
    my $NumFiles = @FastaList;
    open(FASTAOUT, ">> $FastaFileOut");
    
  FILE: foreach (@FastaList)
  { # Begin of foreach loop
      # Temp exit for debug
      #if ($FileNum == 50) { exit; } 
      # Show the name of the file that is being added
      
      # Show the number of the file being processed
      print "Processing ".$FileNum." of ".$NumFiles."\n"; 
      open(FASTAFILE, $_) || ((warn "Can't open file $_\n"), next FILE);
      
      while (<FASTAFILE>) {  # For every line in the FASTA file
	  print FASTAOUT ;   # Print the output to FASTAOUT
      }
      close(FASTAFILE);
  } # End of foreach @BlastList loop
    
    close(FASTAOUT);
    
} # End of CatFasta Subfunction

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 08/20/2006
#  - Modified program from a BLAST concatenation program I 
#    wrote a while ago.
# 04/04/2007
#  - Changed from hard coded variables to command line
# 04/30/2007
#  - Added POD style documentation

#-----------------------------------------------------------+
# TO DO                                                     |
#-----------------------------------------------------------+
#
# 08/20/2006 
#  - Allow this to search for specific extensions
#    or ignore specific extensions.
