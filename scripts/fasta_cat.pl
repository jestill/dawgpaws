#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_cat.pl - Concatenate all fasta files in a directroy |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# STARTED: 08/20/2006                                       |
# UPDATED: 01/30/2008                                       |
#                                                           |
#-----------------------------------------------------------+
# TODO:
# Update to Getopt::Long
use Getopt::Long;
#use Getopt::Std;

#
my $indir;
my $outfile;

# BOOLEANS                                                                                  
my $show_help = 0;             # Show program help                                          
my $show_version = 0;          # Show program version                                       
my $show_man = 0;              # Show program manual page using peldoc                      
my $show_usage = 0;            # Show program usage command                                 
my $quiet = 0;                 # Boolean for reduced output to STOUT                        
my $test = 0;                  # Run the program in test mode                               
my $verbose = 0;               # Run the program in verbose mode 

#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
my $ok = GetOptions(
                    # Required Arguments                                                    
                    "i|indir=s"    => \$indir,
                    "o|outfile=s"  => \$outfile,
                    # Booleans                                                              
                    "verbose"      => \$verbose,
                    "test"         => \$test,
                    "usage"        => \$show_usage,
                    "version"      => \$show_version,
                    "man"          => \$show_man,
                    "h|help"       => \$show_help,
                    "q|quiet"      => \$quiet,);

#-----------------------------+                                                             
# CHECK FOR SLASH IN DIR      |                                                             
# VARIABLES                   |                                                             
#-----------------------------+ 
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

&cat_fasta ($indir, $outfile);

exit;

sub cat_fasta { #  Begin of CatBlast Subfunction
    
    my $indir = $_[0];
    my $fasta_file_path = $_[1];
    
    my $file_num = 0;
    
    # First warn the user if the output file already exists
    my $ArrayRefI = 0;                # Start the arrray ref at 0
    
    # Set the current working directory to the directory that was selected
    #chdir $Directory;
    
    #-----------------------------+                                                             
    # Get the FASTA files from the|                                                             
    # directory provided by the   |                                                             
    # var $indir                  |                                                             
    #-----------------------------+                                                             
    opendir( DIR, $indir ) ||
	die "Can't open directory:\n$indir";
    my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
    closedir( DIR );
    
    my $total_files = @fasta_files;
    
    # Currently this is set to append to an existing file
    #my $NumFiles = @FastaList;
    open(FASTAOUT, ">> $fasta_file_path") ||
	die "Can not open fasta file:\n$fasta_file_path\n";
    
    for my $ind_file (@fasta_files) {
	
	$file_num++;
	my $in_file_path = $indir.$ind_file;

	print "Processing ".$file_num." of ".$total_files."\n"; 
	open(FASTAFILE, $in_file_path) || 
	    die "Can not open input file:\n$in_file_path";
	
	while (<FASTAFILE>) {  # For every line in the FASTA file
	    print FASTAOUT ;   # Print the output to FASTAOUT
	}
	close(FASTAFILE);
	
    } # End of foreach @BlastList loop
    
    close(FASTAOUT);
    
} # End of cat_fasta Subfunction

__END__

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
