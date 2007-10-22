#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           | 
# RunGenscan.pl                                             |
#                                                           |
#-----------------------------------------------------------+
# AUTHORS: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.org                       |
# STARTED: 04/17/2007                                       |
# UPDATED: 04/17/2007                                       |
# DESCRIPTION:                                              |
#  Run the genscan gene prediction program.                 |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                   # Keeps thing running clean
use Getopt::Std;              # Get options from command line

#-----------------------------+
# VARIABLES                   |
#-----------------------------+
my $InFile;                    # Full path to the input file
my $OutGen;                    # Full path to the output file in genscan format
my $OutGff;                    # Full path to the output file in gff format
                               # This will require
my $cmd;                       # The system command that will be run
my $PrintHelp;                 # Boolean to indicate to print the Usage string

# The training file to use for gene prediction
my $Lib = "/usr/local/genome/lib/genscan/Maize.smat";

# The $0 is the name of the program
my $Usage = "$0 -i InFile.fasta -o Outfile.genscan -g Outfile.gff\n";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my %Options;                  
getopts('i:o:g:h', \%Options);      # Get the options from the command line 

$PrintHelp = $Options{h};
if ($PrintHelp)
{
    print $Usage;
    exit;
}

$InFile = $Options{i} ||
    die "\aERROR: An input file must be provided.\n$Usage";
$OutGen = $Options{o} ||
    $InFile.".genscan";
$OutGff = $Options{g} ||
    $OutGen.".gff";

#-----------------------------+
# CHECK EXISTENCE OF INFILE   |
#-----------------------------+
unless (-e $InFile)
{
    print "\n\aERROR: The input file could not be found at:\n$InFile\n\b";
    exit;
}

#-----------------------------+
# RUN THE GENSCAN PROGRAM     |
#-----------------------------+
print "RUNNING GENSCAN\n";
$cmd = "genscan $Lib $InFile -v > $OutGen";
system ($cmd);

print "CONVERTING GENSCAN TO GFF\n";
Genscan2Gff ($OutGen, $OutGff);

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub Genscan2Gff 
{
# Convert a Genscan format output to the GFF format
# Modified from
# http://nucleus.cshl.org/agsa/presentations.html
# as Downloaded by JCE 02/22/2007

    my $In = $_[0];          # Path to genscan format input file
    my $Out = $_[1];         # Path to gff foramt output file

    my %exon_type = ('Sngl', 'Single Exon',
		     'Init', 'Initial Exon',
		     'Intr', 'Internal Exon',
		     'Term', 'Terminal Exon');
    
    open (IN, "<".$In);
    open (OUT, ">".$Out);
    
    while (<IN>) {
	
	# Last line before predictions contains nothing but spaces and dashes
	if (/^\s*-[-\s]+$/)  {
	    while (<IN>)  {
		my %feature; 
		if (/init|term|sngl|intr/i) {
		    
		    my @f  = split;
		    
		    my ($gene, $exon) = split (/\./, $f[0]); 
		    
                    #name must be a number
		    $feature {name} = $gene + ($exon/1000); 
		    #arrange numbers so that start is always < end
		    if ($f[2] eq '+') {
			$feature {'start'}  = $f[3];
			$feature {'end'}    = $f[4];
			$feature {'strand'} = "+";
		    } elsif ($f[2] eq '-') {
			$feature {'start'}  = $f[4];
			$feature {'end'}    = $f[3];
			$feature {'strand'} = "-";
		    }
		    
		    $feature {'score'}   = $f[12];
		    $feature {'p'}       = $f[11];
		    $feature {'type'}    = $exon_type{$f[1]};
		    $feature {'program'} = 'Genscan';
		    $feature {'program_version'} = '1.0';
		    $feature {'primary'} = 'prediction';
		    $feature {'source'}  = 'genscan';
		    
		    # Pring GFF format output
		    print OUT "$gene.$exon\tgenscan\texon\t" .
			$feature{start} . "\t" . 
			$feature{end} . "\t" . 
			$feature{p} . "\t" . 
			$feature{strand} . "\t.\t" . 
			$gene . "\n";
		} elsif (/predicted peptide/i) {
		    last;   
		}
	    } # End of second while statement
	} # End of if seach command
    } # End of first while statment


} #End of Genscan2Gff subfunction

=head1 HISTORY

STARTED: 04/17/2007

UPDATED: 07/31/2007

=cut
#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 04/17/2007
#  - Program started. Takes infile as variable and produces a
#    genscan output.
#  - Added subfunction to convert the genscan output to 
#    gff format.
