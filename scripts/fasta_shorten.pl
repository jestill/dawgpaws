#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_short.pl - Give fasta file shorter names            |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/17/2007                                       |
# UPDATED: 12/10/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#   Change headers in a fasta file to give shorter names    |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                     # Input dir of fasta files to shorten
my $outdir;                    # Output dir for shortened fasta files
my $new_len = 20;              # New length of the header
my $head_new;                  # New, shorter header
my $test_len;                  # Test length of the header
my $cur_len;                   # Current length of the header

# Booleans
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $test = 0;
my $uppercase = 0;            # Convert sequence strigs to uppercase

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required options
		    "i|indir=s"   => \$indir,
                    "o|oudir=s"   => \$outdir,
		    # Additional options
		    "l|length=s"  =>, \$new_len,
		    # Booleans
		    "uppercase"   => \$uppercase,
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
		    "test"        => \$test,
		    "q|quiet"     => \$quiet,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( ($show_usage) ) {
#    print_help ("usage", File::Spec->rel2abs($0) );
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
#    print_help ("help",  File::Spec->rel2abs($0) );
    print_help ("help",  $0 );
}

if ($show_version) {
    print "\n$0:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
print "\n" unless $quiet;

# The test length is the length of the header plus one
#$test_len = $new_len + 1;
$test_len = int($new_len + 1);

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
print "Creating output dir ...\n" unless $quiet;
unless (-e $outdir) {
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# PROCESS EACH INDIVIDUAL     |
# FILE IN THE DIRECTORY       |
#-----------------------------+
for my $ind_file (@fasta_files) {

    my $infile = $indir.$ind_file;
    my $outfile = $outdir.$ind_file;

    open (IN, $infile) ||
	die "Can not open infile:\n$infile\n";
    
    open (OUT, ">".$outfile) ||
	die "Can not open outfile:\n$outfile\n";

    while (<IN>) {
	
	chomp;
	# If the first characeter is a >
	if (/^\>/) {
	    #-----------------------------+
	    # FASTA HEADER PROCESSING     |
	    #-----------------------------+
	    
	    print " processing $_\n";

	    #$cur_len = length($_);
	    $cur_len = int(length($_));

	    print "\tLEN: $cur_len\n";
	    if ( $cur_len > $test_len ) {
		# If the fasta header is longer then the desired length then
		# create a shortened header print to the outfile
		# We will start at 1 instead of zero to ignore
		# the > character
		$head_new = substr ( $_, 1, $new_len );
		print "\tNEW: >$head_new\n" unless $quiet;
		print OUT ">$head_new\n";
	    } 
	    else {
		print OUT "$_\n";
		# print the header to the outfile unchanged
		
	    } # End of if header is too long

	}
	
	else {
	    #-----------------------------+
	    # SEQUENCE STRING PROCESSING  |
	    #-----------------------------+
	    # print the string to the outfile unchanged if this is
	    # not a fasta header line
	    # T
	    #if ($)
	    if ($uppercase) {
		# Convert the sequence string to uppercase
		# I will try this as a bare regexp since 
		# we are workign with $_
		tr/a-z/A-Z/;
		print OUT "$_\n";

	    }
	    else {
		# If not converting the sequence string to uppercase
		# print the output file unchanged.
		print OUT "$_\n";
	    }
	}


    } # End of while IN

} # End of for each file in the directory

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {
    my ($help_msg, $podfile) =  @_;
    # help_msg is the type of help msg to use (ie. help vs. usage)
    
    print "\n";
    
    #-----------------------------+
    # PIPE WITHIN PERL            |
    #-----------------------------+
    # This code made possible by:
    # http://www.perlmonks.org/index.pl?node_id=76409
    # Tie info developed on:
    # http://www.perlmonks.org/index.pl?node=perltie 
    #
    #my $podfile = $0;
    my $scalar = '';
    tie *STDOUT, 'IO::Scalar', \$scalar;
    
    if ($help_msg =~ "usage") {
	podselect({-sections => ["SYNOPSIS|MORE"]}, $0);
    }
    else {
	podselect({-sections => ["SYNOPSIS|ARGUMENTS|OPTIONS|MORE"]}, $0);
    }

    untie *STDOUT;
    # now $scalar contains the pod from $podfile you can see this below
    #print $scalar;

    my $pipe = IO::Pipe->new()
	or die "failed to create pipe: $!";
    
    my ($pid,$fd);

    if ( $pid = fork() ) { #parent
	open(TMPSTDIN, "<&STDIN")
	    or die "failed to dup stdin to tmp: $!";
	$pipe->reader();
	$fd = $pipe->fileno;
	open(STDIN, "<&=$fd")
	    or die "failed to dup \$fd to STDIN: $!";
	my $pod_txt = Pod::Text->new (sentence => 0, width => 78);
	$pod_txt->parse_from_filehandle;
	# END AT WORK HERE
	open(STDIN, "<&TMPSTDIN")
	    or die "failed to restore dup'ed stdin: $!";
    }
    else { #child
	$pipe->writer();
	$pipe->print($scalar);
	$pipe->close();	
	exit 0;
    }
    
    $pipe->close();
    close TMPSTDIN;

    print "\n";

    exit 0;
   
}

1;
__END__

=head1 NAME

fasta_shorten.pl - Change headers in a fasta file to give shorter names.

=head1 VERSION

This documentation refers to fasta_shorten version $Rev$

=head1 SYNOPSIS

=head2 Usage

    fasta_shorten.pl -i InDir -o OutDir

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

This program will take all of the fasta files in a input directory 
and will shorten the name in the fasta header file. This is primarily
developed for the wheat project. The name used in the fasta header
is the name used by apollo as well as RepeatMaksker.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item -l

New length. The fasta header will be shortened to this length.
Default length is 20.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not currently make use of configuration files
or settings in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

No external software is currently required to use this program

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Short names limited to a length

There is currently no feature to support any sort of shorted names
based on a delimiting character. You are therefore limited to setting
a length for making the fasta file shorter.

=back

=head1 SEE ALSO

The fasta_shorten.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html   

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 07/17/2007

UPDATED: 12/10/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/17/2007
# - Program started.
#
# 12/10/2007
# - Moved POD docs to end of script
# - Changed program version reference to SVN revision
# - Added new help subfunction that extracts help and usage
#   text from POD documentation
# - Added check for required arguments
