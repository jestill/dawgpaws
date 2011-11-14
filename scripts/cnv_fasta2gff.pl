#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_get_fasta_size.pl - Get length of fasta files       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/02/2011                                       |
# UPDATED: 08/02/2011                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given an input directory of FASTA files this will        |
#  generate a tab delimited text file that gives the        |
#  size of all sequences in the input directory. Ase        |
#  seqname <tab> seqsize                                    |
#                                                           |
# USAGE:                                                    |
#  batch_get_fasta_size.pl -i fasta_dir/ -o size.txt        |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
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
use Bio::SeqIO;                # Read and write seq files in different formats

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outfile;
my $infile_format = "fasta";

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "format"      => \$infile_format,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

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

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\nbatch_get_fasta_size.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( !$indir ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    
    print_help ("usage", $0 );
}



#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
# Expect input as dir but will accept single files
my @fasta_files;
if (-d $indir) {
    opendir( DIR, $indir ) || 
	die "Can't open directory:\n$indir"; 
    @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
    closedir( DIR );
    
    # add slash
    unless ($indir =~ /\/$/ ) {
	$indir = $indir."/";
    }

}
elsif (-f $indir) {
    print STDERR "Expecting this to be a file.\n";
    push (@fasta_files, $indir);
}


my $count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

# OUTPUT FILE
if ($outfile) {
    open(TXTOUT, $outfile) ||
	die "Can not open output file $outfile\n";
}
else {
    open (TXTOUT, ">&STDOUT") ||
	die "Can not print to STDOUT\n";
}


my $fasta_file_num = 0;
my $seq_count = 0;

for my $ind_seq_file (@fasta_files) {
    
    # Get input file path by adding dir, with exception
    # for situations where indir was a file
    my $infile;
    if (-d $indir) {
	$infile = $indir.$ind_seq_file;
    }
    elsif (-f $indir) {
	$infile = $indir;
    }

    my $seq_in;
    # OPEN THE INPUT FILE
    if ($infile) {
	$seq_in = Bio::SeqIO->new(-file   => "<$infile",
				  -format => $infile_format) 
	    || "Can not connect to input file $infile\n";
    }
    
    # Cycle through each sequence record in the file
    # and report the sequence name and 
    while (my $inseq = $seq_in->next_seq) {
	
#	print TXTOUT $inseq->primary_id."\t".$inseq->length()."\n";

	# Print the output as GFF file for scaffold feature in gbrowse
	# database
	#
	# TO DO .. allow name as chromosome
	# 

	print TXTOUT $inseq->primary_id."\t".    #  1
	    "dawgpaws\t".                        # 2. program
	    "scaffold\t".                        # 3. feature type
	    "1\t".                               # 4. feature start
	    $inseq->length()."\t".               # 5. feature end
	    ".\t".                               # 6. score
	    ".\t".                               # 7. strandgff
	    ".\t".                               # 8. frame
	    "Name=".$inseq->primary_id."\n";    # 9. attribute
	
	$seq_count++;
	
    }

}

exit 0;

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

batch_get_fasta_size.pl - Short program description. 

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    batch_get_fasta_size.pl -i indir -o outfile

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outfile

Path of the output file. If no outfile is given, output will be written to 
STDOUT.

=item -c, --config

Path to a config file. This is a tab delimited text file
indicating the required information for each of the databases to blast
against. Lines beginning with # are ignored.

=back

=head1 OPTIONS

=over 2

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item --verbose

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

=back

=head1 Additional Options

=over 2

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

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of configuration files or variables
set in the user environment.

=head1 DEPENDENCIES

=head2 Required Perl Modules

=over

=item * Bio::SeqIO

The Bio::SeqIO module from bioperl is required for this program.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Known Limitation

This program will currently only find fasta files in the input dir
that have the fasta or fa extension.

=back

=head1 SEE ALSO

The program is part of the DAWGPAWS package of genome
annotation programs. See the DAWGPAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/02/2011

UPDATED: 08/02/2011

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
