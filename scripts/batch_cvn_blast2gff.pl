#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blast2gff.pl Convert BLAST output to gff 
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/06/2006                                       |
# UPDATED: 08/06/2007                                       |
#                                                           |  
# DESCRIPTION:                                              | 
# Convert blast output to a Apollo compatible gff file.     |
#                                                           |
# USAGE:                                                    |
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

 Usage:
  Name.pl -i InFile -o OutFile

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

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

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

The list of error messages that can be generated,
explanation of the problem
one or more causes
suggested remedies
list exit status associated with each error

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=cut

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
use Bio::SearchIO;             # Parse BLAST output

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my $VERSION = "0.1";

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my $BlastFile = "/home/jestill/projects/asgr/Asgr_c_Net_BigCat/CPAN001/".
    "contig02/contig02_CPAN001.blo";
my $GFFOut = "/home/jestill/projects/asgr/Asgr_c_Net_BigCat/CPAN001/".
    "contig02/contig02_CPAN001.gff";
my $max_e = 0.00001;
my $min_len = 100;

# Set variable scope
my $seqname;
my $infile;
my $outfile;

# Booleans
my $verbose = 0;
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_append = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    "n|name=s"    => \$seqname,
		    # OPTIONS
		    "verbose"     => \$verbose,
		    "append"      => \$do_append,
		    "maxe"        => \$max_e,
		    "minlen"      => \$min_len,
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
		    "q|quiet"     => \$quiet,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ($show_usage) {
    print_help("");
}

if ($show_help || (!$ok) ) {
    print_help("full");
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

#-----------------------------+
# CHECK REQUIRED OPTIONS      |
#-----------------------------+
if ( (!$infile) || (!$outfile) || (!$seqname) ) {
    print_help("full");
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
blast2gff ($infile, $outfile, $do_append, $seqname);

exit;



#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub blast2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # blastin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    my ($blastin, $gffout, $append, $seqname) = @_;
    my $blastprog;        # Name of the blast program used (blastn, blastx)
    my $dbname;           # Name of the database blasted
    my $hitname;          # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit

    # Open the BLAST report object
    my $blast_report = new Bio::SearchIO ( '-format' => 'blast',
					   '-file'   => $blastin,
					   '-signif' => $max_e,
					   '-min_query_len' => $min_len) 
	|| die "Could not open BLAST input file:\n$blastin.\n";
    
    # Open file handle to the gff outfile    
    if ($append) {
	open (GFFOUT, ">>$gffout") 
	    || die "Can not open file:\n $GFFOut\n";
    }
    else {
	open (GFFOUT, ">$gffout") 
	    || die "Can not open file:\n $GFFOut\n";
    }
    
    while (my $blast_result = $blast_report->next_result())
    {
	$blastprog = $blast_result->algorithm;
	$dbname = $blast_result->database_name();

    	while (my $blast_hit = $blast_result->next_hit())
	{

	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {

		my $hitname = $blast_hit->name();

		$strand = $blast_hsp->strand('query');
		
		if ($strand =~ "-1") {
		    $strand = "-";
		}
		elsif ($strand =~ "1") {
		    $strand = "+";
		}
		else {
		    die "Error parsing strand\n";
		}
		
		# Make certain that start coordinate is
		# less then the end coordinate
		if ( $blast_hsp->start() < $blast_hsp->end() ) {
		    $start = $blast_hsp->start();         # Start
		    $end = $blast_hsp->end();
		    #$strand = "+";
		    
		}
		else {
		    #
		    $start = $blast_hsp->end();         # Start
		    $end = $blast_hsp->start();
		    #$strand = "-";
		}


		#-----------------------------+
		# PRINT OUTPUT TO TERMINAL    |
		#-----------------------------+
		if ($verbose) {
		    print
			"$seqname\t".                            # Seqname
			"$blastprog:$dbname\t".                  # Source
			"exon\t".                                # Feature type name
			"$start\t".                              # Start
			"$end\t".                                # End
			$blast_hsp->score()."\t".                # Score
			"$strand\t".
			".\t".                                   
			"$hitname\n";                            # Feature name
		}		
		
		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		# WITH BAC DATA               |
		#-----------------------------+
		# Changing BLASTN to the Bac Name appears to allow 
		# these to be drawn on different levels.
		print GFFOUT 
		    #$blast_hit->name()."\t".
		    "$seqname\t".                            # Seqname
		    "$blastprog:$dbname\t".                  # Source
		    "exon\t".                                # Feature type name
		    "$start\t".                              # Start
		    "$end\t".                                # End
		    $blast_hsp->score()."\t".                # Score
		    "$strand\t".                             # Strand
		    ".\t".                                   # Frame
		    "$hitname\n";                            # Feature name

	    } # End of while next hsp
	} # End of while next hit
    } # End of while next result
    
    close GFFOUT;
    
}

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"cnv_blast2gff.pl -i InFile -o OutFile\n";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS::\n".
	"  --append       # Append results to an existing file\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --quiet        # Run program with minimal output\n";
	
    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
}


=head1 HISTORY

STARTED: 08/06/2007

UPDATED: 08/06/2007

VERSION: $Id: $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/06/2007
# -Program started
