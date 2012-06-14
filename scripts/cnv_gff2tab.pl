#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_gff2tab.pl - Converts GFF3 file to table              |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/14/2012                                       |
# UPDATED: 06/14/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert a GFF3 file to a tab delimited table with the    |
#  original data and a separate column for each value in    |
#  the attribute column. This assumes that each value in    |
#  the attribute field has a unique name.                   |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
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

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

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
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "gff-ver=s"   => \$gff_ver,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);



 
#-----------------------------+
# STANDARDIZE GFF VERSION     |
#-----------------------------+
unless ($gff_ver =~ "GFF3" || 
	$gff_ver =~ "GFF2") {
    # Attempt to standardize GFF format names
    if ($gff_ver =~ "3") {
	$gff_ver = "GFF3";
    }
    elsif ($gff_ver =~ "2") {
	$gff_ver = "GFF2";
    }
    else {
	print "\a";
	die "The gff-version \'$gff_ver\' is not recognized\n".
	    "The options GFF2 or GFF3 are supported\n";
    }
}


#-----------------------------+
# PRINT REQUESTED HELP        |
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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

if ($infile) {
    open(INFILE, "<$infile") ||
	die "Can not open infile $infile";
}
else {
    open (INFILE, "<&STDIN") ||
	die "Can not open STDIN for input";
}

if ($outfile) {
    open (OUTFILE, ">$outfile") ||
	die "Can not open outfile $outfile"
}
else {
    open (OUTFILE, ">&STDOUT") ||
	die "Can not STDOUT";
}

# Entie gff file in columns
my @gff;
# list of all keys uses
my @keys;
# Load the entire gff file to @gff
my $i = 0;
while (<INFILE>) {

    chomp;

    # Skip existing comments or pragmas
    next if m/^\#/; 
  

    my @gff_parts = split(/\t/, $_);

    # Get attributes as key value pairs
    my %attributes = ();
    my @attribute_parts = split(/\;/,  $gff_parts[8]);
    my $num_attributes = @attribute_parts;

    print STDERR $_."\n" if $verbose;
    print STDERR "\tatr:".$num_attributes."\n" if $verbose;
    foreach my $attribute (@attribute_parts) {
	my ($key,$value) = split(/\=/, $attribute);
	print STDERR "\t$key->$value\n" if $verbose;
	
	# Push the key used to the entire set of key value pairs
	push(@keys,$key);
	
	# push key value to the attributes hash
	$attributes{ $key } = $value;
	
    }
    
    # Push gff parts to array
    $gff[$i][0] = $gff_parts[0];
    $gff[$i][1] = $gff_parts[1];
    $gff[$i][2] = $gff_parts[2];
    $gff[$i][3] = $gff_parts[3];
    $gff[$i][4] = $gff_parts[4];
    $gff[$i][5] = $gff_parts[5];
    $gff[$i][6] = $gff_parts[6];
    $gff[$i][7] = $gff_parts[7];
    $gff[$i][8] = $gff_parts[8];
    $gff[$i][9] = \%attributes;
    
    # Increment iterator
    $i++;

}

# Get uniqu list of attributes used in the gff file
my @uniq_attributes = keys %{{ map { $_ => 1 } @keys }};

if ($verbose) {
    print STDERR "UNIQUE KEYS:\n";
    foreach my $uniq_key (@uniq_attributes) {
	print "\t".$uniq_key."\n";
    }
}


# print to tab delimited output file
#-----------------------------+
# Print outfile header        |
#-----------------------------+
print OUTFILE "seqid\t".
    "source\t".
    "type\t".
    "start\t".
    "end\t".
    "score\t".
    "strand\t".
    "phase\t".
    "attribute";
foreach my $uniq_key (@uniq_attributes) {
    print OUTFILE"\t".$uniq_key;
}
print OUTFILE "\n";

#-----------------------------+
# Print each row              |
#-----------------------------+
foreach my $row (@gff) {
    my $col_num = 0;
    foreach my $col (@$row) {

	if ($col_num == 9) {
	    foreach my $uniq_key (@uniq_attributes) {
		# If a key value exists 
		if ( exists $col->{ $uniq_key } ) {
		    print OUTFILE "\t".$col->{$uniq_key};
		}
		else {
		    print OUTFILE "\t.";
		}
		$j++;
	    }

	}
	elsif($col_num == 8 ) {
	    print OUTFILE $col;
	}
	else {
	    print OUTFILE $col."\t";
	}
	$col_num++;
	
    }
    print OUTFILE "\n";

}

close (INFILE);
close (OUTFILE);

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

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    Name.pl -i InFile -o OutFile

=head2 Required Arguments

    --infile        # Path to the input file
    --outfie        # Path to the output file

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

=head1 EXAMPLES

The following are examples of how to use this script

=head2 Typical Use

This is a typcial use case.

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 REFERENCE

A manuscript is being submitted describing the DAWGPAWS program. 
Until this manuscript is published, please refer to the DAWGPAWS 
SourceForge website when describing your use of this program:

JC Estill and JL Bennetzen. 2009. 
The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes.
http://dawgpaws.sourceforge.net/

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
