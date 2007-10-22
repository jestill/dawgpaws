#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_tenest2gff.pl - Convert TENest Output to GFF          |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 08/27/2007                                       |
# UPDATED: 08/27/2007                                       |
#                                                           |  
# DESCRIPTION:                                              | 
# Convert blast output to a Apollo compatible gff file.     |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#                                                           |
#-----------------------------------------------------------+
# TODO: May want to consider sources as 
#               tenest:component
#               tenest:nest_group
#               tenest:nest_level
#
=head1 NAME

cnv_tenest2gff.pl - Convert TENest output to GFF

=head1 VERSION

This documentation refers to program version 1.0

=head1 SYNOPSIS
    
  USAGE:
    cnv_tenest2gff.pl -i InFile -o OutFile

=head1 DESCRIPTION

Converts TE Nest output to gff format output.

=head1 COMMAND LINE ARGUMENTS

=over 2 Required Arugments

=item -i,--infile

Path of the input file to convert.

=item -o,--outfile

Path of the gff formatted output file that .

=back

=head2 Additional Information

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
my $VERSION = "1.0";

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+

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
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    "n|name=s"    => \$seqname,
		    # OPTIONS
		    "verbose"     => \$verbose,
		    "append"      => \$do_append,
		    # BOOLEANS
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
    print "ERROR: Input file was not specified\n" if (!$infile);
    print "ERROR: Output file was not specified\n" if (!$outfile);
    print "ERROR: Sequence name was not specified\n" if (!$seqname);
    print_help("full");
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
tenest2gff ($infile, $outfile, $do_append, $seqname);

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub tenest2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # tenestin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    my ($tenestin, $gffout, $append, $seqname) = @_;
    my $tename;           # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit

    my $num_te_data;      # Length of the te data array

    my $i;                # Array index val
    my $j;                # Array index val

    # Initialize counters
    my $numfrag = 0;
    my $numpair = 0;
    my $numsolo = 0;
    my $numnltr = 0;
    my $pair_line = 0;

    # BOOLEANS
    my $in_frag = 0;
    my $in_pair = 0;
    my $in_solo = 0;
    my $in_nltr = 0;

    # SOLO VARS
    # These need to be available outside of the
    # initial observation
    my $sol_entry;
    my $sol_number;
    my $sol_type;
    my $sol_dir;
    my $sol_nest_group;
    my $sol_nest_order;
    my $sol_nest_level;
    
    # PAIR VARS
    my $pair_number;
    my $pair_type;
    my $pair_dir;
    my $pair_bsr;
    my $pair_nest_group;
    my $pair_nest_order;
    my $pair_nest_level;

    # FRAG VARS
    my $frag_number;
    my $frag_type;
    my $frag_dir;
    my $frag_nest_group;
    my $frag_nest_order;
    my $frag_nest_level;

    # NLTR VARS
    my $nltr_number;
    my $nltr_type;
    my $nltr_dir;
    my $nltr_nest_group;
    my $nltr_nest_order;
    my $nltr_nest_level;

    # Open the BLAST report object
    open (TENESTIN,"<$tenestin")
	|| die "Could not open TE-NEST input file:\n$tenestin.\n";
    
    # Open file handle to the gff outfile    
    if ($append) {
	open (GFFOUT, ">>$gffout") 
	    || die "Can not open file:\n $gffout\n";
    }
    else {
	open (GFFOUT, ">$gffout") 
	    || die "Can not open file:\n $gffout\n";
    }
    
    
    while (<TENESTIN>) {
	chomp;                   # Remove line endings
	my @te_data = split;   # Split by spaces
	
	# Four annotation types
	#SOLO (solo LTRs), 
	#PAIR (full LTR retrotransposons), 
	#FRAG (fragmented TEs of all types), 
	#NLTR (Full length non-LTR containing TEs)
 
	# Temp show the line being processed
	#print STDERR "$_\n";

	#-----------------------------------------------------------+
	# ADDITIONAL FEATURE DATA                                   |
	#-----------------------------------------------------------+

	#-----------------------------+
	# ADDITIONAL SOLO DATA        |
	#-----------------------------+
	if ($in_solo) {
      
	    $in_solo = 0;
	    $num_te_data = @te_data;

	    if ($verbose) {
		print STDERR $te_data[0]."\n";
		print STDERR "\tNum in array:\t$num_te_data\n";
	    }

	    # i is used to index in the parent array
	    # j is used to index the sol_coords array where
	    # 1,2,3,4 is seq_start, seq_end, te_start, te_end
	    #
	    my @sol_coords = ();
	    $j = 0;
	    for ($i=1; $i<$num_te_data; $i++) {
		$j++;
		$sol_coords[$j] = $te_data[$i];

		if ($verbose) {
		    print STDERR "\t$i:\t".$te_data[$i]."\n";
		}

		if ($j == 4) {

		    # Appending sol_num to sol type will give a
		    # unique name for every occurrence in the gff file
		    my $sol_name = "solo_".$sol_type."_".$sol_number;
		    # Print output to gff file
		    print GFFOUT 
			"$seqname\t".                # Seqname
			"tenest\t".                  # Source
			"exon\t".                    # Feature type name
			$sol_coords[1]."\t".         # Start
			$sol_coords[2]."\t".         # End
			".\t".                       # Score
			# $sol_dir may give strand
			".\t".                 # Strand
			".\t".                       # Frame
			"$sol_name\n";                # Feature name
		    
		    $j = 0;
		    #$sol_coords=();

		} # End of if $j==4
	    } # End of for $i
	} # End of if $in_solo

	#-----------------------------+
	# ADDITIONAL PAIR DATA        |
	#-----------------------------+
	elsif ($in_pair) {

	    # IS BSR Substitution rate ?

	    # Pair has three additional lines of info
	    $pair_line++;

	    # Get additional info
	    
	    #-----------------------------+
	    # LEFT LTR                    |
	    #-----------------------------+
	    if ($te_data[1] =~ 'L') {
		
		$num_te_data = @te_data;

		if ($verbose) {
		    print STDERR $te_data[0]."\n";
		    print STDERR "\tNum in array:\t$num_te_data\n";
		}

		my @l_pair_coords = ();
		$j = 0;
		for ($i=2; $i<$num_te_data; $i++) {
		    $j++;
		    $l_pair_coords[$j] = $te_data[$i];
		    
		    if ($verbose) {
			print STDERR "\t$i:\t".$te_data[$i]."\n";
		    }
		    
		    if ($j == 4) {
			
			my $l_start;
			my $l_end;
			
			# Make start less then end
			if ($l_pair_coords[1] < $l_pair_coords[2]) {
			    $l_start = $l_pair_coords[1];
			    $l_end = $l_pair_coords[2];
			}
			else {
			    $l_start = $l_pair_coords[2];
			    $l_end = $l_pair_coords[1];
			}

			# Appending sol_num to sol type will give a
			# unique name for every occurrence in the gff file
			my $pair_name = "pair_".$pair_type."_".$pair_number;
			# Print output to gff file
			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "tenest\t".                  # Source
			    "exon\t".                    # Feature type name
			    "$l_start\t".                # Start
			    "$l_end\t".                  # End
			    ".\t".                       # Score
			    # $pair_dir may give strand
			    ".\t".                       # Strand
			    ".\t".                       # Frame
			    "$pair_name\n";              # Feature name
			
			$j = 0;
			#$sol_coords=();
			
		    } # End of if $j==4
		} # End of for $i
		
		
	    }

	    #-----------------------------+
	    # RIGHT LTR                   |
	    #-----------------------------+
	    elsif ($te_data[1] =~ 'R') {


		$num_te_data = @te_data;

		if ($verbose) {
		    print STDERR $te_data[0]."\n";
		    print STDERR "\tNum in array:\t$num_te_data\n";
		}		

		my @r_pair_coords = ();
		$j = 0;
		for ($i=2; $i<$num_te_data; $i++) {
		    $j++;
		    $r_pair_coords[$j] = $te_data[$i];
		    
		    if ($verbose) {
			print STDERR "\t$i:\t".$te_data[$i]."\n";
		    }		    

		    if ($j == 4) {
			
			my $r_start;
			my $r_end;
			
			# Make start less then end
			if ($r_pair_coords[1] < $r_pair_coords[2]) {
			    $r_start = $r_pair_coords[1];
			    $r_end = $r_pair_coords[2];
			}
			else {
			    $r_start = $r_pair_coords[2];
			    $r_end = $r_pair_coords[1];
			}

			# Appending sol_num to sol type will give a
			# unique name for every occurrence in the gff file
			my $pair_name = "pair_".$pair_type."_".$pair_number;
			# Print output to gff file
			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "tenest\t".                  # Source
			    "exon\t".                    # Feature type name
			    "$r_start\t".                # Start
			    "$r_end\t".                  # End
			    ".\t".                       # Score
			    # $pair_dir may give strand
			    ".\t".                       # Strand
			    ".\t".                       # Frame
			    "$pair_name\n";              # Feature name
			
			$j = 0;
			#$sol_coords=();
			
		    } # End of if $j==4
		} # End of for $i
	

	    }
	    
	    #-----------------------------+
	    # MIDDLE                      |
	    #-----------------------------+
	    elsif ($te_data[1] =~ 'M') {
			
		$num_te_data = @te_data;
		if ($verbose) {
		    print STDERR $te_data[0]."\n";
		    print STDERR "\tNum in array:\t$num_te_data\n";
		}		

		my @m_pair_coords = ();
		$j = 0;
		for ($i=2; $i<$num_te_data; $i++) {
		    $j++;
		    $m_pair_coords[$j] = $te_data[$i];
		   
		    if ($verbose) {
			print STDERR "\t$i:\t".$te_data[$i]."\n";
		    }

		    if ($j == 4) {
			
			my $m_start;
			my $m_end;
			
			# Make start less then end
			if ($m_pair_coords[1] < $m_pair_coords[2]) {
			    $m_start = $m_pair_coords[1];
			    $m_end = $m_pair_coords[2];
			}
			else {
			    $m_start = $m_pair_coords[2];
			    $m_end = $m_pair_coords[1];
			}

			# Appending sol_num to sol type will give a
			# unique name for every occurrence in the gff file
			my $pair_name = "pair_".$pair_type."_".$pair_number;
			# Print output to gff file
			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "tenest\t".                  # Source
			    "exon\t".                    # Feature type name
			    "$m_start\t".                # Start
			    "$m_end\t".                  # End
			    ".\t".                       # Score
			    # $pair_dir may give strand
			    ".\t".                       # Strand
			    ".\t".                       # Frame
			    "$pair_name\n";              # Feature name
			
			$j = 0;
			
		    } # End of if $j==4
		} # End of for $i


	    }


	    # END OF PAIR DATA
	    if ($pair_line == 3) {
		$in_pair = 0;
		$pair_line = 0;

		# Print output to gff
	    }

	}

	#-----------------------------+
	# ADDITIONAL FRAG DATA        |
	#-----------------------------+
	elsif ($in_frag) {
	    $in_frag = 0;

	    # Get additinal info
      	    $num_te_data = @te_data;

	    if ($verbose) {
		print STDERR $te_data[0]."\n";
		print STDERR "\tNum in array:\t$num_te_data\n";
	    }

	    # i is used to index in the parent array
	    # j is used to index the sol_coords array where
	    # 1,2,3,4 is seq_start, seq_end, te_start, te_end
	    #
	    my @frag_coords = ();
	    $j = 0;
	    for ($i=1; $i<$num_te_data; $i++) {
		$j++;
		$frag_coords[$j] = $te_data[$i];

		if ($verbose) {
		    print STDERR "\t$i:\t".$te_data[$i]."\n";
		}

		if ($j == 4) {

		    # Appending sol_num to sol type will give a
		    # unique name for every occurrence in the gff file
		    my $frag_name = "frag_".$frag_type."_".$frag_number;
		    # Print output to gff file
		    print GFFOUT 
			"$seqname\t".                # Seqname
			"tenest\t".                  # Source
			"exon\t".                    # Feature type name
			$frag_coords[1]."\t".         # Start
			$frag_coords[2]."\t".         # End
			".\t".                       # Score
			# $sol_dir may give strand
			".\t".                 # Strand
			".\t".                       # Frame
			"$frag_name\n";                # Feature name
		    
		    $j = 0;

		} # End of if $j==4
	    } # End of for $i

	} # End of $in_frag

	#-----------------------------+
	# ADDITIONAL NLTR DATA        |
	#-----------------------------+
	elsif ($in_nltr) {
	    $in_nltr = 0;

	    # Get additinal info
      	    $num_te_data = @te_data;

	    if ($verbose) {
		print STDERR $te_data[0]."\n";
		print STDERR "\tNum in array:\t$num_te_data\n";
	    }

	    # i is used to index in the parent array
	    # j is used to index the sol_coords array where
	    # 1,2,3,4 is seq_start, seq_end, te_start, te_end
	    #
	    my @nltr_coords = ();
	    $j = 0;
	    for ($i=1; $i<$num_te_data; $i++) {
		$j++;
		$nltr_coords[$j] = $te_data[$i];

		if ($verbose) {
		    print STDERR "\t$i:\t".$te_data[$i]."\n";
		}		

		if ($j == 4) {

		    # Appending sol_num to sol type will give a
		    # unique name for every occurrence in the gff file
		    my $nltr_name = "frag_".$nltr_type."_".$nltr_number;
		    # Print output to gff file
		    print GFFOUT 
			"$seqname\t".                # Seqname
			"tenest\t".                  # Source
			"exon\t".                    # Feature type name
			$nltr_coords[1]."\t".         # Start
			$nltr_coords[2]."\t".         # End
			".\t".                       # Score
			# $sol_dir may give strand
			".\t".                 # Strand
			".\t".                       # Frame
			"$nltr_name\n";                # Feature name
		    
		    $j = 0;

		} # End of if $j==4
	    } # End of for $i

	} # End of $in_nltr



	#-----------------------------------------------------------+
	# NEW RECORD STARTING                                       |
	#-----------------------------------------------------------+

	#-----------------------------+
	# SOLO                        |
	#-----------------------------+
	if ($te_data[0] =~ 'SOLO') {
	    $numsolo++;
	    $in_solo = 1;              # Flip boolean to true

	    #$sol_entry = $te_data[0]; 
	    $sol_number = $te_data[1];
	    $sol_type = $te_data[2];
	    $sol_dir = $te_data[3];
	    $sol_nest_group = $te_data[4];
	    $sol_nest_order = $te_data[5];
	    $sol_nest_level = $te_data[6];

	}

	#-----------------------------+
	# PAIR                        |
	#-----------------------------+
	elsif ($te_data[0] =~ 'PAIR') {
	    $numpair++;
	    $in_pair = 1;              # Flip boolean to true

	    $pair_number = $te_data[1];
	    $pair_type = $te_data[2];
	    $pair_dir = $te_data[3]; 
	    $pair_bsr = $te_data[4];
	    $pair_nest_group = $te_data[5];
	    $pair_nest_order = $te_data[6]; 
	    $pair_nest_level = $te_data[7];

	}

	#-----------------------------+
	# FRAG                        |
	#-----------------------------+
	elsif ($te_data[0] =~ 'FRAG') {
	    $numfrag++;
	    $in_frag = 1;              # Flip boolean to true

	    $frag_number = $te_data[1]; 
	    $frag_type = $te_data[2]; 
	    $frag_dir = $te_data[3];
	    $frag_nest_group = $te_data[4]; 
	    $frag_nest_order = $te_data[5]; 
	    $frag_nest_level = $te_data[6];

	}
	
	#-----------------------------+
	# NLTR                        |
	#-----------------------------+
	elsif ($te_data[0] =~ 'NLTR') {
	    $numnltr++;
	    $in_nltr = 1;              # Flip boolean to true

	    #Non-LTR entry, 
	    $nltr_number = $te_data[1]; 
	    $nltr_type = $te_data[2]; 
	    $nltr_dir = $te_data[3]; 
	    $nltr_nest_group = $te_data[4]; 
	    $nltr_nest_order = $te_data[5]; 
	    $nltr_nest_level = $te_data[6];

	}

    } # End of while TESTIN
    

    # Show summary of counts if vebose
    if ($verbose) {
	print STDERR "\n";
	print STDERR "NUM SOLO:\t$numsolo\n";
	print STDERR "NUM PAIR:\t$numpair\n";
	print STDERR "NUM FRAG:\t$numfrag\n";
	print STDERR "NUM NLTR:\t$numnltr\n";
    }

    close GFFOUT;
    
}

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"cnv_tenest2gff.pl -i InFile -o OutFile -n SeqName\n";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"  --name         # Name of the reference sequence\n".
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

STARTED: 08/27/2007

UPDATED: 08/27/2007

VERSION: $Id: $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/27/2007
# -Program started
