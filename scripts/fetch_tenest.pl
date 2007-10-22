#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fetch_tenest.pl - Fetch TE Nest LTR and SVG files         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/29/2007                                       |
# UPDATED: 08/29/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Workflow to fetch the tenest program results and covert  |
#  the output to gff format for use in the DAWG-PAWS        |
#  pipeline.                                                |
#                                                           |
# VERSION:                                                  |
# $Id:: fetch_tenest.pl 85 2007-08-29 14:29:27Z JamesEst#$: |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TODO: Update POD
#       Fix help function
#       Add dir to the tenest2gff subfunction

=head1 NAME

fetch_tenest.pl - Fetch TE Nest LTR and SVG files

=head1 VERSION

This documentation refers to $Rev$

=head1 SYNOPSIS

  USAGE:
    fetch_tenest.pl -i SeqList.txt -o OutDir [--gff]

=head1 DESCRIPTION

Given a tab delimited input file of theFetches the TE Nest results using wget.
The results are converted to gff and copied to the gff dir if requested.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. This is a tab delimited input file that must be
in the following format:

    HEX0350E24     1188334554_18873
    HEX0358K19     1188334578_19208
    HEX0411K13     1188334604_19821

Where the first column is the name of the sequence that was analyzed
by TE Nest and the second column is the id of the TE Nest job. This job 
id is returned by TE Nest when the job is initially submitted.

=item -o,--outdir

Path of the output directory. This is the base dir that the analysis 
results are stored in. 

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

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # Sanity checker
use File::Copy;                # Copy files 
use Getopt::Long;              # Get options from the command line
use LWP::Simple;               # Used to get files

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+

my @file_list = ();                 # Files to process
                                    # Index starts at zero

my $infile;
my $outdir;

# Booleans
my $help = 0;
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_gff_convert = 0;

# Counters
my $i = 0;                         # Index for file_list array
my $max_i;                         # Max valu in the file_list array
my $num_files;                     # The number of files to process

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# Required variables
		    "i|infile=s"  => \$infile,
                    "o|outdir=s" => \$outdir,
		    # Additional options
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    "gff"         => \$do_gff_convert,
		    # Additional information
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);


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
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}


#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print STDERR "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

#-----------------------------+
# LOAD DATA TO FETCH FROM     |
# TAB DELIM TEXT FILE         |
#-----------------------------+
open (INFILE, "<$infile")
    || die "Can not open the input file:\n$infile\n";

while (<INFILE>) {
    chomp;

# Use the following if an error crops up later
#    print "IN: $_\n";

    # Ignore comment lines starting with #
    unless (m/^\#/) {
	my @in_line = split(/\t/, $_);
	my $num_in_line = @in_line; 
	$file_list[$i][0] = $in_line[0];
	$file_list[$i][1] = $in_line[1];
	
# Use the following if a error crops up later
#	if ($verbose) {
#	    print STDERR "\tProcessed as:\n";
#	    print STDERR "\t".$file_list[$i][0]."\n";
#	    print STDERR "\t".$file_list[$i][1]."\n";
#	}
	
	$i++;
    }

} # End of while INFILE

$num_files = $i;
$max_i = $i-1;

close INFILE;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INFILE    |
#-----------------------------+
if ($num_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No files were found in the input file\n".
	"$infile\n".
    exit;
}

print STDERR "\n$num_files files will be processed\n" if $verbose;

#-----------------------------------------------------------+
# PROCESS EACH FILE IN @file_list                           |
#-----------------------------------------------------------+
for ($i=0; $i<$num_files; $i++) {

    my $file_num = $i+1;
    print STDERR "Processing $file_num of $num_files\n" if $verbose;
    
    # Load array data to variables
    my $name_root = $file_list[$i][0] ||
	die "ERROR: file_list array error in name for file $file_num\n";
    my $tenest_id = $file_list[$i][1] ||
	die "ERROR: file_list array error in tenest id for file $file_num\n";
    
    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE TENEST OUTDIR        |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $tenest_dir = $name_root_dir."tenest/";
    unless (-e $tenest_dir) {
	mkdir $tenest_dir ||
	    die "Could not create genscan out dir:\n$tenest_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $gff_dir = $name_root_dir."gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create genscan out dir:\n$gff_dir\n";
    }

    #-----------------------------+
    # FETCH TE NEST - LTR Results |
    #-----------------------------+
    # This is the text file produced by TE Nest
    my $ltr_url = "http://www.plantgdb.org/tmp/qry-$tenest_id.LTR";
    my $ltr_local = "$tenest_dir$tenest_id.LTR";
    my $ltr_rename = "$tenest_dir$name_root.LTR";

    print STDERR "\tFetching: $ltr_url\n\tTo:$ltr_local\n" if $verbose;
    my $ltr_res = getstore($ltr_url, $ltr_local);
    print STDERR "\tResult:$ltr_res\n" if $verbose;

    # If the http response is 200 (ok) then work with the local copy
    if ($ltr_res =~ "200") {
	# Make a copy of the ltr file using the root name
	# This will facilitate fetching the name later

	copy ( $ltr_local, $ltr_rename ) ||
	    print STDERR "Error copying file $ltr_local\n";
	
	# Convert output to gff if requested
	if ($do_gff_convert) {

	    my $do_append = 0;
	    my $ltr_gff_out = $tenest_dir.$name_root."_tenest.gff"; 
	    tenest2gff ( $ltr_rename, $ltr_gff_out, $do_append, $name_root);

	    my $ltr_gff_copy = $gff_dir.$name_root."_tenest.gff";	    
	    copy ( $ltr_gff_out, $ltr_gff_copy ) ||
		print STDERR "Error copying file $ltr_gff_out\n";    
	}

    }
    else {
	print STDERR "Error fetching file: $ltr_url\nERR:$ltr_res\n";
    }
    
    #-----------------------------+
    # FETCH TE NEST - SVG Results |
    #-----------------------------+
    # This is the image produced by TE Nest
    my $svg_url = "http://www.plantgdb.org/tmp/qry-$tenest_id.svg";
    my $svg_local = "$tenest_dir$tenest_id.svg";
    my $svg_rename = "$tenest_dir$name_root.svg.xml";

    print STDERR "\tFetching: $svg_url\n\tTo:$svg_local\n" if $verbose;
    my $svg_res = getstore($svg_url, $svg_local);
    print STDERR "\tResult:$svg_res\n" if $verbose;

    # If the http response is 200 (ok) then work with the local copy
    if ($svg_res =~ "200") {

	copy ($svg_local, $svg_rename) ||
	    print STDERR "Error copying file $svg_local\nto$svg_rename";

    }
    else {
	# http responses listed at
	# http://www.w3.org/Protocols/rfc2616/rfc2616-sec10.html
	print STDERR "Error fetching file $svg_url\nERR:$ltr_res";
    }

} # End of increment $i

print STDERR "Process completed\n" if $verbose;

# Create tenest dir if it does not exist

# COPY FILES TO NEW NAME
# Use copy instead of move to keep track of what the original
# file id was

# 

# Convert te_nest output to gff format

# copy gff file to gff dir

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"fetch_tenest.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS:\n".
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


sub tenest2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # tenestin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    my ($tenestin, $gffout, $append, $seqname) = @_;
    my $tename;          # Name of the hit
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


=head1 HISTORY

STARTED: 08/29/2007

UPDATED: 08/29/2007

VERSION: $Id: fetch_tenest.pl 85 2007-08-29 14:29:27Z JamesEstill $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/29/2007
# - Program started
# - Added tenest2gff from cnv_tenest2gff.pl
