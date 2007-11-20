#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_fgenesh2gff.pl - Convert FGENESH output to gff foramt |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill, Renyi Lu                        |
# STARTED: 11/28/2006                                       |
# UPDATED: 10/30/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  A parser of the Fgenesh contig-by-contig results.        |
#  Jamie modified this from something that Renyi had.       |
#  Input is a text file output from the FGENESH gene        |
#  prediction program from softberry. Output is a gff file  |
#  of the predicted gene-spans and gene models.             |
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;

#-----------------------------+
# VARIABLES                   |
#-----------------------------+
my $insertDB = "";
my $PrintGeneSpan = "";
my $printScreen = 1;
my $debug = 0;

my $usage = "$0 input_directory output_directory\n";

my $gffOutPath = "/home/jestill/projects/wheat_annotation/Hex/".
    "fgenesh_out/GffText.gff.txt";
my $inputDir = "/home/jestill/projects/wheat_annotation/Hex/fgenesh_out";

#-----------------------------+
# USE INPUT FROM COMMAND LINE |
# IF THE INPUT DIR IS         |
# SPECIFIED FROM THE COMMAND  |
# LINE                        |
#-----------------------------+
# This will override any predefined variables from above.
if(@ARGV > 0)
{
    $inputDir = $ARGV[0];
    $gffOutPath = $ARGV[1];
}

#-----------------------------+
# OPEN THE OUTPUT FILES       |
#-----------------------------+
open (GFFOUT, ">".$gffOutPath) ||
    die "Can not open output file ".$gffOutPath."\n";

#-----------------------------+
# OPEN THE INPUT DIR          |
#-----------------------------+
chdir $inputDir;
opendir (DIR, $inputDir) || 
    die "Cannot open input directory:\n $inputDir\n $!";

my @files = grep /txt/, readdir DIR;

my $num_files = @files;
print STDERR "$num_files files to process\n";

my $counter = 0;

foreach  my $file (@files) 
{

    $counter++;

    print STDERR "$counter\t$file\n";

    #-----------------------------+
    # INITIALIZE VARIABLES        |
    #-----------------------------+
    my $features = {};
    my $proteins = {};
    my $d = {};

    $d->{parsing_predicted} = 0;
    $d->{parsing_protein} = 0;
    $d->{parsing_mRNA} = 0;

    #-----------------------------+
    # OPEN THE FILE FOR READING   |
    #-----------------------------+ 
    open (IN,"$file") || 
	die "Cannot open $file: $!";
    
    while (my $line = <IN>) 
    {
	chomp $line;

	next if ($line eq ''); # skip the empties

	$d->{clone}     = ($line =~ /Seq\sname:\s(\S+)\s(.*)/) ? $1 : $d->{clone};
	
	#-----------------------------+
	# LENGTH OF SEQUENCE          |
	#-----------------------------+
	if ($line =~ /Length\sof\ssequence:\s(\d+)/) 
	{
	    $d->{length} = $1;
	    #$d->{gc}     = 0;
	}

	#-----------------------------+
	# GENE NUMBER INFORMATION     |
	#-----------------------------+
	if ($line=~/.*predicted\sgenes\s(\d+)\sin\s\+chain\s(\d+)\sin\s-chain\s(\d+).*/) 
	{
	    if($debug){
		#	print STDERR "Feature line: $line\n";
	    }
	    $d->{total_genes}  = $1;
	    $d->{total_genes_plus}  = $2;
	    $d->{total_genes_minus} = $3;
	}
	

	if ($line=~/.*predicted\sexons\s(\d+)\sin\s\+chain\s(\d+)\sin\s-chain\s(\d+).*/) 
	{
	    if($debug){
		print STDERR "Feature line: $line\n";
	    }
	    $d->{total_exons}  = $1;
	    $d->{total_exons_plus}  = $2;
	    $d->{total_exons_minus} = $3;
	}
	
	#-----------------------------+
	# DO THE FOLLOWING IF WE HAVE |
	# READ THE HEADER LINE THAT   |
	# INDICATES THE START OF DATA |
	#-----------------------------+
	if ($d->{flag}) 
	{
	    # Parse each line - store multiple information together by gene
	    # Store each predicted protein in a hash of arrays of arrays
	    
	    if ($line =~ /\s*(\d+)\s+([+-])\s+(\d+)\s+(\w+)\s+(\d+)\s+-\s+(\d+)\s+(\d+\.\d+)\s+(\d+)\s+-\s+(\d+)\s+(\d+).*/) 
	    {
		if($debug)
		{
		    # print STDERR "Feature list: $line\n";
		}
		
		# This parses correctly for the most part
		
		my $id = $1;
		# Each feature for this gene should be in a hash
		# order: clone,strand,exon_id, feature,start,end,score,orf_start,orf-end,length)
		push (@{$features->{$id}},[ $d->{clone},$2,$3,$4,$5,$6,$7,$8,$9, $10 ]);
	    } elsif ($line =~ /\s*(\d+)\s+([\+-])\s+(\w+)\s+(\d+)\s+(\d+\.\d+)/) {
		if($debug)
		{
		    print STDERR "Feature list: $line\n";
		}
		# some records look like this
		# 1 -     PolA     617                1.05  
		my $id = $1;
		# Each feature for this gene should be in a hash
		# order: clone,strand,exon_id, feature,start,end,score,orf_start,orf-end,length)
		push (@{$features->{$id}},[ $d->{clone},$2,'', $3,$4,'',$5,'','','' ]);
		
	    } else {}
	}
	

	#-----------------------------+
	# DO THE FOLLOWING IF WE      |
	# ARE PARSING PREDICTED       |
	# PROTEINS                    |
	#-----------------------------+
	if ($d->{parsing_predicted}) 
	{
	    # first line is the identifier line
	    if($line =~ /.*FGENESH:\s+(\d+)\s+(\d+)\s+exon\s+\(s\)\s+(\d+)\s+-\s+(\d+)\s+(\d+).*chain\s+([\+-])/)
	    {
		
		if($debug)
		{
		    print STDERR "protein line $line\n";
		}

		$d->{parsing_protein} = 1;
		$d->{parsing_mRNA} = 0;
		
		# Parse out each of the fields
		# method id exons position aa strand
		# still need to account for positions, aa, and strand
		# >FGENESH   1   8 exon (s)    203  -   3137    601 aa, chain +
		
		# This regex isn't working
		#$line =~ />FGENESH\s+(\d+)\s+(\d+)\sexon.*(\d+)\s+-\s+(\d+)\s+(\d+)\saa,\s+chain\s([+-])/;
		
		my $id = $1;
		$proteins->{$id}->{clone}      = $d->{clone};
		$proteins->{$id}->{method}     = 'FGENESH';
		$proteins->{$id}->{fgenesh_id} = $1;
		$proteins->{$id}->{exons}      = $2;
		$proteins->{$id}->{start}      = $3;
		$proteins->{$id}->{stop}       = $4;
		$proteins->{$id}->{plength}    = $5;
		$proteins->{$id}->{strand}     = $6;
		
		# Need to store the ID for later - storing sequences
		$d->{current_id} = $id;
		
	    } 
	    
	    #-----------------------------+
	    # FOR SEQUENCE LINES          |
	    #-----------------------------+
	    if($d->{parsing_protein} && (!$d->{parsing_mRNA}) && ($line !~ /FGENESH/))
	    {  		
                # keep concatenating sequence lines
		$proteins->{$d->{current_id}}->{sequence} .= $line;
	    }
	    
	    #-----------------------------+
	    # FLAG IF WE ARE PARSING mRNA |
	    #-----------------------------+
	    if($line =~ /FGENESH:\s*\[mRNA\]/)
	    {
		$d->{parsing_protein} = 0;
		$d->{parsing_mRNA} = 1;
	    }
	    
	} # End of parsing predicted proteins
	
	#-----------------------------+
	# REMAINING LINES ARE THE     |
	# SEQUENCE                    |
	#-----------------------------+
	# Am I done with the data set?
	if ($line =~ /Predicted\s+protein/) 
	{
	    $d->{parsing_predicted} = 1;
	    delete $d->{flag};
	}
	
	#-----------------------------+
	# STR FEATURES                |
	#-----------------------------+
	# Extracting the salient information is going to be a pain
	# Need to find out what exactly all of the fields are.
	if ($line =~ /G\s+Str\s+Feature\s+/) 
	{
	    # We've arrived at the features.  Set a flag.
	    $d->{flag}++;
	    if($debug){
		print STDERR "Entering feature list section.\n";
	    }
	}

    } # next line
    
    # Now let's print everything:
    # First, the easy stuff.  Clone info.  
    #delete $d->{parsing_predicted};
    #delete $d->{current_id};
    #delete $d->{parsing_protein};
    #delete $d->{parsing_mRNA};
    
    # Insert the contig info
    my @bac_info = ($file, $d->{clone}, $d->{length}, $d->{total_genes}, 
		    $d->{total_genes_plus}, $d->{total_genes_minus}, $d->{total_exons},
		    $d->{total_exons_plus}, $d->{total_exons_minus});
#    if($printScreen)
#    {
#	print "Insert bac_info_fgenesh: @bac_info\n";
#    }
    
#    if($insertDB)
#    {
#	$bac_sth->execute(@bac_info);
#    }
 
    foreach my $id (sort keys %$proteins) 
    {
	my $ref = $proteins->{$id};
	#foreach my $f (keys %{$proteins->{$id}}) {
	#  $vals{$f} = $proteins->{$id}->{$f};
	#}
	
	my @vals = ($file, $ref->{clone}, $ref->{fgenesh_id}, $ref->{method},
	            $ref->{strand}, $ref->{exons}, $ref->{plength}, $ref->{start},
		    $ref->{stop}, $ref->{sequence});

	if($printScreen)
	{
	    print "Insert protein_fgenesh @vals\n";
	}

	#-----------------------------+
	# PRINT DATA TO THE GFF       |
	# OUTPUT FILE                 |
	# 11/28/2006 - JCE            |
	#-----------------------------+
	# I believe that this is the overall gene spans at this point
	# and not the actual full set of exons
       
	if ($PrintGeneSpan)
     	{       
	    print GFFOUT $ref->{clone}."\t";                      # SEQNAME
	    print GFFOUT $ref->{method}."\t";                     # SOURCE
	    print GFFOUT "GeneSpan\t";      # FEATURE
	    print GFFOUT $ref->{start}."\t";                      # START
	    print GFFOUT $ref->{stop}."\t";                       # END
	    print GFFOUT ".\t";                                   # SCORE
	    print GFFOUT $ref->{strand}."\t";                     # STRAND
	    print GFFOUT ".\t";                                   # FRAME
	    print GFFOUT "FGENESH-".$ref->{fgenesh_id}."-GeneSpan\t";      # ATTRIBUTE
	    print GFFOUT "\n";
	}

	my $pred_id = 0;

#	if($insertDB)
#	{
#	    $protein_sth->execute(@vals);
#	    $pred_id = $dbh->{'mysql_insertid'};
#	}
	
	# Processing features
	#my @cols = (qw/clone strand feature start stop 
	#	score orfstart orfstop flength/ );

	foreach my $array (@{$features->{$id}}) 
	{
	    #my $count = 0;
	    #my %featvals;
	    #$featvals{predicted_id} = $pred_id;
	    #$featvals{fgenesh_id} = $id;
	    #foreach my $col (@cols) {
	    #$featvals{$col} = $array->[$count];
	    #$count++;
	    # }
	    # Do the insert
	    #$sql->insert(-table=>'features',
	    #	   -field_vals=>\%featvals);
	    my @fvals = @$array;
	    my @JamieVals = @$array;
	    my $clone = shift @fvals;

	    #-----------------------------+
	    # ADD RELEVANT INFORMATION TO |
	    # THE BEGINNING OF THE ARRAY  |
	    #-----------------------------+
	    # We just set the pred_id to zero above, may want to move
	    # that down to a more useful location.
	    unshift @fvals, ($pred_id, $clone, $id);


	    if($printScreen)
	    {
		print "Insert features_fgenesh @fvals\n";
	    }
	    
	    #-----------------------------+
	    # PRINT THE GENE MODELS TO    |
	    # THE GFF FILE                |
	    #-----------------------------+
	    # The following prints all available variables
	    # This will be left here, but commented out
	    # for future reference. 
	    # 11/29/2006 - JCE
#	    print GFFOUT $fvals[0]."\t";             # PRED_ID
#	    print GFFOUT $fvals[1]."\t";             # CLONE
#	    print GFFOUT $fvals[2]."\t";             # ID
#	    print GFFOUT "STRN:".$fvals[3]."\t";     # STRAND
#	    print GFFOUT $fvals[4]."\t";             # 
#	    print GFFOUT "FEAT:".$fvals[5]."\t";     # FEATURE
#	    print GFFOUT "STRT:".$fvals[6]."\t";     # START
#	    print GFFOUT "STOP:".$fvals[7]."\t";     # STOP
#	    print GFFOUT "SCOR:".$fvals[8]."\t";     # SCORE
#	    print GFFOUT "OSTR:".$fvals[9]."\t";     # ORF-START 
#	    print GFFOUT "OSTP:".$fvals[10]."\t";    # ORF-STOP
#	    print GFFOUT "FLEN:".$fvals[11]."\t";    # F-LENGTH
#	    print GFFOUT "\n";                       #

	    #-----------------------------+
	    # PRINT THE GENE MODELS       |
	    #-----------------------------+
	    # Does not print the PolA results
	    unless ($fvals[5] =~ "PolA")
	    {
		#print GFFOUT $fvals[1]."\t";             # SEQNAME
		# Need to use the following for Apollo
		print GFFOUT "FGENESH-".$ref->{fgenesh_id}."-CDS\t";
		print GFFOUT "FGENESH\t";                # SOURCE
		print GFFOUT $fvals[5]."\t";             # FEATURE
		print GFFOUT $fvals[9]."\t";             # ORF-START 
		print GFFOUT $fvals[10]."\t";            # ORF-STOP
		print GFFOUT $fvals[8]."\t";             # SCORE
		print GFFOUT $fvals[3]."\t";             # STRAND
		print GFFOUT ".\t";                      # FRAME
		print GFFOUT "FGENESH-".$ref->{fgenesh_id}."-CDS\t";
		#print GFFOUT $fvals[4]."\t";             # ATTRIBUTE
		print GFFOUT "\n";
	    }
	    
	    

#	    if($insertDB)
#	    {
#		$feature_sth->execute(@fvals);
#	    }

	}
    }
    if($debug && $counter >= 1){
	die "That's it for test!\n";
    }
} # next file


close GFFOUT;


# Three tables:
# General clone info
# Predicted proteins
# unique id features

#-----------------------------------------------------------+
# NOTES                                                     |
#-----------------------------------------------------------+
# 11/28/2006
# - Many of the notes in the code above are not from J. Estill
# - It is successfully printing the output to the screen
# - Now I just need to get things in the proper format
#   for a gff dump.
# - J.Estill working on modifying the text parser for use
#   with the FGenesh program.
# - Fixed GFF output to just produce CDS information
#
# 10/30/2007
# - Attempting to modify the existing code to accept arguments
# - changed name to cnv_fgenesh2gff.pl from ParseFgenesh.pl

#-----------------------------------------------------------+
# TO DO
#-----------------------------------------------------------+
# 11/28/2006
# -It may be necessary to switch the start and stop numbers
#  to make sure that start is always less then stop.
# -Should allow for the ability to print the predicted
#  proteins from FGENESH to an output file.
