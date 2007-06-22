#!/usr/bin/perl -w
#-----------------------------------------------------------+
#
# GFF2Game.pl
#
#-----------------------------------------------------------+
#
#  AUTHOR: James C. Estill
# STARTED: 02/09/2007
# UPDATED: 02/12/2007
# DESCRIPTION:
#  Converts gff data tracks from gff format to the game
#  xml format for use in the Apollo Genome Annotation 
#  Curation program. This is for use with output from the
#  TriAnnot pipeline.
# 
# DEPENDENCIES:
# *Apollo
#    Requires that apollo be installed on the local machine
#    since apollo is being used as the engine to do the
#    conversion between formats.
#-----------------------------------------------------------+ 

# TO DO:
# The following three should be variables that can be passed
# at the command line if desired.


# THE FOLLOWING ARE FOR TESTING ON THE JLB10 MACHINE
#my $InFastaFile = "/home/jestill/projects/wheat_annotation/wed/Wed/".
#    "DQ537336/DQ537336.fasta";

#my $InputRoot = "/home/jestill/projects/wheat_annotation/wed/Wed/".
#    "DQ537336/JamesEstill17103038/";

#my $OutputRoot = "/home/jestill/projects/wheat_annotation/wed/Wed/".
#    "DQ537336/game/";

# THE FOLLOWING IS FOR TESTING ON THE ATEAM MAC MACHINE
#my $InFastaFile = "/Users/ateam/Wed/DQ537336/DQ537336.fasta";
#my $InputRoot = "/Users/ateam/Wed/DQ537336/gffs/";
#my $OutputRoot = "/Users/ateam/Wed/DQ537336/game/";


# THE FOLLOWING IS FOR TESTING ON THE ATEAM MAC MACHINE
#my $InFastaFile = "/Users/ateam/Wed/DQ537335/DQ537335.fasta";
#my $InputRoot = "/Users/ateam/Wed/DQ537335/gff/";
#my $OutputRoot = "/Users/ateam/Wed/DQ537335/game/";

# NEEDS GFF FILE AND GAME DIR FOR OUTPUT
my $InFastaFile = "/Users/ateam/DirtyDozen/HEX1057J19/HEX1057J19.fasta.masked";
my $InputRoot = "/Users/ateam/DirtyDozen/HEX1057J19/gffs/";
my $OutputRoot = "/Users/ateam/DirtyDozen/HEX1057J19/game/";

#-----------------------------+ 
# CHECK FOR EXISTENCE OF      |
# INPUT FILES                 |
#-----------------------------+
unless (-e $InFastaFile)
{
    "The fasta file does not exit\n:".$InFastaFile."\n";
    exit;
}


# THE FOLLOWING ARE TYPICAL FOR THE
# TRIANNOT PIPELINE
#MY $FILEROOT = "1blastnRpb";       # RepBase

#
# TO DO:
# * Load this array using the command to list the gff files in the
#   directory
# * Also use the ls command to get the fasta file to use with this
# 

my @Files2Convert = (
#		     "1blastnTIGRRep",    # TIGR Repeat Database
#		     "1blastnTrp",        # Dont Know
#		     "1blastnTtot",       # TREP total
#		     "1blastxTpr",        # TREP protein
# Parsers for the files above are not available
# TEST OF SINGLE GFF
# This would need to change how the software is aware of what parser
# to use and will probably not be worth the effort given that I
# want to quit using TriAnnot Pipeline
#		     "gffcat"
# START OF WORKING
		     "1trf",              # Trandem Repeat Finder
		     "2fGh",              # Fgenesh - Gene Model Prediction
		     "2gID",              # GeneID - Gene Model Prediction
		     "2gmHv",             # GeneMark Hordeum vulgare
		     "2gmOs",             # GeneMark Oryza sative
		     "2gmTa",             # GeneMark Triticum aestevum
		     "2gmZm",             # GeneMark Zea mays 
# END OF WORKING
# Parsers for the files below are not available
#		     "3blastnOsU",        # Blastn Oryza sativa Unigenes
#		     "3blastnTaU",        # Blastn Triticum aestivm Unigenes
#		     "3blastnZmU",        # Blastn Zea mays Unigenes
#		     "3blastxSprot",      # Blastx
#		     "3blastxTrembl",     # Blastx
# 		     "3tblastxTAta"       # tBLASTx
		     );


foreach my $FileRoot (@Files2Convert)
{   

    my $TriGffPath = $InputRoot.$FileRoot.".gff";      # TriAnnot Gff
    my $ApGffPath = $InputRoot.$FileRoot.".ap.gff";     # Apollo Gff
    my $OutPath = $OutputRoot.$FileRoot.".game.xml";   # Game XML output 
    
    #-----------------------------+
    # CONVERT TRIANNOT GFF TO     |
    # APOLLO USABLE GFF           |
    #-----------------------------+
    # This will need to be a different convert for each
    # of the input files

    #---------------+
    # FGENESH       |
    #---------------+
    if ($FileRoot =~ "2fGh")
    {
	&ModFGENESH ( $TriGffPath, $ApGffPath );
    }
    
    #---------------+
    # GeneMark.hmm  |
    # Ta Matrix     |
    #---------------+
    elsif ($FileRoot =~ "2gmTa")
    {
	&ModGenMark ( $TriGffPath, $ApGffPath );
    }

    #---------------+
    # GeneMark.hmm  |
    # Os Matrix     |
    #---------------+
    elsif ($FileRoot =~ "2gmOs")
    {
	&ModGenMark ( $TriGffPath, $ApGffPath );
    }

    #---------------+
    # GeneMark.hmm  |
    # Hv Matrix     |
    #---------------+
    elsif ($FileRoot =~ "2gmHv")
    {
	&ModGenMark ( $TriGffPath, $ApGffPath );
    }

    #---------------+
    # GeneMark.hmm  |
    # Zm Matrix     |
    #---------------+
    elsif ($FileRoot =~ "2gmZm")
    {
	&ModGenMark ( $TriGffPath, $ApGffPath );
    }

    #---------------+
    # GeneId        |
    #---------------+
    elsif ($FileRoot =~ "2gID")
    {
	&ModGeneId ( $TriGffPath, $ApGffPath );
    }

    #---------------+
    # TRF           |
    #---------------+
    # Tandem Repeat Finder
    elsif ($FileRoot =~ "1trf")
    {
	&ModTRF ( $TriGffPath, $ApGffPath );
    }

    #-----------------------------+
    # Convert each gff file game  |
    # xml using the ApolloConvert |
    # subfunction                 |
    #-----------------------------+
    # Will comment this out while I am working
    # on the FGENESH and GenMark conversion subfunctions
    &ApolloConvert ($ApGffPath, "gff", $OutPath, 
		    "game", $InFastaFile, "NULL");

}

# Add bioperl dependent code here to concatenate the game xml files
# to a single XML file name for the root name in the fasta file or
# a root name given by the user.

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub ApolloConvert
{
#-----------------------------+
# CONVERT AMONG FILE FORMATS  |
# USING THE APOLLO PROGRAM    |
#-----------------------------+
# Converts among the various data formats that can be used 
# from the command line in tbe Apollo program. For example
# can convert GFF format files into the game XML format.
# NOTES:
#  - Currently assumes that the input file is in the correct
#    coordinate system.
#  - GFF files will require a sequence file
#  - ChadoDB format will require a db password


    # ApPath - the path of dir with the Apollo binary
    #          Specifying the path will allow for cases
    #          where the program is not in the PATHS
    # ApCmd  - the apollo commands to run

    my $InFile = $_[0];        # Input file path
    my $InForm = $_[1];        # Output file format:
                               # game|gff|gb|chadoxml|backup
    my $OutFile = $_[2];       # Output file path
    my $OutForm = $_[3];       # Ouput file foramt
                               # chadoDB|game|chadoxml|genbank|gff|backup
    my $SeqFile = $_[4];       # The path of the sequence file
                               # This is only required for GFF foramt files
                               # When not required this can be passed as na
    my $DbPass = $_[5];        # Database password for logging on to the 
                               # chado database for reading or writing.
    my ( $ApPath, $ApCmd );


# The following path is for the jlb10 machine
#    $ApPath = "/home/jestill/Apps/Apollo/Apollo";

    $ApPath = "/Applications/Apollo/bin/apollo";

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
    $ApCmd = $ApPath." -i ".$InForm." -f ".$InFile.
	" -o ".$OutForm." -w ".$OutFile;

    # Make sure that that input output formats are in lowercase
    # may need to add something here to avoid converting chadoDB
    $InForm = lc($InForm);
    $OutForm = lc($OutForm);
    
    # Determine the proper command to use based on the input format
    # since GFF file also require a sequence file
    if ($InForm =~ "gff" )
    {
	$ApCmd = $ApCmd." -s ".$SeqFile;
    }
    
    if ($InForm =~ "chadodb")
    {
	$ApCmd = $ApCmd." -D ".$DbPass;
    }

    # Do the apollo command
    system ( $ApCmd );

}


sub ModFGENESH
{
#-----------------------------+
# MODIFY THE FGENESH GFF FILE |
# TO WORK WITH APOLLO         |
#-----------------------------+
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debu
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$Infile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "exon")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/_/, $Sp[8]);

	    # Print to screen for debug
	    print $LineCount."\t".$Sp[2]."\n";
	    print "\t".$Id[1]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    print OUT $Id[1]."\n";

	}
	

    }

    close IN;
    close OUT;

}



sub ModGenMark
{
#-----------------------------+
# MODIFY THE FGENESH GFF FILE |
# TO WORK WITH APOLLO         |
#-----------------------------+
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debu
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$Infile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "exon")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/;/, $Sp[8]);
	    my @ModId = split(/_/, $Id[0]);
	    # Print to screen for debug
	    #print "\t".$Id[0]."\n";
	    print $LineCount."\t".$ModId[2]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    print OUT $ModId[2]."\n";

	}
	

    }

    close IN;
    close OUT;

}

sub ModGeneId
{
#-----------------------------+
# MODIFY THE FGENESH GFF FILE |
# TO WORK WITH APOLLO         |
#-----------------------------+
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debu
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$Infile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "exon")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/_/, $Sp[8]);

	    # Print to screen for debug
	    print $LineCount."\t".$Sp[2]."\n";
	    print "\t".$Id[1]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    print OUT $Id[1]."\n";

	}
	

    }

    close IN;
    close OUT;

}



sub ModTRF
{
#-----------------------------+
# MODIFY THE TRF GFF FILE     |
# TO WORK WITH APOLLO         |
#-----------------------------+
# TRF: Tandem Repeat Finder
# InFile is the full path to a valid gff infile
# OutFile is the full path for the modified file

    my $InFile = $_[0];
    my $OutFile = $_[1];
    my $LineCount = 0;  # LineCount used for debug
    
    open (IN, "<$InFile") ||
	die "Can not open infile:\n$Infile\n";
    open (OUT, ">$OutFile") ||
	die "Can not open outfile:\n$OutFile\n";
    
    while (<IN>)
    {
	$LineCount++;
	chomp;                  # Remove the newline char
	my @Sp = split(/\t/);   # Split by the tab character

	if ($Sp[2] =~ "tandem_repeat")
	{

	    # Split out the name to use for the
	    # exon
	    my @Id = split(/;/, $Sp[8]);
	    # SmID for Small ID, Additional information
	    # can be parsed from the TriAnnot pipeline GFF file
	    # including the length and repeat motif
	    my @SmId = split(/_/, $Id[0]);

	    # Print to screen for debug
	    print $LineCount."\t".$Sp[2]."\n";
	    print "\tTRF".$SmId[2]."\n";

	    # PRINT OUTPUT TO GFF FILE
	    print OUT $Sp[0]."\t";
	    print OUT $Sp[1]."\t";
	    print OUT $Sp[2]."\t";
	    print OUT $Sp[3]."\t";
	    print OUT $Sp[4]."\t";
	    print OUT $Sp[5]."\t";
	    print OUT $Sp[6]."\t";
	    print OUT $Sp[7]."\t";
	    # Small Id is just TRF####
	    print OUT "TRF".$SmId[2]."\n";
	}
	

    }

    close IN;
    close OUT;

} # End of ModTRF



=head1
HISTORY
=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 
# 02/09/2007
#  - Program started basic funnctionality of main body
#  - Added the ApolloConvert subfunction from pre-existing
#    source
#
# 02/12/2007
#  - Added ModFGENESH subfunction to convert the TriAnnot
#    GFF file to a format that Apollo can use.
#  - Added ModGenmark subfunction to convert the GenMark
#    files to the format Apollo can ues.
#  - Added code to the body of the script to use the
#    ModFGENESH and ModGenmark subfunctions
#
# 02/13/2007
#  - Added ModTRF subfunction. The names are being
#    lost when this GFF file is parsed, but leaving
#    this alone for now.
#
# 02/19/2007
# - Working on a modification to read fasta files
#   and load to a Files to convert. This will need
#   to also be aware of where the gff files are
#   for the fasta file of interest.
#
