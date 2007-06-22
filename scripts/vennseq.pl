#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# vennseq.pl - DAWG-PAWS Venn Diagram Program               |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 03/06/2007                                       |
# UPDATED: 05/21/2007                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Creates a Venn Diagram comparing the overlap of sequence |
#  features along a query sequence. The purpose is to allow |
#  me to visualize the overlap of Repeat Databases or gene  |
#  models along a BAC.                                      |
#                                                           |
# DEPENDENCIES:                                             |
#  -BioPERL                                                 |
#  -NCBI BLAST                                              |
#  -Venn Master                                             |
#                                                           |
# USAGE:                                                    |
#                                                           |
#                                                           |
# NOTES:                                                    |
#  -All array indices start at one.                         |
#   This will facilitate parsing information from sequence  |
#   features that are all indexed from 1 to L               |
#  -q is quiet, Q is super quiet                            |
#-----------------------------------------------------------+
# To do. If no fasta file then the sequence length may
#        be passed as a variable from the command line.
=head1 NAME

vennseq.pl - DAWG-PAWS Venn Diagram Program

=head1 SYNOPSIS

    vennseq.pl -i InputFastaFile.fasta -o OutputFile -d FeatureFileDir
               -f tab
           
=head1 DESCRIPTION

Creates a Venn Diagram comparing the overlap of sequence
features along a query sequence. The purpose is to allow
me to visualize the overlap of Repeat Databases or gene
models along a BAC.            

=head1 ARGUMENTS

=over 2

=item -i InputFastaFile

=item -o OutputFile

=item -d FeatureFileDir

=item -f tab

The feature file format. This must be one of ( tab | blast | gff ).
Currently only tab delim blast hits are supported.

=back

=head1 REQUIREMENTS

Requires the VennMaster program 
L<http://www.informatik.uni-ulm.de/ni/staff/HKestler/vennm/doc.html>.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=cut

package DAWGPAWS::VennSeq;

print "The VennSeq Program has started.\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Std;               # Get options from the command line
use Bio::SearchIO;             # Parse BLAST output
use Bio::SeqIO;                # Parse fasta sequence input file
use GD;                        # Draw using the GD program
                               # In case I want to draw a coverage map

#-----------------------------+
# SET VARIABLE SCOPE          |
#-----------------------------+
my $Infile;                    # Path to sequence file in fasta format
my $Out;                       # Output file path for Venn text file
my $SvgOut;                      # Output file path for Venn SVG image output
my $FeatDir;                   # Directory containing sequence features
my $Format;                    # Format of the sequence feature files
                               #  - blast or b
                               #  - tab or t
                               #  - gff or g


my $SeqLen;                    # Full length of the sequence that the
                               # features are being mapped onto
my @FeatFiles;                 # Array to hold the path of the 
                               # feature files.
my @FeatFilesList;             # The original list of potential feature
                               # files. The zero length files will
                               # be ignored leaving only files with data
my @EmptyFiles;                # Files with no hits (0 length files)
my $NumFeatFiles;              # The number of feature files
my @FeatMatrix;                # Count of feature occurrences
my @BinFeatMatrix;             # Presence/abscence binary form of the
                               # feature matrix
my @FeatLabels;                # Category labels for the features
                               # this will take the 
my @CrossTab;                  # Cross table matrix
my $FileTestNum = 0;           # Initialize counter of file test number
my $EmptyFileNum = 0;          # Initialize counter of empty files

#-----------------------------------------------------------+
# COMMAND LINE VARIABLES                                    |
#-----------------------------------------------------------+
getopts('i:o:s:d:f:x:y:C:q', \%Options);
# x and y will possibly be used later to draw coverage map
# using x and y scale factor where x and y are integers
# greater then 1.
# may also use a "MAX" varialble
my $Usage = "\nVennSeq.pl -i InputFastaFile -o OutputFile\n".
	"-d FeatureFileDirecotry \n";
my $Config = $Options{C};
my $quiet = $Options{q};

if ($Config)
{
    print "Varialbes will be loaded from the config file";
    #-----------------------------+
    # GET CONFIG FILE OPTIONS     |
    #-----------------------------+
    $Infile = &ParseConfigFile($Config,"SeqIn")
        || die "Config file must include a fasta file\n$Usage";
    $Out =  &ParseConfigFile($Config,"Out")
        || die "Config file must include an output file\n$Usage";
    $SvgOut = &ParseConfigFile($Config, "VSvg")
	|| die "Config file must included an VSVG file path \n$Usage";
    $FeatDir = &ParseConfigFile($Config,"FeatDir")
        || die "Config file must include a sequence feature directory\n$Usage";
    $Format = &ParseConfigFile($Config,"Format")
        || die "Config file must indicate a feature file format\n$Usage";

}else{
    #-----------------------------+
    # GET COMMAND LINE OPTIONS    |
    #-----------------------------+
    $Infile = $Options{i} ||
	die "Command line must include  an input fasta file path.\n$Usage\n";
    $Out = $Options{o} ||
	die "Command line must include and output file path.\n$Usage\n";
    $FeatDir = $Options{d} ||
	die "Command line must include a sequence feature directory\n$Usage\n";
    $Format = $Options{f} ||
	die "Command line must include a feature file format\n$Usage\n";
    $SvgOut = $Options{s} ||
	die "Command line must include a svg file output path.\n$Usage\n";
}

#-----------------------------+
# CHECK FILE EXISTENCE        |
#-----------------------------+
unless (-e $Infile)
{
    print "ERROR:The input sequence file could not be found:\n$Infile\n";
}

#-----------------------------+
# OPEN THE OUTPUT FILE        |
#-----------------------------+
open (OUT, ">$Out") ||
    die "Could not open output file for writing at:\n$Out\n";

#-----------------------------+
# GET LENGTH OF THE SEQUENCE  |
# FROM THE INPUT FILE         |
#-----------------------------+
# This will be used to determine the length of the array
# This includes an internal check for only one sequence
# in the input file.
my $seq_in = Bio::SeqIO->new('-file' => "<$Infile",
                             '-format' => 'fasta') ||
    die "Can not open input file:n\$Infile\n";

my $NumSeq = 0;
while (my $inseq = $seq_in->next_seq)
{
    $NumSeq++;
    $SeqLen = $inseq->length();
    if ($NumSeq >> 1)
    {
	print "ERROR. Input sequence file contains more then one sequence.\n";
	exit;
    }
}

#-----------------------------------------------------------+
# GET SEQUENCE FEATURES                                     |
#-----------------------------------------------------------+
if ($Format =~ "tab")
{
    #-----------------------------------------------------------+
    # TAB DELIM BLAST OUTPUT                                    |
    #-----------------------------------------------------------+
    # Currently only supporting the blo extension
    opendir(INDIR, $FeatDir) ||
	die "Can not open input directory:\n$FeatDir: $!";
    @FeatFilesList = grep /blo$/, readdir INDIR;

    #-----------------------------+
    # TEST FOR FILE LENGTH        |
    #-----------------------------+

    # Check for files with no hits here by checking for zero length
    # 
    # Can load indexes to delete ??
    foreach my $FileTest (@FeatFilesList)
    {
	$FileTestNum++;
	
	#LOAD ZERO LENGTH FILE TO EmptyFiles array
	# Otherise push them to the FeatFiles array
	if (-z $FeatDir.$FileTest) # If the file is zero length
	{
	    $EmptyFileNum++;
#	    $EmptyFiles[$EmptyFileNum][1] = $EmptyFileNum;
#	    $EmptyFiles[$EmptyFileNum][2] = $FileTestNum;
#	    $EmptyFiles[$EmptyFileNum] = &GetShortName($FileTest);
#	    $EmptyFiles[$EmptyFileNum][3] = $FileTest;
	    # Load the short name of the database to the 
	    # empty files array, if the GetShort Name does
	    # not return a database then I will use the File Name
	    push (@EmptyFiles, &GetShortName($FileTest) || $FileTest);
	}else{
	    push (@FeatFiles, $FileTest); 
	}# End of if FileSize = 0
    } # End of for each file in File test
    
    #-----------------------------+
    # PRINT OUT THE EMPTY FILES   |
    #-----------------------------+
    # This is a place to drop empty files from the array
    if ($EmptyFileNum > 0)
    {
	print "The following files had no BLAST hits:\n";
	foreach my $IndEmptyFile (@EmptyFiles)
	{
	    print "\t".$IndEmptyFile."\n";
	}
    }
    
    $NumFeatFiles = @FeatFiles;
    print "\n$NumFeatFiles files have hits\n";

    # Temp exit while I test working with file size
#    exit;
    

    if ($NumFeatFiles == 0) 
    {
	print "No feature files of type \"$Format\" were found in:\n".
	    "$FeatDir\n";
	exit;
    }else{
	print "$NumFeatFiles of type $Format were found in\n".
	    "$FeatDir\n";
	# temp exit for debug
	#exit;
    }

    #-----------------------------+
    # FILL THE FEATURE MATRIX     |
    # WITH ZEROS                  |
    #-----------------------------+
    print "Initializing feature matrix and binary matrix with zeros\n";
    for (my $row=1; $row<=$NumFeatFiles; $row++)
    {
	my $NumCol=0;              # Varialbe to count number of cols
	                           # for debug purpsoses
	print "\tInitializing row $row...";
	for (my $col=1; $col<=$SeqLen; $col++)
	{
	    $FeatMatrix[$row][$col] = 0;
	    $BinFeatMatrix[$row][$col] = 0;
	    $NumCol++;
	}
	print "$NumCol cols.\n";
    }

    #-----------------------------+
    # FOR EACH OF THE FEATURE     |
    # FILES INCREMENT THE PROPER  |
    # ROW IN @FeatMatrix          |
    #-----------------------------+
    # Use uppercase Row and Col here to avoid
    # any problems with var names
    my $Row = 0;                   # Reset row count to zero
                                   # Each feature has a unique row

    print "Loading data from feature files.\n";
    foreach my $IndFile (@FeatFiles)
    {
	$Row++;                    # Increment row for each file

	$FeatLabels[$Row] = &GetShortName($IndFile);

	print "\tLoading data from feature file $FeatLabels[$Row]\n";

	my $FeatPath = $FeatDir.$IndFile;
	my $NumLines = 0;                  # Number of lines in the feat file
	my $BinFeatCount = 0;              # Counter for the "binary" features

	open (FEATIN, "<$FeatPath") ||
	    die "Can not open the feature file:\n$FeatPath\n";
	while (<FEATIN>)
	{
	    chomp;
	    $NumLines++;

	    # Use unless to ignore comment lines from -m 9 output
	    unless (m/^\#.*/)       
	    {
		# TEST THAT THE LENGTH OF THE ARRAY IS 
		# WHAT IS EXPECT
		my @TabSplit = split(/\t/);
		my $TestLen = @TabSplit;
		if ($TestLen > '12')
		{
		    print "ERROR: The BLAST file has an unexpected\n".
			"number of tab delimitedy columns.";
		    #exit;
		    # Currently just skip this record and move forward
		}

		my ($QryId, $SubId, $PID, $Len,
		    $MisMatch, $GapOpen,
		    $QStart,$QEnd, $SStart, $SEnd,
		    $EVal, $BitScore) = split(/\t/);
		
		#-----------------------------+
		# PRINT SOME FEATURE VALUES   |
		# FOR DEBUG PURPOSES          |
		#-----------------------------+
		# Just do this for the first few hits
		unless ($quiet)
		{
		    if ($NumLines < '3')
		    {
			print "\n\t\t  LINE: $NumLines\n";
			print "\t\t   QRY: $QryId\n";
			print "\t\t   SUB: $SubId\n";
			print "\t\tQSTART: $QStart\n";
			print "\t\t  QEND: $QEnd\n\n";
		    } # End of If NumLines less then three
		} # End of unless quiet

		for (my $Col=$QStart; $Col<=$QEnd; $Col++)
		{
		    $BinFeatCount++;
		    $FeatMatrix[$Row][$Col] = $FeatMatrix[$Row][$Col] + 1;
		} # End of for my $col from Query Start to End
	    } # End of unless # (unless comment line)
	} # End of while FEATIN
	close FEATIN;

	#-----------------------------+
	# REPORT NUMBER OF FEATURES   |
	# THAT WERE PROCESSED         |
	#-----------------------------+
	# For tab delimited blast the number of features
	# is roughly equal to the number of lines == number of HSPS
	#print "\t\tFILE: $FeatPath\n";
	unless ($quiet)
	{
	    print "\t\tFILE: $IndFile\n";
	    print "\t\tFEAT: $NumLines\n";
	    print "\t\t BIN: $BinFeatCount\n";
	} # End of unless quiet
    } # End of For each $IndFile in @FeatFiles

}else{
    print "A proper sequence feature file format was not selected.";
    exit;
} # END OF IF $Format

#-----------------------------------------------------------+
# FILL THE BINARY MATRIX AND PRINT OUTPUT FOR VENN MASTER   |
#-----------------------------------------------------------+

print "Filling the binary matrix.\n";
for ($row=1; $row<=$NumFeatFiles; $row++)
#for ($row=1; $row<=3; $row++)
{
    print "\tLoading row $row\n";

    my $BinFeatLen = 0;                    # Binary feature length
    my $BinFeatOccur = 0;                  # Occurrence of features

    for ($col=1; $col<=$SeqLen; $col++)
    {
	$BinFeatLen++;
	# If the feature matrix contains data
	if ($FeatMatrix[$row][$col] >> 0)
	{
	    $BinFeatOccur++;
	    $BinFeatMatrix[$row][$col] = 1;
	    # Print to screen first to see if this actually
	    # does work
	    print OUT "$col\t$FeatLabels[$row]\n";
#	    print "$col\t$FeatLabels[$row]\n";
	} 
    } # End of for $col

    #-----------------------------+
    # REPORT THE PERCENT COVERAGE |
    # OF THE FEATURE AFTER        |
    # TESTING FOR DIVIDE BY ZERO  |
    #-----------------------------+
    if ($BinFeatOccur == 0)
    {
	print "\t\tERROR: The feature has zero length.\n";
    } else {
	# Calculate feature coverage
	my $BinFeatCov = $BinFeatOccur/$BinFeatLen;
	# Convert feature coverage to percent
	$BinFeatCov = sprintf("%.4f",$BinFeatCov);
	print "\t\tLEN: $BinFeatLen\tCOV: $BinFeatCov\n"
    }

} # End of for row

# close the VENN Master text file to make it available
close OUT;

#-----------------------------------------------------------+
# GENERATE COMPARISON MATRIX @CrossTab                      |
#-----------------------------------------------------------+
# This should be made and option, ie. don't do under quiet
# Generates the NxN CrossTab matrix 
# Where N is the number of sequence features being comparedh zeros
print "Generating the Cross tab\n";

print "\n";

#-----------------------------+
# PRINT TABLE HEADER          |
#-----------------------------+


#-----------------------------+
# CREATE SEPERATOR STRING AND |
# PRINT TOP OF CROSS TABLE    |
#-----------------------------+
my $Sep = "+";
for (my $i=1; $i<=$NumFeatFiles ; $i++)
{$Sep = $Sep."--------+";}
print $Sep."\n";

#-----------------------------+
# PRINT BODY OF CROSS TABLE   |
#-----------------------------+
for (my $i=1; $i<=$NumFeatFiles ; $i++)
{
    print "|"; # Print the left hand side of the cross tab

    for (my $j=1; $j<=$NumFeatFiles ; $j++)
    {
	my $SharedCount=0;    # Initialize the count of shared
	my $iCount=0;         # Initialize the i feature count
	my $CrossRatio; # Set scope of the cross ratio

	if ($i != $j)
	{

	    #-----------------------------+
	    # CALCULATE CROSS TABLE RATIO |
	    #-----------------------------+
	    # k is the sequence position
	    for (my $k=1; $k<=$SeqLen; $k++)
	    {
		if ($BinFeatMatrix[$i][$k] == 1)
		{
		    # Increment the i feature count
		    $iCount++;
		    if ($BinFeatMatrix[$j][$k] == 1)
		    {
			# Increment the shared feature count
			$SharedCount++;
		    } # End of if j = 1 in binary matrix
		} # End of if i = 1 in binary matrix
	    } # End of for k
	    
	    #-----------------------------+
	    # LOAD RESULT TO CROSS TABLE  |
	    # AFTER DIVIDE BY ZERO TEST   |
	    #-----------------------------+
	    if ($iCount >0)
	    {
		$CrossRatio = $SharedCount/$iCount;
		$CrossRatio = sprintf("%.4f",$CrossRatio);
	    }else{
		# If if = 0 this is the NULL case
		$CrossRatio = "  NULL  ";
	    } # End of if $iCount > 0
	    
	    # Load to the Cross Tab Array
	    $CrossTab[$i][$j]=$CrossRatio;	    
	    # To print in table format
	    print " ".$CrossRatio." |";
	}else{
	    # If i = j then print the data source
	    # label in the table, this will need to
	    # be formatted to fit in the box
	    # - sign below indicates left justification
	    my $DatSrcLab = sprintf ('%-*s', 8, $FeatLabels[$i]);
	    print $DatSrcLab."|";
	}
	
    } # End of for j (Seq feature set two)

    print "\n"; # Print new line for end of row in cross tab
    print $Sep."\n";
} # End of for i (Seq feature set one)

#print $Sep."\n";
print "\n"; # Print final new row at the very end of cross tab

#-----------------------------------------------------------+
# OPEN VENN MASTER TO DISPLAY VENN DIAGRAMS                 |
#-----------------------------------------------------------+
print "Running the VennMaster program.\n";

my $VMasterDir = "/home/jestill/Apps/VennMaster/VennMaster-0.33.6/";
my $JavaPath = "/usr/java/jre1.6.0/bin/java";

# Working path is:
# /usr/java/jre1.6.0/bin/java -Xms256m -Xmx256m -jar 
#   /home/jestill/Apps/VennMaster/VennMaster-0.33.6/venn.jar

# Without the At we have
my $VCmd = "$JavaPath -Xms256m -Xmx256m -jar $VMasterDir".
    "venn.jar --list $Out --svg $SvgOut";

system ($VCmd) ||
    die "Could not run the VennMaster program with cmd:\n$VCmd\n";

#-----------------------------------------------------------+
# CREATE COVERAGE DIAGRAM AS A HEAT MAP                     |
#-----------------------------------------------------------+

# This can me made to be an option

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub ParseConfigFile
{

#-----------------------------------------------------------+
# This will parase a configuration text file to set         |
# variables that would normally be set at the command line  |
# This is not the fastest way to parse the config file      |
# but this makes the subfuction reusable.                   |
# If the variable name occurs more then once in the text    |
# file, the last occurrence will be used.                   |
#-----------------------------------------------------------+
    my $ConfigFilePath = $_[0];
    my $VarName = $_[1];
    my $VarValue;


    open (CONFILE, $ConfigFilePath) ||
        die "Could not open config file:\n\t$ConfigFilePath";

    while (<CONFILE>)
    {
        chomp;                 # Remove newline character
        unless (m/\#.*/)       # Ignore comment lines
        {
            my @SplitLine = split;
            if ($SplitLine[0] =~ $VarName){
                $VarValue = $SplitLine[1];}
        }
    }
    close CONFILE;
    return $VarValue;

}

sub GetShortName
{
#-----------------------------------------------------------+
# CONVERT THE  NAME TO A SHORTENED FORM                     |
#-----------------------------------------------------------+ 
# This will be used to write short names for the CrossTab
# table or any other use of the name in the program. The short
# form of the name will be limited to 8 characters
# NOTE:
# This is a kluge that will only work for the BLAST databases
# that I am working with in wheat (03/12/2007 - JCE)
# I should definited convert this to as hash table

    my $InName = $_[0]; #  Input name string
    my $OutName;        #  Set scope for the name to return


    #-----------------------------------------------------------+
    # REPEAT DATABASES                                          |
    #-----------------------------------------------------------+
    if ($InName =~ m/TREP9_nr/)
    {
	$OutName = "TREP9nr";
    }elsif ($InName =~ m/TREP9_total/){
	$OutName = "TREP9tot";
    }elsif ($InName =~ m/mips_REdat/){
	$OutName = "MIPS";
    }elsif ($InName =~ m/TIGR_Gram/){
	$OutName = "TIGRGram";
    }elsif ($InName =~ m/RB_pln/){
	$OutName = "RB_pln";
    }elsif ($InName =~ m/TATrans/){
	$OutName = "TATran";
    }elsif ($InName =~ m/Wessler/){
	$OutName = "Wessler ";
    }elsif ($InName =~ m/SanMiguel/){
	$OutName = "SanMig";
    #-----------------------------------------------------------+
    # EST DATABASES                                             |
    #-----------------------------------------------------------+
    # TIGR GENE INDICES	
    }elsif ($InName =~ m/TaGI_10/){
	$OutName = "TaGI_10";    	
    }elsif ($InName =~ m/AtGI_13/){
	$OutName = "AtGI_13";    	
    }elsif ($InName =~ m/ZmGI_17/){
	$OutName = "ZmGI_17";
    }elsif ($InName =~ m/SbGI_8/){
	$OutName = "SbGI_8";
    }elsif ($InName =~ m/OsGI_17/){
	$OutName = "OsGI_17";
    }elsif ($InName =~ m/HvGI_9/){
	$OutName = "HvGI_17";
    #-----------------------------+ 
    # EST: TIGR TAs               |
    #-----------------------------+
    }elsif ($InName =~ m/AcTA_1/){
	$OutName = "AcTA_1";
    }elsif ($InName =~ m/AsTA_1/){
	$OutName = "AsTA_1";
    }elsif ($InName =~ m/AsTA_1/){
	$OutName = "AsTA_1";
    }elsif ($InName =~ m/AvenSatTA_2/){
	$OutName = "AvenSaTa";
    }elsif ($InName =~ m/CdTA_2/){
	$OutName = "CeTA_2";
    }elsif ($InName =~ m/FaTA_2/){
	$OutName = "FaTA_2";
    }elsif ($InName =~ m/HvTA_2/){
	$OutName = "HvTA_2";
    }elsif ($InName =~ m/OsTA_2/){
	$OutName = "OsTA_2";
    }elsif ($InName =~ m/PgTA_2/){
	$OutName = "PgTA_2";
    }elsif ($InName =~ m/SbTA_2/){
	$OutName = "SbTA_2";
    }elsif ($InName =~ m/ShTA_2/){
	$OutName = "ShTA_2";
    }elsif ($InName =~ m/SoTA_2/){
	$OutName = "SoTA_2";
    }elsif ($InName =~ m/SpTA_2/){
	$OutName = "SpTA_2";
    }elsif ($InName =~ m/TaTA_2/){
	$OutName = "TaTA_2";
    }elsif ($InName =~ m/TmTA_2/){
	$OutName = "TmTA_2";
    }elsif ($InName =~ m/TtTA_1/){
	$OutName = "TtTA_1";
    }elsif ($InName =~ m/ZmB73TA_1/){
	$OutName = "ZmB71TA_1";
    #-----------------------------+	
    # EST: NCBI UniGenes	  |
    #-----------------------------+
    }elsif ($InName =~ m/AtUniGene_56/){
	$OutName = "AtUn_56";
    }elsif ($InName =~ m/HvUniGene_47/){
	$OutName = "HvUn_47";
    }elsif ($InName =~ m/OsUniGene_65/){
	$OutName = "OsUn_65";
    }elsif ($InName =~ m/SbUniGene_21/){
	$OutName = "SbUn_21";
    }elsif ($InName =~ m/SoUniGene_8/){
	$OutName = "SoUn_8";
    }elsif ($InName =~ m/TaUniGene_46/){
	$OutName = "TaUn_46";
    }elsif ($InName =~ m/ZmUniGene_61/){
	$OutName = "ZmUn_61";
    #-----------------------------+
    # EST: PLANT GDB PUTS         |
    #-----------------------------+
    }elsif ($InName =~ m/AsPUT_157/){
	$OutName = "AsPUT157";
    }elsif ($InName =~ m/AtPUT_157/){
	$OutName = "AtPUT157";
    }elsif ($InName =~ m/HvPUT_157/){
	$OutName = "HvPUT157";
    }elsif ($InName =~ m/OsIndPUT_157/){
	$OutName = "OsIndPUT";
    }elsif ($InName =~ m/OsJap_157/){
	$OutName = "OsJapPUT";
    }elsif ($InName =~ m/OsPut_157/){
	$OutName = "OsPUT157";
    }elsif ($InName =~ m/SbPUT_157/){
	$OutName = "SbPUT157";
    }elsif ($InName =~ m/SpPUT_157/){
	$OutName = "SpPUT157";
    }elsif ($InName =~ m/TaPUT_157/){
	$OutName = "TaPUT157";
    }elsif ($InName =~ m/TmPUT_157/){
	$OutName = "TmPUT157";
    }elsif ($InName =~ m/ZmPUT_157/){
	$OutName = "ZmPUT157";
    #-----------------------------------------------------------+
    # PLANT PROTEIN DATABASES                                   |
    #-----------------------------------------------------------+
    #-----------------------------------------------------------+
    # FINISHED GENOMES PROTEINS BLASTX                          |
    #-----------------------------------------------------------+
    }elsif ($InName =~ m/TAIR6/){
	$OutName = "TAIR6";
    }elsif ($InName =~ m/Os_5_cds_pep/){
	$OutName = "OS_5_cds_pep";
    }elsif ($InName =~ m/Os_5_RAP1_pep/){
	$OutName = "OS_5_cds_pep";
    #-----------------------------------------------------------+
    # UNIPROT                                                   |
    #-----------------------------------------------------------+
    }elsif ($InName =~ m/Os_5_RAP1_pep/){
	$OutName = "OS_5_cds_pep";
    #-----------------------------------------------------------+
    # COMBINATION DATABASES                                     |
    #-----------------------------------------------------------+
    }elsif ($InName =~ m/EST_DB/){
	$OutName = "EST_DB";
    }elsif ($InName =~ m/TE_DB/){
	$OutName = "TE_DB";
    #-----------------------------------------------------------+
    # UNKNOWN DATABASES                                         |
    #-----------------------------------------------------------+
    # Use the input name
    # It may be better to flag this as UNK
    }else{

	$OutName = "UNK";
	#$OutName = $InName;
    }

    return $OutName;

}

=head1 HISTORY

Started: 03/06/2007

Updated: 05/21/2007

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 03/06/2007
# - Program started
# - Included ParseConfigFile from previous work
# - Main body of program written
#   Reading input fasta file and sequence features and
#   loading to the FeatMatrix 2d array.
# - Printing output to screen for testing purposes
# - Ability to read tab delim blast output
#
# 03/10/2007
# - Working with test tabl delim blast output
# - Adding a coverage calculation and report
# - Adding quiet flag to turn of debug print information
# - Fixed problem parsing the tab delimited BLAST files:
#   I was using The regular expression: (m/\#.*/) to ignore
#   comment lines. However this will ignore any line that 
#   CONTAINS the # symbol. Changing the reg exp to
#   (m/^\#.*/) fixes the problem. This is a particular issue with
#   BLAST databases that have been formatted in the 
#   RepBase format. (ie EnSpm1_TD#EnSpm)
# - Added CrosTab function
#
# 03/12/2007
# - Added GetShortName subfunction
# - Got GetShortName working for Repeat databases and
#   the for all other BLAST databases used for the 
#   Wheat Annotation project 
# - Added the EmptyFiles array. This will prevent
#   wasting time doing counts for zero length files
# - Added print of the Empty files to the end of the array
#
# 03/14/2007
# - Working on the automated startup of Venn Master
#
# 05/21/2007
# - Added Pod Documentation


#-----------------------------------------------------------+
# TO DO LIST                                                |
#-----------------------------------------------------------+
# - Use GD to draw image of feature coverage and richness

# - TEST CMD IS:
#   ./vennseq.pl -i /home/jestill/projects/wheat_annotation/VennTest/HEX0075G21.fasta.masked -o /home/jestill/projects/wheat_annotation/VennTest/TestOutTwo.txt -d /home/jestill/projects/wheat_annotation/VennTest/ -f tab -s /home/jestill/projects/wheat_annotation/VennTest/HEX0075G21.test.svg


