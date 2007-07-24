#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# BLAST PAWS                                                |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 05/10/2006                                       |
# UPDATED: 04/24/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run NCBI blast for the wheat BACs against a set of       |
#  repeast databases and EST contigs of interest.           |
#                                                           |
# USAGE:                                                    |
#  Blast_PAWS.pl -d DatabaseSet -t                          |
#   -d Database option. Valid values are limited to:        |
#       - embl --> EMBL (UniProt/Swiss-Prot) BLASTX         |
#       - estn --> EST contig databases BLASTN              |
#       - estx --> EST contig databases TBLASTX             |
#       - prot --> Protein databases                        |
#       - repn --> Repetitive element databases blastn      |
#       - repx --> Repetitive element databases blastx      |
#       - genx --> Finished Genomes proteins (CDS from      |
#                  arabidopsis,rice ..) BLASTX              |
#       - cont --> Potenetial contamination. BLASTN         |
#                  E.coli, vector (human ?)                 |
#   -t Test (optional). This will Test that all input files |
#      and will create the directories if needed.           |
#                                                           |
#-----------------------------------------------------------+

# TODO:
#  - Add output dir option that is different then then input
#    dir path
print "The Blast_PAWS program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Std;               # Allows to get options from the command line
use Text::Wrap;                # Allows word wrapping and hanging indents
                               # for more readable output for long strings.


#-----------------------------+
# FILE OPTIONS                |
#-----------------------------+
my $QDir;                     # Base dir for qry sequences  
my $DbDir;                    # Base dir for BLAST databases
#$QDir = "/home/jestill/projects/wheat_annotation/200702/DirtyDozen/";
$DbDir = "/home/jestill/blast/paws/";

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my @Db;                       # Array of the databases to search
my $BlSuf;                    # Blast suffix
my $BlProg;                   # Blast program to use 
my $QryRoot;                  # Root name of the QRY fasta file
my $OutDirPath;               # Full path of the output directory
my $OutExt;                   # The extension to add to the output file
my $MakeGff;                  # Create a GFF output file
my ( $BlastCmd, $QryPath, $DbPath, $OutPath ); 
my $ProcNum = 0;

my $usage = "  Blast_PAWS.pl -d DatabaseSet -t -g\b".
    "-i Directory containing fasta files for the BLAST query.\n".
    "-d Database option. Valid values are limited to:\n".
    "   -d embl--> EMBL (TrEmbl/Swiss-Prot) databases blastx\n".
    "   -d estn ---> EST contig databases blastn\n".
    "   -d estx ---> EST contig databases tblastx\n".
    "   -d prot --> Protein databases\n".
    "   -d repn --> Repetitive element databases blastn\n".
    "   -d repx --> Repetitive element protein database blastx\n".
    "   -d cont --> Potential contamination\n".
    "   -d genx --> BLASTX against finished genome CDS models\n".
    "-g Create gff output BLAST".
    "-t Test (optional). This will Test that all input files\n".
    "   and output dirs exist. The output directories will be\n ".
    "   created if needed.\n";

#-----------------------------+
# COMMAND LINE VARIABLES      |
#-----------------------------+
# Options hash to hold the options from the command line
my %Options;                  
getopts('i:d:q:tg', \%Options);      # Get the options from the command line
my $test;
my $DbSet = $Options{d} ||      # Database set to use
    die "ERROR. A database set must be identified.\n$usage";
$test = $Options{t};
$QDir = $Options{i};
$MakeGff = $Options{g};

#-----------------------------+
# LOAD THE LIST OF INPUT      |
# FILES TO THE @Qry ARRAY     |
#-----------------------------+
# Assumes that all files have been masked and are
# named *.fasta.masked
opendir(IN, $QDir) || 
    die "Can not open input directory:\n$QDir: $!";
my @Qry = grep /fasta\.masked$/, readdir IN;
closedir IN;
  
#-----------------------------+
# BLAST QUERY DATABASES       |
#-----------------------------+
if ($DbSet =~ "estn")
{

    #---------------------+
    # BLAST OPTIONS       |
    #---------------------+
    $BlProg = "blastn";
    # -m 8 ---> Produce tab delimited BLAST output
    # Added the m8 option
    $BlSuf = "-e 0.00001 -a 2 -U -m8";
#    $BlSuf = "-e 0.00001 -a 2 -U";
    $OutExt = "blo";
    
    @Db = ( 
	    #---------------------+
	    # TIGR GENE INDICES   |
	    #---------------------+
	    "TaGI_10",
	    "AtGI_13",
	    "ZmGI_17",
	    "SbGI_8",
	    "OsGI_17",
	    "HvGI_9",
	    #---------------------+
	    # TIGR TAs            |
	    #---------------------+
	    "AcTA_1",
	    "AsTA_1",
	    "AtTA_2",
	    "AvenSatTA_2",
	    "CdTA_2",
	    "FaTA_2",
	    "HvTA_2",
	    "OsTA_2",
	    "PgTA_2",
	    "SbTA_2",
	    "ShTA_2",
	    "SoTA_2",
	    "SpTA_2",
	    "TaTA_2",
	    "TmTA_2",
	    "TtTA_1",
	    "ZmB73TA_1",
	    "ZmTA_3",
	    #---------------------+
	    # NCBI UniGenes       |
	    #---------------------+
	    "AtUniGene_56",        # Arabidopsis thaliana
	    "HvUniGene_47",        # Hordeum vulgare
	    "OsUniGene_65",        # Oryza sativa
	    "SbUniGene_21",        # Sorghu bicolor
	    "SoUniGene_8",         # Saccharum
	    "TaUniGene_46",        # Triticum aestivum
	    "ZmUniGene_61",        # Zea mays
	    #---------------------+
	    # Plant GDB PUTS      |
	    #---------------------+
	    "AsPUT_157",
	    "AtPUT_157",
	    "HvPUT_157",
	    "OsIndPUT_157",
	    "OsJapPUT_157",
	    "OsPUT_157",
	    "SbPUT_157",
	    "SpPUT_157",
	    "TaPUT_157",
	    "TmPUT_157",
	    "ZmPUT_157"
	    );
    
}
#-----------------------------+
# TBLASTX TO EST DATABASES    |
#-----------------------------+
elsif ($DbSet =~ "estx")
{
    #---------------------+
    # BLAST OPTIONS       |
    #---------------------+
    # TO DO Changed the sequence CDS
    # to an actual protein database
    $BlProg = "tblastx";
#    $BlSuf = "-e 0.00001 -U";
    $BlSuf = "-e 0.00001 -U -m 8";
    $OutExt = "blx";

    @Db = ( 
	    #---------------------+
	    # TIGR GENE INDICES   |
	    #---------------------+
	    "TaGI_10",
	    "AtGI_13",
	    "ZmGI_17",
	    "SbGI_8",
	    "OsGI_17",
	    "HvGI_9",
	    #---------------------+
	    # TIGR TAs            |
	    #---------------------+
	    "AcTA_1",
	    "AsTA_1",
	    "AtTA_2",
	    "AvenSatTA_2",
	    "CdTA_2",
	    "FaTA_2",
	    "HvTA_2",
	    "OsTA_2",
	    "PgTA_2",
	    "SbTA_2",
	    "ShTA_2",
	    "SoTA_2",
	    "SpTA_2",
	    "TaTA_2",
	    "TmTA_2",
	    "TtTA_1",
	    "ZmB73TA_1",
	    "ZmTA_3",
	    #---------------------+
	    # NCBI UniGenes       |
	    #---------------------+
	    "AtUniGene_56",        # Arabidopsis thaliana
	    "HvUniGene_47",        # Hordeum vulgare
	    "OsUniGene_65",        # Oryza sativa
	    "SbUniGene_21",        # Sorghu bicolor
	    "SoUniGene_8",         # Saccharum
	    "TaUniGene_46",        # Triticum aestivum
	    "ZmUniGene_61",        # Zea mays
	    #---------------------+
	    # Plant GDB PUTS      |
	    #---------------------+
	    "AsPUT_157",
	    "AtPUT_157",
	    "HvPUT_157",
	    "OsIndPUT_157",
	    "OsJapPUT_157",
	    "OsPUT_157",
	    "SbPUT_157",
	    "SpPUT_157",
	    "TaPUT_157",
	    "TmPUT_157",
	    "ZmPUT_157"
	    );

}

# FINISHED GENOMES PROTEINS BLASTX
elsif ($DbSet =~ "genx")
{
    $BlProg = "blastx";
#    $BlSuf = "-e 0.00001 -U";
    $BlSuf = "-e 0.00001 -U -m 8";
    $OutExt = "blx";
    
    @Db = (
	   "TAIR6",             # Arabidopsis thaliana (TAIR)
	   "Os_5_cds_pep",      # Oryza sativa TIGR v 5. (lots of TEs) 
	   "Os_RAP1_pep"        # O.sativa japonica, better protein set
	   );
}

# POTENTIAL CONTAMINATION
elsif ($DbSet =~ "cont")
{
    # CURRENTLY THESE ARE NOT MASKED
    $BlProg = "blastn";
#    $BlSuf = "-e 0.00001";
    $BlSuf = "-e 0.00001 -m 8";
    $OutExt = "blo";
    
    @Db = (
	   "contam_Ecoli",         # E. coli
	   "contam_UniVec"         # NCBI UniVec
	   );
}

elsif ($DbSet =~ "prot")
{
    #---------------------+
    # BLAST OPTIONS       |
    #---------------------+
    # TO DO Changed the sequence CDS
    # to an actual protein database
    $BlProg = "tblastx";
#    $BlSuf = "-e 0.00001 -U";
    $BlSuf = "-e 0.00001 -U -m 8";
    $OutExt = "blx";

    @Db = ( 
	    #---------------------+
	    # CDS MODEL SEQS      |
	    # FROM TIGR           |
	    #---------------------+
	    "ATH1_cds",
	    "Os_5_cds"
	    );

}
#-----------------------------+
# REPEAT PROTEIN DATABASES    |
#-----------------------------+
elsif ($DbSet =~ "repx")
{
    #---------------------+
    # BLAST OPTIONS       |
    #---------------------+
    # BLASTn with no lowercase filtering 
    $BlProg = "blastx";
#    $BlSuf = "-e 0.00001";
    $BlSuf = "-e 0.00001 -m 8";
    $OutExt = "blx";
    
    @Db = ( 
            #---------------------+
            # REPEAT DATABASES    |
	    #---------------------+
	    "PTREP_9"              # TREP9 PROTEINS
	    );
}
#-----------------------------+
# UNIPROT PROTEIN DATABASES   |
#-----------------------------+
elsif ($DbSet =~ "embl")
{
    #---------------------+
    # BLAST OPTIONS       |
    #---------------------+
    # BLASTx with lowercase filtering 
    $BlProg = "blastx";
#    $BlSuf = "-e 0.00001 -U";
    $BlSuf = "-e 0.00001 -U -m 8";
    $OutExt = "blx";
    
    @Db = ( 
            #---------------------+
            # UniProt DATABASES   |
	    #---------------------+
	    "UniProt_SProt_34_7",  # 
	    "UniProt_TrEMBL_34_7"  # Translated EMBL, not verified
	    );
}
elsif ($DbSet =~ "repn")
{
    #---------------------+
    # BLAST OPTIONS       |
    #---------------------+
    # BLASTn with no lowercase filtering 
    $BlProg = "blastn";
#    $BlSuf = "-e 0.00001";
    $BlSuf = "-e 0.00001 -m 8";
#    $BlSuf = "-e 0.00001";
    $OutExt = "blo";
    
    @Db = ( 
            #---------------------+
            # REPEAT DATABASES    |
	    #---------------------+
	    "TATrans_200702",      # TriAnnotation Transposons
	    "mips_REdat_4_3",      # MIPS REDat repeats
	    "RB_pln",              # RepBase Plants
	    "SanMiguel_200610",    # SanMiguel
	    "TIGR_Gram_3_3",       # TIGR Gramineae
	    "TREP9_nr",            # TREP 9 Nonredundant
	    "TREP9_total",         # TREP 9 Total
	    "Wessler"              # Wessler dataset
	    );
}else{
    print "A valid database set option has not been used:\n";
    print "Valid options are:\n";
    print "\t-d estn\n";
    print "\t-d estx\n";
    print "\t-d prot\n";
    print "\t-d repn\n";
    print "\t-d repx\n";
} # END OF $DbSet option loading

#-----------------------------+
# DO MATH TO DETERMINE THE    |
# NUMBER OF PROCESSES         |
#-----------------------------+
my $LenQry = @Qry;
if ($LenQry == 0)
{
    print "There are no qry sequences in the directory:\n".
	$QDir."\n";
    exit;
}

my $LenDb =  @Db;
if ($LenDb == 0)
{
    print "There are no databases associated with".
	" the database set $DbSet\n";
    exit;
}

my $NumProc = $LenQry * $LenDb;

for $IndQry (@Qry)
{

    #-----------------------------+
    # GET BASE NAME OF FASTA FILE |
    # AND MAKE OUTPUT DIR         |
    #-----------------------------+
    if ($IndQry =~ /(.*)\.fasta/)
    {
	$QryRoot = $1;
    }elsif($IndQry =~ /(.*)\.fasta\.masked/)
    {
	$QryRoot = $1;
    }else{
	print "INPUT FILE NOT IN EXPECTED FORMAT.\n";
	print "EXPECTED FILE EXTENSION:\n";
	print "\tfasta.masked\n";
	die;
    }
    
    #-----------------------------+
    # CREATE OUTPUT DIR           |
    #-----------------------------+
    $OutDirPath = $QDir.$QryRoot."/";
    unless (-e $OutDirPath)
    {mkdir($OutDirPath, 0755) 
	 || die "Cannot mkdir:\n$OutDirPath: $!";}
    
    for $IndDb (@Db)
    {
	
	$ProcNum++; # Increment the process number

	print "BLAST Process ".$ProcNum." of ".$NumProc."\n";
	
	$QryPath = $QDir.$IndQry;
	$DbPath = $DbDir.$IndDb;
	$TestFile = $DbPath.".nhr";

	# TODO
	# Options here can make use of the database
	# set that is being blasted against (ie. blastx
	# can use the *.blx output while blastn can use *.bln

	$OutPath = $OutDirPath.$QryRoot."_".$IndDb.".".$OutExt;

	#-----------------------------+
	# DOES THE BLAST DB EXIST     |
	#-----------------------------+
	# CHANGE THIS TO AN UNLESS STATEMENT
	# TODO:
	# THIS CHECK ONLY WORKS FOR DNA DATABASES
	# THIS SHOULD BE CHANGED TO A REG EXP
	# OR ADD THE OPTION TO WORK ON PROTEIN NAME
	# EXTENSIONS
#	unless (-e $TestFile ) 
#	{
#	    die "Can not find database:\n$IndDb\n"; 
#	}
	
	#-----------------------------+
	# DOES THE BLAST QRY EXIST    |
	#-----------------------------+
	unless (-e $QryPath)
	{
	    die "Can not find qry file:\n$QryPath\n";
	}

	#------------------------------+
	# PRINT THE BLAST COMMAND      |
	#------------------------------+
	$BlastCmd = "blastall -p ".$BlProg." -i $QryPath -d ".
	    "$DbPath -o $OutPath $BlSuf";
	print wrap("\t", "\t", $BlastCmd );
	print "\n";

	#------------------------------+
	# RUN THE BLAST COMMAND        |
	#------------------------------+
	if (! $test)  # If this is not a test then run the BlastCmd
	{
# Declawed here to make the GFF files 
#	    system ($BlastCmd);
	}
	

	#-----------------------------+
	# CONVERT BLAST OUTPUT TO GFF |
	#-----------------------------+
	# 04/18/2007
	# Check that the BLAST input file exists
	# and is not zero length
        # I think that -s will ignore 0 size files
	if ($MakeGff)
	{
	    if (-s $OutPath)
	    {
		print "CONVERTING BLAST TO GFF\n";
		my $GffOutFile = $OutDirPath.$QryRoot."_AllBlast.gff";
		my $GffCmd = "Blast2Gff.pl -i $OutPath -o $GffOutFile".
		    " -p $BlProg -d $IndDb -s $QryRoot -a";
		print "GFF CMD:$GffCmd\n";
		system ($GffCmd);
	    } # End of if $OutPath file size not zero (-s)
	} # End of if MakeGff
	
    } # End of for each database loop

} # End of for each query seq loop

exit;


#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


=head1 HISTORY

Started: 05/10/2006

Updated: 06/04/2007

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 05/10/2006
# - Program was started to create an easy way to blast all
#   of the query datasets of interest by all of the repeat
#   databases of interest. The BLAST reports should then
#   be parsed and uploaded to the dbASGR database.
#
# 02/13/2007
# - Existing program was modified for the PAWS program
#
# 02/17/2007
# - Added the search for fasta files instead of hard coded
#   list of query files.
#
# 02/18/2007
# - Added the d parameter and changed code to load
#   different parameters for the different database set
# - General code cleanup
#
# 02/19/2007
# - Added tblastx option for the est databases
# - Added blastx for repeat protein databases (PTREP_9)
# - Changed vars for input to the -d flag to allow for both
#   BLASTN, TBLASTX, and BLAST as appropriate
#
# 02/20/2007
# - Added option for the BLASTX against the UNIPROT database
#   The version of uniprot I am using is dated 02/20/2007
#   and listed as version 34.6
#   This option is -d embl
#
# 03/01/2007
# - Added finished genome CDS protein models, (A BLASTX) 
#   this is the variable -d genx. The BLAST databases include
#    - TAIR6, OS_TIGR, OS_RAP1
# - Added blastn against potential contamination
#    - E.coli, UniVec
# 03/07/2007
# - Added the query flag -q
# 
# 04/18/2007
# - Converted all blast output to m8 output. This will make the
#   program compatible with the Blast2Gff PERL script I wrote
# - Adding -g flag, and the ability to convert the output to
#   to gff foramt.
# - Changing the GFF output to append to a single GFF file
#
# 04/24/2007
# - Trying to fix problems with GFF parsing of blastx results
#
# 06/04/2007
# - Adding POD documentation
