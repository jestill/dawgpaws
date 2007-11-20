#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# HMMER FOR FASTA                                           |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 09/21/20060                                      |
# UPDATED: 09/21/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a fasta file, create a individiual fasta file      |
#  for each of the records and then run and parse HMMER     |
#  to discover Transposable Elements for which there        |
#  is a HMM profile.                                        |
#                                                           |
#                                                           |
# USAGE:                                                    |
#  HMM4FASTA.pl                                             |
#                                                           |
#-----------------------------------------------------------+
#$test = "true";
#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;               # Required to work with FASTA files 
use Getopt::Std;               # Allows to get options from the command line
use Bio::Tools::HMMER::Results;
use Text::Wrap;                # Allows word wrapping and hanging indents
                  # for more readable output for long strings.

# OPEN A METAOUTPUT FILE FOR TESTING
open (METAOUT, ">/home/jestill/projects/asgr/MetaNewTest.html") ||
    die "Can not open metatout file\n";

open (TABOUT, ">/home/jestill/projects/asgr/200609/Mites_Mules_Again.txt") ||
    die "Can not open tabout file\n";
print TABOUT "NAME\tCLASS\tSTART\tEND\tLENGTH\t".
    "BIT-SCORE\tE-VALUE\tREP-NAME\n";

# FOR TEST PPAN058CON001
# This will get vars like ($CurCon, $OutPath, $OutDir, $Class);
# Where $CurCon is the name of the PAN CONTIG
#       $OutPaht is the fasta file to analyze
#       $OutDir is the directory holding the output
#         A subdir to this will need to be created


#-----------------------------------------------------------+
# OPEN THE FASTA FILE                                       |
#-----------------------------------------------------------+ 
$BaseDir = "/home/jestill/projects/asgr/200609/ind/";
$ClustFile = "/home/jestill/projects/asgr/200609/asgr_con/asgr_con.fasta";
$informat = "fasta";
my $RecordNum = 0;
my $OutFormat = "fasta";
my $seq_in = Bio::SeqIO->new('-file' => "<$ClustFile",
                             '-format' => $informat);

while (my $inseq = $seq_in->next_seq)
{
    $RecordNum++;
#    if ($RecordNum == 100){
#	print "\a\a\a";
#	print "ALARM ONE\n";
#	print "\a\a\a";
#	sleep 1;
#	print "ALARM TWO\n";
#	print "\a\a\a";
#	sleep 1;
#	print "ALARM THREE\n";
#	print "\a\a\a";
#	exit;
#    }
    my $SeqName = $inseq->primary_id();
    print "PROCESSING ".$RecordNum."\t".$SeqName."\n";

     
    #-----------------------------+
    # CREATE AN OUTPUT DIR OT HOLD|
    # THE RESULTS IF NEEDED       |
    #-----------------------------+
    my $SeqDir = $BaseDir.$SeqName;
    print "\tDIR:".$SeqDir."\n";
    mkdir $SeqDir, 0777 unless (-e $SeqDir); # set permissions    
    
    my $SeqPath = $BaseDir.$SeqName."/".$SeqName.".fasta";
    
    my $seq_out = Bio::SeqIO->new('-file' => ">$SeqPath",
				  '-format' => $OutFormat)
	|| die "Can not open:\n$SeqPath\n";
    
    #-----------------------------+
    # WRITE THE RECORD THE OUT    |
    # FILE PATH                   |
    #-----------------------------+
    $seq_out->write_seq($inseq) ||
	die "Can not write out sequence:\n".$SeqPath."\n"; 


    #-----------------------------+
    # RUN HMMER ON THE MITE MODELS|
    #-----------------------------+
    print "\tANALYZING THE MITE MODELS\n";
    &RepHmmerRun ($SeqName,        # Name of the query sequence
		  $SeqPath,        # Path of the sequence file
		  $SeqDir,         # Path of the dir to store HMMER output 
		  "MITE");         # CLASS of TE for the            

    #-----------------------------+
    # RUN HMMER ON THE MULE MODELS|
    #-----------------------------+    
    print "\tANALYZING THE MULE MODELS\n";
    &RepHmmerRun ($SeqName,        # Name of the query sequence
		  $SeqPath,        # Path of the sequence file
		  $SeqDir,         # Path of the dir to store HMMER output 
		  "MULE");         # CLASS of TE for the            

}

# SOUND THE ALARM FOR FINISHED
print "\a\a\a";
sleep 1;
print "\a\a\a";
sleep 1;
print "\a\a\a";

exit;


#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+



sub RepHmmerRun
#-----------------------------+
# RUN THE HMMER PROGRAM       |
#-----------------------------+ 
{
    #-----------------------------+
    # VARS PASSED TO THE SUBFUN   |
    #-----------------------------+
    my $PanName = $_[0];           # The name of the individual query ie PPAN058CON001
    my $QryPath = $_[1];           # The path to to the qry fasta file
    my $WorkDir = $_[2];           # The work dir,
                                   # Will need to make a hmm_class dir here
    my $class = $_[3];             # The class of repeast to search
                                   # vars can include MITE,MULE

    #-----------------------------+
    # VAR SCOPE AND INITIALIZE    |
    #-----------------------------+
    $ProcNum = 0;                  # Process number starts at zero
    my ( $HmmCmd, $DbPath, $OutPath );  # Declare scope for varaiables used later
		  
    #-----------------------------+
    # HMMER MODELS                |
    #-----------------------------+    
    # FIRST DETERMINE THE APPROPRIATE DIR GIVEN
    # THE MODEL SET THAT IS BEING SEARCHED (MITE/MULE)
    #-----------------------------+
    # HMMER MODELS                |
    #-----------------------------+    
    # FIRST DETERMINE THE APPROPRIATE DIR GIVEN
    # THE MODEL SET THAT IS BEING SEARCHED (MITE/MULE)
    if ($class =~ "MITE")
    {
	$ModDir = "/home/jestill/HMMData/db/hmm/mite_models/";
    }elsif ($class =~ "MULE"){
	$ModDir = "/home/jestill/HMMData/db/hmm/mule_models/";
    }elsif ($class =~ "TPASE"){
	$ModDir = "/home/jestill/HMMData/db/hmm/tpase_models/";
    }elsif ($class =~ "PFAM"){
	$ModDir = "/home/jestill/HMMData/db/hmm/pfam/";
    }else{
	$ModDir = "/home/jestill/HMMData/db/hmm/mite_models/"; # Default is to just use MITES
    }
    
    # Open the appropriate dir and load files 
    opendir( MODDIR, $ModDir );
    my @Mod = grep !/^\.\.?$/, readdir MODDIR ;
    closedir( MODDIR );    

    #-----------------------------+
    # CREATE A HMMER OUTPUT DIR   |
    #-----------------------------+
    my $HmmOutDir = $WorkDir.$class;
    mkdir $HmmOutDir, 0777 unless (-e $HmmOutDir); # set permissions
    
    # Determine the total of HMMER queries that will be run
    # Since a single sequence is passed to the subfun this
    # will just be the same as the number of models to test
    # against.
    my $LenMod =  @Mod;
    my $NumProc = $LenMod;
    
    for $IndMod (@Mod)
    {
	
	$ProcNum++; # Increment the process number

	
	# HMMER MODEL PATH
	$ModPath = $ModDir.$IndMod;
	$ModFile = $ModPath;

	$OutPath = $HmmOutDir."/".$PanName."_".$IndMod.".hmmout";

	#-----------------------------+
	# DOES THE HMM MODEL EXIST    |
	#-----------------------------+
	if (-e $ModFile ) 
	{
	    #print "DB: $IndDb exists\n";
	}
	else 
	{die "Can not find model:\n$IndMod\n"; }
	
	#-----------------------------+
	# DOES THE HMM QRY SEQ EXIST  |
	#-----------------------------+
	if (-e $QryPath)
	{
	    #print "QRY: $QryPath exists\n";
	}
	else
	{die "Can not find qry file:\n$QryPath\n";}
	
	#------------------------------+
	# PRINT THE HMMER COMMAND      |
	#------------------------------+
	$HmmCmd = "hmmsearch " . 
	    "--domT 2 $ModFile $QryPath >$OutPath";
	
	#-----------------------------+
	# PRINT STATUS OF EACH HMM    |
	# PROCESS                     |
	#-----------------------------+	
	#print "HMM Process ".$ProcNum." of ".$NumProc."\n";
	#print wrap("\t", "\t", $HmmCmd );
	#print "\n";
	
	#------------------------------+
	# RUN THE HMMER COMMAND        |
	#------------------------------+
	if (! $test)  # If this is not a test then run the BlastCmd
	{
	    system ($HmmCmd);
	}
	
	
    } # End of for each database loop


    # ONCE THE ENTIRE SET HAS RUN IT IS TIME TO PARSE
    print "\tPARSING THE HMM_".$PanName."_".$class." OUTPUT\n";
    &RepHmmerParse ( $PanName, $WorkDir, $HmmOutDir, $class );


} # END OF RepHmmerRun

sub RepHmmerParse
{
#-----------------------------+
# PARSE THE OUTPUT FROM A     |
# HMMER RUN AGAINST A SET OF  |
# HMM PROFILES                |
#-----------------------------+ 

    #-----------------------------+
    # VARIABLES                   |
    #-----------------------------+
    # THE FOLLOWING SHOULD BE PASSED TO THE SUBFUNCTION
    my $PanName = $_[0];
    my $WorkDir = $_[1];
    my $HmmDir = $_[2];
    my $class = $_[3];


    # File to write the parsed output to
    my $HmmOutFile = $WorkDir.$PanName."_".$class.".txt";

    if ($class =~ "MITE")
    {
	$dataset = "hmm_mite";
    }
    elsif ($class =~ "MULE")
    {
	$dataset = "hmm_mule";
    }else{
	$dataset = "UNK";
    }   
    
    my $BestName = "";             # Name of the best hit
    my $FilePath;                  # The file path for the inidividual HMM output
    my $TotRes = '0';              # Total number of results for the sequence
    my $BestBit = '0';             # Best BIT Score
    my $BestHit;                   # The name of the Best Hit
    my $BestEval = 'null'; # Set scope of this variable

    #-----------------------------+
    # LOAD FILES TO PARSE INTO    |
    # THE @HmmFiles ARRAAY        |
    #-----------------------------+
    opendir( HMMDIR, $HmmDir );
    my @HmmFiles = grep !/^\.\.?$/, readdir HMMDIR ;
    closedir( HMMDIR );
    
    # Open up an Output file
    open ( OUT, ">".$HmmOutFile);
    print OUT "SEQNAME      \tSTART\tEND\tSCORE\tEVAL\tHITNAME\n";
    print OUT "=============\t=====\t===\t=====\t====\t=======\n";
    
    #-----------------------------------------------------------+
    # FOR EACH FILE IN THE ARRAY OF FILES                       |
    #-----------------------------------------------------------+
    for $IndFile (@HmmFiles)
    {
	
	$FilePath = $HmmDir."/".$IndFile;
	
        # OPEN THE HMM RESULT AS HMMER RESULTS OBJECT
	my $HmmRes = new Bio::Tools::HMMER::Results ( -file => $FilePath ,
						      -type => 'hmmsearch') 
	    || die "Could not open file\n$FilePath\n";
	
	my $NumRes = $HmmRes->number;
	$TotRes = $TotRes + $NumRes;
	
	# ONLY PRINT OUTPUT FOR QUERIES WITH MATCHES
	if ($NumRes >> 0) #08/07/2006
	{	              #08/07/2006
	    foreach $seq ( $HmmRes->each_Set ) 
	    {
		foreach $domain ( $seq->each_Domain ) 
		{		
		    #my $CurBit = $domain->bits;
		    my $CurName = $domain->hmmname;
		    $CurName =~ m/.*\/(.*)\.hmm/;  # Returns what I want
		    $CurName = $1;
		    # RECORD THE NAME AND SCORE OF THE
		    # BEST HIT
		    if ($domain->bits > $BestBit)
		    {
			# ASSUMES BIT SCORE AND E VALUE
			# HAVE THE SAME RANK
			$BestBit = $domain->bits;
			$BestName = $CurName;
			$BestEval = $domain->evalue || "NULL";
		    }

		    # Determine the length of the match of the
		    # element length
		    my $ElmLen = $domain->end - $domain->start;
		    
		    print OUT $seq->name."\t";
		    print OUT $domain->start."\t";
		    print OUT $domain->end."\t";
		    print OUT $domain->bits."\t";
		    print OUT $domain->evalue."\t";
		    print OUT $CurName."\n";

		    #-----------------------------+
		    # PRINT RESULTS TO FILE THAT  |
		    # CONTAINS ALL POSITIVE HITS  |
		    #-----------------------------+
		    print TABOUT $seq->name."\t";
		    print TABOUT $class."\t";       # Show if this is a MITE or a MULE
		    print TABOUT $domain->start."\t";
		    print TABOUT $domain->end."\t";
		    print TABOUT $ElmLen."\t";
		    print TABOUT $domain->bits."\t";
		    print TABOUT $domain->evalue."\t";
		    print TABOUT $CurName."\n";
		    

		    		    
		} # End of for each domain in the HMM output file
	    } # End of if greater then zero hits 08/07/2006
	    
	} # End of for each seq in the HMM output file
	
    }
    
    print "\n\t===================================\n";
    print "\t$class HMMER RESULT\n";
    print "\t===================================\n";
    print "\tPAN:       \t".$PanName."\n";
    print "\tDATASET:   \t".$dataset."\n";
    print "\tCLASS:     \t".$class."\n";
    print "\tTOTAL:     \t". $TotRes."\n";
    print "\tBEST HIT:  \t".$BestName."\n";
    print "\tBEST BIT:  \t".$BestBit."\n";
    print "\tBEST EVAL: \t".$BestEval."\n";
    print "\t===================================\n";

    
    # PRINT APPROPRIATE OUTPUT TO THE METAFILE IF
    # HITS WERE FOUND IN THE DATABASE
    if ($TotRes > 0 )
    {
	# SOUND THE ALARM TO INDICATE A CONTIG WITH A POSITIVE HIT
	print "\a";

	print METAOUT "<TR>".
	    "<TD align=right>".$dataset."</TD>".      # Repeat Database
	    "<TD>".$class."</TD>".                    # Repeat class
	    "<TD>".$BestName."</TD>".                 # Repeat name
	    "<TD>".$TotRes."</TD>".      # Number of total hits
	    "<TD>".$BestBit."</TD>".
	    "<TD>".$BestEval."</TD>".
	    "</TR>\n";

    }
    
    
    close OUT;

    #-----------------------------------------------------------+
    # TO DO
    #-----------------------------------------------------------+
    # - Figure out if there are multiple MITES on a single
    #   contig. This is very likely given the short length
    #   of MITES. Currently just have to manually look at the
    #   output.
 

} # END OF THE REP HMMER PARSE SUBFUNCTION



#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
# 09/21/2006
#  - Program started; modified from other HMMER parsing
#    programs I wrote for the RepMiner pipeline.
#
