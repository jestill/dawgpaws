#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_rephmmer.pl - Run repeat hmmer searches in batch    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/17/2007                                       |
# UPDATED: 09/17/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  Run hmmer against repeat hmmer models in batch mode,     |
#  produce results in gff format if requested.              |
#                                                           | 
# VERSION: $Rev$                                            |
#                                                           |
#-----------------------------------------------------------+
# CONFIG FILE CAN INCLUDE
# 
# Name for the paramter set
# Folder containing the models for this set 
#   ie Mites, Mules, LTRS, etc.
# -A Limit to n best domains
# -E Evalue cutoff
# -T T bit threshold
# -Z # seqs for E-Value calc





package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
use Bio::Tools::HMMER::Results;# Bioperl HMMER results parser
use Text::Wrap;                # Allows word wrapping and hanging indents
                               # for more readable output for long strings.



#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;
my $config_file;

# Array of rephmmer parameters
my @rh_params = ();            # Array of rephmmer parameters

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"  => \$infile,
                    "o|outdir=s" => \$outfile,
		    "c|config=s" => \$config_file,
		    # ADDITIONAL OPTIONS
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
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

# Throw error if required arguments not present
if ( (!$indir) || (!$outdir) || (!$config_file) ) {
    print "\a";
    print STDERR "ERROR: An input directory must be specified" if !$indir;
    print STDERR "ERROR: An output directory must be specified" if !$outdir;
    print STDERR "ERROR: A config file must be specified" if !$config_file;
    print_help("full");
    exit;
}



#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}


#-----------------------------+
# LOAD THE CONFIG FILE        |
#-----------------------------+
$i=0;
my $config_line_num=0;

open (CONFIG, "<$config_file") ||
    die "ERROR Can not open the config file:\n $config_file";

while (<CONFIG>) {
    chomp;
    $config_line_num++;
    unless (m/^\#/) {
       	my @in_line = split;           # Implicit split of $_ by whitespace
	my $num_in_line = @in_line; 
	
	if ($num_in_line == 10) { 
	    $rh_params[$i][0] = $in_line[0] || "NULL";  # Name
	    $rh_params[$i][1] = $in_line[1] || "NULL";  # MIN-MEM
	    $rh_params[$i][2] = $in_line[2] || "NULL";  # MIN-MEM-DIST
	    $rh_params[$i][3] = $in_line[3] || "NULL";  # MAX-MEM-DIST
	    $rh_params[$i][4] = $in_line[4] || "NULL";  # MAX-MEM-GAP
	    $rh_params[$i][5] = $in_line[5] || "NULL";  # MIN-LEN-LTR 
	    $rh_params[$i][6] = $in_line[6] || "NULL";  # MAX-LEN-LTR
	    $rh_params[$i][7] = $in_line[7] || "NULL";  # RANGE-BIN 
	    $rh_params[$i][8] = $in_line[8] || "NULL";  # MIN LEN ORF
	    $rh_params[$i][9] = $in_line[9] || "NULL";  # MAX E Value HMM
	    $i++;
	} # End of if $num_in_line is 10
	else {
	    print "\a";
	    print STDERR "WARNING: Unexpected number of line in config".
		" file line $config_line_num\n$config_file\n";
	}

   } # End of unless comment line
} # End of while CONFIG file
close CONFIG;

# Number of parameter sets specified in the config file
my $num_par_sets = $i;

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $num_files = @fasta_files;

my $num_proc_total = $num_files * $num_par_sets;


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"MyProg.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS::\n".
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



sub RepHmmerRun
#-----------------------------+
# RUN THE HMMER PROGRAM       |
#-----------------------------+ 
{
    #-----------------------------+
    # VARS PASSED TO THE SUBFUN   |
    #-----------------------------+
    my $seq_name = $_[0];           # The name of the individual query ie PPAN058CON001
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

	$OutPath = $HmmOutDir."/".$seq_name."_".$IndMod.".hmmout";

	#-----------------------------+
	# DOES THE HMM MODEL EXIST    |
	#-----------------------------+
	if (-e $ModFile ) {
	    #print "DB: $IndDb exists\n";
	}
	else {
	    die "Can not find model:\n$IndMod\n"; 
	}
	
	#-----------------------------+
	# DOES THE HMM QRY SEQ EXIST  |
	#-----------------------------+
	if (-e $QryPath) {
	    #print "QRY: $QryPath exists\n";
	}
	else {
	    die "Can not find qry file:\n$QryPath\n";
	}
	
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
    print "\tPARSING THE HMM_".$seq_name."_".$class." OUTPUT\n";
    &RepHmmerParse ( $seq_name, $WorkDir, $HmmOutDir, $class );


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
    my $seq_name = $_[0];
    my $WorkDir = $_[1];
    my $HmmDir = $_[2];
    my $class = $_[3];


    # File to write the parsed output to
    my $HmmOutFile = $WorkDir.$seq_name."_".$class.".txt";

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
    print "\tPAN:       \t".$seq_name."\n";
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


=head1 NAME

batch_rephmmer.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InFile -o OutFile

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 COMMAND LINE ARGUMENTS

=head 2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 Additional Options

=over

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

=head1 HISTORY

STARTED: 09/17/2007

UPDATED: 09/17/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/17/2007
# - Program started
# - 
