#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# BLAST TO GFF FILE
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/06/2006                                       |
# UPDATED: 07/06/2006                                       |
#                                                           |  
# DESCRIPTION:                                              | 
# Parse blast from single large contig against a blast db   |
# of multiple reads. Output is a single gff format file     |
# suitable for viewing in Apollo.                           |
#                                                           |
# USAGE:                                                    |
#                                                           |
#-----------------------------------------------------------+

print "The program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Std;               # Allows to get options from the command line
use Bio::SearchIO;             # Parse BLAST output
use DBI();                     # Connect to the database

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my $BlastFile = "/home/jestill/projects/asgr/Asgr_c_Net_BigCat/CPAN001/".
    "contig02/contig02_CPAN001.blo";
my $GFFOut = "/home/jestill/projects/asgr/Asgr_c_Net_BigCat/CPAN001/".
    "contig02/contig02_CPAN001.gff";
my $MaxE = 0.00001;
my $MinQryLen = 100;

#-----------------------------+
# DB RELATED VARS             |
#-----------------------------+ 
my $DbName = "dbASGR";             # Database name to connect to
my $tblSrcTable = "tblSeqData";    # Table with src sequence data
my $tblDatCat = "tblDatCat";       # Table with information about data
                                   # categories. Can include color
                                   # and symbol for drawing.
my $DbUserName = "jestill";
my $Species = "cen";               # Species to look up id number for
#my $Species = "pen";

#-----------------------------+
# COMMAND LINE VARIABLES      |
#-----------------------------+
my %Options;                  # Options hash to hold the options from the command line
getopts('q', \%Options);      # Get the options from the command line
my $quiet = $Options{q};

#-----------------------------+
# SET VARIABLE SCOPE          |
#-----------------------------+
my $HitName;                   # Name of the hit

#-----------------------------+
# GET USER PASSWORD           |
#-----------------------------+
print "\nPassword for $DbUserName\n";
system('stty', '-echo') == 0 or die "can't turn off echo: $?";
$DbUserPassword = <STDIN>;
system('stty', 'echo') == 0 or die "can't turn on echo: $?";
chomp $DbUserPassword;

#-----------------------------+
# OPEN OUTPUT FILE            |
#-----------------------------+
open (GFFOUT, ">$GFFOut") ||
    die "Can not open file:\n $GFFOut\n";

my $dbh = DBI->connect("DBI:mysql:database=$DbName;host=localhost",
		       $DbUserName, $DbUserPassword,
		       {'RaiseError' => 1});


# OUTPUT IN FORM
#LTR/gypsy	RepeatMasker: AllRep	LTR/gypsy	27676	27772	268	.	.

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

my $BlastReport = new Bio::SearchIO ( '-format' => 'blast',
				      '-file'   => $BlastFile,
				      '-signif' => $MaxE,
				      '-min_query_len' => $MinQryLen) ||
    die "Could not open BLAST input file:\n$BlastFile.\n";

while ($BlastResult = $BlastReport->next_result())
{
    	while ( $BlastHit = $BlastResult->next_hit())
	{
	    while ($BlastHSP = $BlastHit->next_hsp())
	    {

		my $HitName = $BlastHit->name();
		my $Bac = &GetBacName( $HitName, $Species);

		#-----------------------------+
		# PRINT OUTPUT TO TERMINAL    |
		#-----------------------------+ 
		print $BlastHit->name()."\n";
		print "\t".$BlastHSP->start('query')."\n";
		print "\t".$BlastHSP->end('query')."\n";
		print "\t".$Bac."\n";

		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		#-----------------------------+
		#print GFFOUT $BlastHit->name()."\t"."BLASTN\tBLASTNHIT\t".
		#    $BlastHSP->start('query')."\t".
		#    $BlastHSP->end('query')."\t".
		#    $BlastHSP->evalue()."\t".
		#    ".\t.\n";

		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		# WITH BAC DATA               |
		#-----------------------------+
		#print GFFOUT $BlastHit->name()."\t"."BLASTN\tBLASTNHIT\t".
		#    $BlastHSP->start('query')."\t".
		#    $BlastHSP->end('query')."\t".
		#    $BlastHSP->evalue()."\t".
		#    ".\t.\t".
		#    "$Bac\n";

		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		# WITH BAC DATA               |
		#-----------------------------+
		# Changing BLASTN to the Bac Name appears to allow 
		# these to be drawn on different levels.
		print GFFOUT $BlastHit->name()."\t".
		    "$Bac\t".
		    "$Bac\t".
		    $BlastHSP->start('query')."\t".
		    $BlastHSP->end('query')."\t".
		    #$BlastHSP->evalue()."\t".
		    # Apollo seems to work better with the 
		    # bit score then with the evalue
		    $BlastHSP->score()."\t".
		    ".\t.\t".
		    "$Bac\n";

	    } # End of while next hsp

	} # End of while next hit
} # End of while next result


exit;



#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub GetBacName
{
#-----------------------------+
# GET THE BAC ID FOR THE      |
# ID NUMBER AND SPECIES NAME  |
# PASSED TO THE SUBFUNCTION   |
#-----------------------------+

    my $SpecNum = $_[0];
    my $SpecName = $_[1];
    my $result;

    my $SelectSQL = "SELECT bac FROM ".$tblSrcTable.
	" WHERE spec_num = '".$SpecNum."'".
	" AND species = '".$SpecName."'";

    if (! $quiet)
    {
	print "\tSQL STATEMENT\n$SelectSQL\n";
    }

    $dbh->do($SelectSQL);
    $cur = $dbh->prepare($SelectSQL);
    $cur->execute();
    @row=$cur->fetchrow;
    $result=$row[0] || "SQL ERROR";
    $cur->finish();
    return $result;

}


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/06/2006
# - Program written to parse BLAST file and send out to
#   Apollo compatible GFF formatted file. 
#
# 07/07/2006
# - Added GetBacName subfunction
# - Changed score used to the BIT Score. This works better
#   with apollo



