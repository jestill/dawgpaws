#!/usr/bin/perl -w

###############################################################################
#                                                                             #
# MySeqMan :: GenBankTaxUpload                                                #
# Author: Jamie Estill                                                        #
# Started: 3/1/2003                                                           #
# Last Updated: 5/18/2005                                                     #
# Automates construction of an NCBI database. The data are downloaded from    #
# NCBI and unzipped usiing the gzip utility. The data are then transformed to #
# a pipe delimited text file suitable for upload to a MySQL database. Data    #
# are then uploaded to the appropriate tables in a MySQL database.            #
#                                                                             #
###############################################################################

# CURRENT NEEDS 1/16/2005
# Add an index to the tax_id field of tblNodes and tblTax to allow for 
# quicker look up of information from these tables.
# This will speed up queries searching for this information
# 2/16/2006
# Need to create the taxonomy directory on demand if it does not already exist

# 5/18/2005
# - Fixed the subfunction that checks for the existence of tables
# 5/22/2005
# - Cleaned up code
# 6/8/2005
# - Added quiet mode
# - Added command line options
# - Changeing to a single dbh connection
#
#  NEEDS
# Need to change this to a single dbh connection

# THIS ONE INCLUDES THE FIX/KLUDGE FOR THE BUG IN THE NODES FIELD
# This problem was due to there being too many fields for the database structure

# LOAD REQUIRED LIBRARIES
use Net::FTP;
use DBI;
use Bio::SeqIO;
use Getopt::Std;

my %Options;                                               # Declare the Options hash to hold the options from the comma\
getopts('u:p:q', \%Options);

# SET VARIABLES
$TaxDir = "/home/jestill/NCBI/taxonomy";
$TimeStartAll = time;

my $quiet;
$DbUserName = $Options{u}|| die "ERROR: A user name for the database must be supplied";
$DbUserPassword = $Options{p}|| die "ERROR: A password must be supplied";
$quiet = $Options{q};

# Check to make sure the directory for storing the downloaded files exists
# The following command is not working, 
#mkdir($TaxDir,'755') || print "$TaxDir Already Is Present";


# MAKE THE TABLE tblNodes with a subset of the available fields

if (! $quiet) {print "Making the table tblNodes. \n"};

my $dbh = DBI->connect("DBI:mysql:database=dbNCBI;host=localhost",
                           $DbUserName, $DbUserPassword,
		       {'RaiseError' => 1});

# RUN A DBI TRACE
# This is used for error checking
# Only need to run this if there are problems with the DBI interface
# DBI->trace(5);

# GET DATE AND TIME USING THE GET TIME SUBFUNCTION
if (! $quiet) {print "Getting the time in a usable format.\n";}
($RealMonth, $Day, $Year, $Hour, $Minute, $Second, $ampm) = GetTime(0,1,2,3,4,5,6);

if (! $quiet) {print "Opening the log file.\n";}

open(NCBILOG, ">NCBILog");                       # NCBI Log file use to keep track of what is going on

print NCBILOG "\n", "======================================================\n";
print NCBILOG "START OF GenBankTaxonomyUpload ";
printf NCBILOG ("%02d/%02d/%04d", $RealMonth, $Day, $Year);
print NCBILOG " ";
printf NCBILOG ("%02d:%02d:%02d $ampm", $Hour, $Minute, $Second);
print NCBILOG "\n";
print NCBILOG "======================================================\n\n";

# DOWNLOAD THE TAXONOMIC DATA FROM NCBI
if (! $quiet) {print "Downloading Taxonomic Data.\n";}
chdir $TaxDir;
&GetTaxData($TaxDir, "gi_taxid_nucl.dmp.gz");
&GetTaxData($TaxDir, "taxdump.tar.Z");

if (! $quiet) {print "Untarring the Taxonomic Data.\n";}
# May need to add a quiet form of tar
# IF NOT QUIET USE VERBOSE TAR
if (! $quiet) 
{
  system ("tar -xvf ".$TaxDir."/taxdump.tar");
}

# OTHERWISE USE QUIET TAR
if ($quiet)
{
  system ("tar -xf ".$TaxDir."/taxdump.tar");
}


system ("rm taxdump.tar");

# MAKE THE tblNames Table
if (! $quiet) {print "Making the table tblNames. \n";}
chdir $TaxDir;
$MakeThis = "CREATE TABLE tblNames(tax_id CHAR(6) NOT NULL, name_txt TEXT, unique_name TEXT, name_class TEXT )";
$DoThis = "LOAD DATA LOCAL
  INFILE '\/home\/jestill\/NCBI\/taxonomy\/names.dmp'
  INTO TABLE tblNames
  FIELDS TERMINATED BY '\t|\t'
  LINES TERMINATED BY '\t|\n'";
&MakeTaxTable($TaxDir, "tblNames", $MakeThis, $DoThis, $dbh);



# MAKE THE TABLE tblNodes with a subset of the available fields
if (! $quiet) {print "Making the table tblNodes. \n";}

chdir $TaxDir;
$MakeThis = "CREATE TABLE tblNodes
(
  tax_id INT(10) NOT NULL,
  parent_tax_id INT(10),
  rank CHAR(20),
  division_id ENUM('0','1','10','2','3','4','5','6','7','8','9') NOT NULL,
  inherited_div_flag CHAR(1),
  left_id INT(10),
  right_id INT(10)
)";
&MakeNodesTable($TaxDir, "tblNodes", $MakeThis, $dbh);

#$a=1;

open(NODES, "<$TaxDir/nodes.dmp") or die "Couldn't open data file $TaxDir/nodes.dmp: $\n";
while (<NODES>) {
    # The following removed on 4/23/03
    #$a = 1 + $a;
    #print "Record: ", $a, "\n";
    
# oringinally was a push command as
    my @data = split(/\s*\|\s*/o, $_);
    #
    #
    # THE FOLLOWING ARE FOR ERROR CHECKING OF THE SQL COMMAND
    #$TheCommand = "INSERT INTO tblNodes VALUES ("
    #                  .$data[0].
    #              ", ".$data[1].
    #              ", '".$data[2].
    #              "', ".$data[4].
    #              ", ".$data[5].
    #              ", ".$data[6].
    #              ", ".$data[8].
    #              ")";                  # End of VALUES   
    #print $TheCommand, "\n";
    #$dbh->do($TheCommand);
    
    $dbh->do("INSERT INTO tblNodes VALUES ( "
                      .$data[0].
                  ", ".$data[1].
                  ", '".$data[2].
                  "', ".$data[4].
                  ", ".$data[5].
                  ", '', ''".
                  " )"                  # End of VALUES
                  );                    # End of dO
}
close(NODES);

# ADD INDEXES TO THE TABLE tblNodes - added 4/14/03
if (! $quiet) {print "Indexing division_id on the table tblNodes \n";}
$dbh->do("ALTER TABLE tblNodes ADD INDEX ( division_id )");
if (! $quiet) {print "Indexing tax_id on the table tblNodes \n";}
$dbh->do("ALTER TABLE tblNodes ADD INDEX ( tax_id )");
if (! $quiet) {print "Indexing parent_tax_id on the table tblNodes \n";}
$dbh->do("ALTER TABLE tblNodes Add INDEX ( parent_tax_id)");   # This added 6/27/2003

#SET THE ROOT NODE PARENT ID TO NULL, THIS WILL ALLOW THE TREE WALKING ALGORITHM TO WORK 6/28/2003
$dbh->do("UPDATE `tblNodes` SET `parent_tax_id` = NULL WHERE `tax_id` = '1'");

#
# MAKE THE TABLE tblDivision
#
chdir $TaxDir;
if (! $quiet) {print "Making the table tblDivision \n";}
$MakeThis = "CREATE TABLE tblDivision
(
  division_id ENUM('0','1','10','2','3','4','5','6','7','8','9') NOT NULL,
  division_cde ENUM('BCT','INV','MAM','PHG','PLN','PRI','ROD','SYN','UNA','VRL','VRT') NOT NULL,
  division_name TEXT,
  comments TEXT
)";
$DoThis = "LOAD DATA LOCAL 
     INFILE '\/home\/jestill\/NCBI\/taxonomy\/division.dmp'
     INTO TABLE tblDivision 
     FIELDS TERMINATED BY '\t|\t'
     LINES TERMINATED BY '\t|\n'";
&MakeTaxTable($TaxDir, "tblDivision", $MakeThis, $DoThis, $dbh);

# MAKE THE TABLE tblGenCode
chdir $TaxDir;
if (! $quiet) {print "Making the table tblGenCode \n";}
$MakeThis = "CREATE TABLE tblGenCode
( genetic_code_id CHAR(10),
  abbreviation TEXT,
  name TEXT,
  cde TEXT,
  starts TEXT
)";
$DoThis = "LOAD DATA LOCAL
     INFILE '\/home\/jestill\/NCBI\/taxonomy\/gencode.dmp'
     INTO TABLE tblGenCode
     FIELDS TERMINATED BY '\t|\t'
     LINES TERMINATED BY '\t|\n'";
&MakeTaxTable($TaxDir, "tblGenCode", $MakeThis, $DoThis, $dbh);

# MAKE THE TABLE tblDelNodes
chdir $TaxDir;
if (! $quiet) {print "Making the table tblDelNodes \n";}
$MakeThis = "CREATE TABLE tblDelNodes
(
     tax_id CHAR(6)
)";
$DoThis = "LOAD DATA LOCAL
     INFILE '\/home\/jestill\/NCBI\/taxonomy/delnodes.dmp'
     INTO TABLE tblDelNodes 
     FIELDS TERMINATED BY '\t|\t'
     LINES TERMINATED BY '\t|\n'";
&MakeTaxTable($TaxDir, "tblDelNodes", $MakeThis, $DoThis, $dbh);

# MAKE THE TABLE tblMerged
chdir $TaxDir;
if (! $quiet) {print "Making the table tblMerged \n";}
$MakeThis = "CREATE TABLE tblMerged
(
  old_tax_id CHAR(6),
  new_tax_id CHAR(6)
)";
$DoThis = "LOAD DATA LOCAL
     INFILE '\/home\/jestill\/NCBI\/taxonomy\/merged.dmp'
     INTO TABLE tblMerged 
     FIELDS TERMINATED BY '\t|\t'
     LINES TERMINATED BY '\t|\n'";
&MakeTaxTable($TaxDir, "tblMerged", $MakeThis, $DoThis, $dbh);

# MAKE THE TABLE tbl_Gi_TaxID_Nucl
chdir $TaxDir;
if (! $quiet) {print "Making the table tbl_Gi_TaxID_Nucl. \n";}
$MakeThis = "CREATE TABLE tbl_Gi_TaxID_Nucl
(
  gi INT(8) NOT NULL,
  tax_id INT(6) NOT NULL,
  PRIMARY KEY (gi)
)";
$MakeThisToo = "ALTER TABLE tbl_Gi_TaxID_Nucl
      Max_Rows=1000000000
      AVG_ROW_Length=7";
$DoThis = "LOAD DATA LOCAL
     INFILE '\/home\/jestill\/NCBI\/taxonomy\/gi_taxid_nucl.dmp'
     INTO TABLE tbl_Gi_TaxID_Nucl
     FIELDS TERMINATED BY '\t'
     LINES TERMINATED BY '\n'";
&MakeGiTable($TaxDir, "tbl_Gi_TaxID_Nucl", $MakeThis, $MakeThisToo, $DoThis, $dbh);
# ADD INDEX TO THE TABLE tbl_Gi_TaxID_Nucl - added 4/14/03
# this was removed on 9/17/2003 .. for some reason this stopped working
# this was changed on 9/23/3003 to add an index for the tax_id

if (! $quiet) {print "Adding an index for tax_id in the table tbl_gi_TaxID_Nucl. \n";}
$dbh->do("ALTER TABLE tbl_Gi_TaxID_Nucl ADD INDEX ( tax_id )");

# SEE HOW LONG THIS TOOK

# 7/13/2004
# Remove the flat text files that were dowloaded from genbank
chdir $TaxDir;  # Make sure that we are working in the directory that contains the taxonomy text files that were downloaded from Genbank


# TEMP REMOVE THIS I WANT TO KEEP THE TEXT FILES
#print "Deleting the text files ... \n";
#system ("rm *.dmp" );
#system ("rm *.prt");
#system ("rm readme.txt");
#print ".. completed", "\n";


$TimeEndAll = time;
$TotalSeconds = $TimeEndAll - $TimeStartAll;
$TotalMinutes = $TotalSeconds/60;

print NCBILOG "\nThe total time to complete GenBankTaxonomyUpload was ".$TotalMinutes." minutes\n";

($RealMonth, $Day, $Year, $Hour, $Minute, $Second, $ampm) = GetTime(0,1,2,3,4,5,6);

print NCBILOG "\nEND OF GenBankTaxonomyUpload ";
printf NCBILOG ("%02d/%02d/%04d", $RealMonth, $Day, $Year);
print NCBILOG " ";
printf NCBILOG ("%02d:%02d:%02d $ampm", $Hour, $Minute, $Second);
print NCBILOG "\n";

exit;



##################################################
# SUBFUNCTION TO FTP TAXONOMIC DATA FROM         #
# GENBANK FTP SERVER                             #
##################################################

sub GetTaxData
{
  # GET VARIABLES
  $WorkingDir = $_[0];
  $RequestedFile = $_[1];
  $TimeStartDownload = time;

  # LOG ON TO NCBI FTP SERVER
  print NCBILOG "Logging on to NCBI file server .....";
  my $ftp = Net::FTP->new('ftp.ncbi.nlm.nih.gov');
  $ftp->login('anonymous', 'anonymous');
  $ftp->cwd('/pub/taxonomy');
  $ftp->binary;
  print NCBILOG ".. completed connection and directory change","\n";

  print NCBILOG "Downloading ", $RequestedFile;
  $ftp->get($RequestedFile, $WorkingDir."/".$RequestedFile);
  print NCBILOG "completed","\n";

  # LOG OFF FROM THE NCBI FTP SERVER
  $ftp->quit();
  $TimeFinishDownload = time;
  $TotalTimeToDownload = $TimeFinishDownload - $TimeStartDownload;
  print NCBILOG "The time to download ".$RequestedFile." was ".$TotalTimeToDownload." seconds \n";

  # UNPACK FILE
  print NCBILOG "Unzipping ", $RequestedFile, ".....";
  system ("gunzip -f ".$WorkingDir."/".$RequestedFile);
  print NCBILOG "completed","\n";

}

##################################################
# SUBFUNCTION TO MAKE MYSQL TABLES IF NEEDED     #
##################################################

sub MakeTaxTable
{
  $WorkingDir = $_[0];                             # Source directory of the file to be uploaded
  $TableName = $_[1];                              # The name of the table to be created
  $TableMakeCommand = $_[2];                       # Command used to create the table
  $TableLoadCommand = $_[3];                       # Command used to load data into the table
  my $dbh = $_[4];

  #DELETE THE TABLE IF IT ALREADY EXISTS
  if (&does_table_exist( $dbh, $TableName ))
  {
    print NCBILOG "The table ".$TableName." already exists and must be deleted ..... ";
    $dbh->do("DROP TABLE ".$TableName);
    print NCBILOG "completed.", "\n";
  }

  print NCBILOG "Making the table ".$TableName." ....";
  $dbh->do($TableMakeCommand);
  $dbh->do($TableLoadCommand);
  print NCBILOG ".. completed \n";

  # Indexes added 5/22/2005 to speed up queries that look up names
  # Add index to tax_id
  # $dbh->do("ALTER TABLE ".$TableName." ADD INDEX ( tax_id )");
  
  # Add index to name_class, a number must be specified to indicate how much of
  # the text field will be index, choose 10 at random this may be able to be smaller
  # $dbh->do("ALTER TABLE ".$TableName." ADD INDEX ( name_class(10) )");

}

##################################################
# SUBFUNCTION TO MAKE NODES TABLES IF NEEDED     #
##################################################

sub MakeNodesTable
{
  $WorkingDir = $_[0];                             # Source directory of the file to be uploaded
  $TableName = $_[1];                              # The name of the table to be created
  $TableMakeCommand = $_[2];                       # Command used to create the table
  $dbh = $_[3];

  #DELETE THE TABLE IF IT ALREADY EXISTS
  if (&does_table_exist( $dbh, $TableName ))
  {
    print NCBILOG "The table ".$TableName." already exists and must be deleted ..... ";
    $dbh->do("DROP TABLE ".$TableName);
    print NCBILOG "completed.", "\n";
  }

  print NCBILOG "Making the table ".$TableName." ....";
  $dbh->do($TableMakeCommand);
  print NCBILOG ".. completed \n";

}

##################################################
# SUBFUNCTION TO MAKE THE TABLE Gi_taxid         #
##################################################

sub MakeGiTable
{
  $WorkingDir = $_[0];                             # Source directory of the file to be uploaded
  $TableName = $_[1];                              # The name of the table to be created
  $TableMakeCommand = $_[2];                       # Command used to create the table
  $TableMakeTooCommand = $_[3];
  $TableLoadCommand = $_[4];                       # Command used to load data into the table
  $dbh = $_[5];

  #DELETE THE TABLE IF IT ALREADY EXISTS
  if (&does_table_exist( $dbh, $TableName ))
    {
    print NCBILOG "The table ".$TableName." already exists and must be deleted .....";
    $dbh->do("DROP TABLE ".$TableName);
    print NCBILOG "completed.", "\n";
  }

  print NCBILOG "Making the table ".$TableName." ....";
  $dbh->do($TableMakeCommand);
  $dbh->do($TableMakeTooCommand);
  $dbh->do($TableLoadCommand);
  print NCBILOG ".. completed \n";
}

##################################################
# SUBROUTINE TO CHECK TO SEE IF THE              #
# MYSQL TABLE ALREADY EXISTS                     #
##################################################
sub does_table_exist
{
  my ($dbh,$whichtable) = @_;
  my ($table,@alltables,$found);
  @alltables = $dbh->tables();
  $found = 0;
  foreach $table (@alltables) 
    {
      $found=1 if ($table eq "`".$whichtable."`");
    }
  return $found;
}

##################################################
# SUBROUTINE TO SEE HOW MANY RECORDS             #
# EXIST IN THE MYSQL TABLE                       #
##################################################
sub how_many_records
{
  my ($dbh,$whichtable) = @_;
  my ($result,$cur,$sql,@row);
  
  $sql = "select count(*) from $whichtable";
  $cur = $dbh->prepare($sql);
  $cur->execute();
  @row=$cur->fetchrow;
  $result=$row[0];
  $cur->finish;
  return $result;
}

##################################################
# SUBROUTINE TO GET THE CURRENT TIME IN A        #
# USABLE FORM                                    #
##################################################

sub GetTime
{
  ($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time);
  
  my $RealMonth = $Month + 1;
  
  if($RealMonth < 10)
    {
      $RealMonth = "0" . $RealMonth; # add a leading zero to one-digit months
    }
  
  if($Day < 10)
    {
      $Day = "0" . $Day; # add a leading zero to one-digit days
    }
  
  if ($Hour < 12) {$ampm = "AM"}
  else {$Hour = $Hour - 12; $ampm = "PM"}
  
  if ($Hour == 0) {$Hour = $Hour + 12}
  
  $Year = $Year + 1900;
  return $RealMonth, $Day, $Year, $Hour, $Minute, $Second, $ampm, $DayOfYear, $IsDST, $WeekDay
    
}
