#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# HMMER PARSE                                               |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 08/07/2006                                       |
# UPDATED: 08/07/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Parse the results of an HMMSEARCH run to find the        |
#  positions.                                               |
#                                                           |
# USAGE:                                                    |
#  HmmParse.pl                                              |
#                                                           |
#-----------------------------------------------------------+

# Currently just using this as a way to work out the variables
# that are available from the HMMER Parser in BioPerl
# Can run this with a test flag to just run through
# the set to make sure that all of the databases and input
# files exist.

#-----------------------------+
# FROM BIOPERL DOCS           |
#-----------------------------+
# HMMER reports two different scores, a per sequence total score 
# (and evalue) and a per domain score and evalue. This object 
# represents a collection of the same domain with the sequence 
# bits score and evalue. (these attributes are also on the per 
# domain scores, which you can get there).

print "The program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Std;               # Allows to get options from the command line
use Text::Wrap;                # Allows word wrapping and hanging indents
                               # for more readable output for long strings.
#use Bio::SearchIO;
#use '/usr/lib/perl5/site_perl/5.8.0/Bio/Tools/HMMER';
use Bio::Tools::HMMER::Results;

my $FilePath =  '/home/jestill/projects/asgr/hmm_mite/seqs/PPAN058CON001/'.
    'PPAN058CON001_Stowaway3_rc.hmm.hmmout';


# HMM RESULT
my $HmmRes = new Bio::Tools::HMMER::Results ( -file => $FilePath ,
					      -type => 'hmmsearch') ||
    die "Could not open file\n$FilePath\n";

print "Number of Results: ".$HmmRes->number."\n";

foreach $seq ( $HmmRes->each_Set ) 
{
    print "Seq BITS   ".$seq->bits."\n";
    print "Seq EVale  ".$seq->evalue."\n";
    print "Seq Name   ".$seq->name."\n";

    my $Desc = $seq->desc || "No Description";
    print "Seq Name   ".$Desc."\n";

    my $Acc = $seq->accession || "NO ACC";
    print "Seq Acc".$Acc."\n";

    foreach $domain ( $seq->each_Domain ) 
    {
	print " Dom Start ".$domain->start."\n";
	print " Dom End   ".$domain->end."\n";
	print " Dom BITS  ".$domain->bits."\n";
	print " Dom E     ".$domain->evalue."\n";
	print " Dom name   ".$domain->hmmname."\n";
	print " Seq BITS  ".$domain->seqbits."\n";
	# Get the Name Start End. This is useful to 
	# create a unique name
	print " Seq NSE   ".$domain->get_nse."\n";
    }
}

exit;


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/07/2006
# - Program started as a way to figure out the vars that
#   are available for a HMM Parse using BioPerl. Goal is
#   to parse the results of the MITE HMM Runs to figure
#   out how well these can be used to find MITEs in the 
#   data set from Tifton.
