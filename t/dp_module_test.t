#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_module_test.t - Are required modules present           |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 12/18/2007                                       |
# UPDATED: 12/18/2007                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test to see if modules required by DAWG-PAWS are         |
#  present.                                                 |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Test::More tests => 10;

print "Testing required modules ...\n";

#-----------------------------+
# Getopt::Long                |
#-----------------------------+
# Used for command line interface
use_ok( 'Getopt::Long' ); 

#-----------------------------+
# OS INDEPENDENCE             |
#-----------------------------+
use_ok( 'File::Copy' );        # Copy files without system cmds
use_ok( 'Cwd' );               # Get the current working directory

#-----------------------------+
# MODULES FOR print_help      |
#-----------------------------+
# Modules used for print help subfunction 
use_ok( 'Pod::Select' );       # Print subsections of POD documentation
use_ok( 'Pod::Text' );         # Print POD doc as formatted text file
use_ok( 'IO::Scalar' );
use_ok( 'IO::Pipe' );
use_ok( 'File::Spec' );

#-----------------------------+
# BIOPERL                     |
#-----------------------------+
# Bioperl packages
use_ok( 'Bio::SearchIO' ); 
use_ok( 'Bio::SeqIO' ); 


