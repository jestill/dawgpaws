DAWGPAWS
========

Description
-----------

DAWGPAWS is a tool for discovering the location of the genes and transposable elements in eukaryotic genomic sequences. It is distributed as a suite of command line programs that are designed to assist a Distributed Annotation Working Group (DAWG) in the annotation of genomic sequence contigs. These programs generate the multiple tracks of annotation evidence that can be computationally combined or manually curated to generate gene models and richly annotated transposable element predictions. The computation of the annotation evidence tracks can be distributed across nodes in a high-performance computing environment providing a scalability that makes this a Pipeline to Annotate Whole-genome Sequences (PAWS). The flexibility of evidences that can be generated by DAWGPAWS allows it to be applied to any eukaryotic genome annotation effort.


Features
--------

DAWGPAWS uses a command line interface for all programs. The command line interface was chosen since computational annotation processes will usually be performed remotely on high-performance computing clusters. Since a standardized command line interface has been written for all programs, no prior knowledge of Perl is required to use this suite.

The current features of DAWGPAWS include:

* FASTA file manipulation to prepare sequences for the annotation pipeline
* Annotation of gaps in the sequence assemblies
* Executing ab initio gene annotation and transposable element annotation programs in a high throughput pipeline
* BLAST pipeline suitable for use in a cluster computing framework
* Conversion of output from annotation software to the GFF file format
* GFF2 and GFF3 are both supported output formats
* Conversion of GFF files to the game.xml format
* Help is available in all command line programs using command line flags. Help info and full manual files are still under development for some programs. To access help where progname is the name of the individual program:
  * progname --usage 
    * generate a basic usage message with info on the required arguments
  * progname --help
    * generate a more extensive help message that includes required arguments and all options
  * progname --man
    * view the full manual for the program

User Manual and Software Documentation
---------------------------------------

A comprehensive user manual describing the use of DAWGPAWS for genome annotation is available from the DAWGPAWS web site on SourceForge. This manual is also included in the release 1.0 version of the DAWGPAWS program.

How to Cite DAWGPAWS
---------------------

If you use DAWGPAWS in your annotation efforts, you should cite the DAWGPAWS manuscript published in Plant Methods:

>JC Estill and JL Bennetzen. 2009. "The DAWGPAWS pipeline for the annotation of genes and transposable elements in plant genomes." Plant Methods. 5:8.

You can download the open access pdf of at [http://www.plantmethods.com/content/5/1/8](link.address.here).