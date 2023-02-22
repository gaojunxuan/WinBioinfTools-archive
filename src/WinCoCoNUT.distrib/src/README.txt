*********************************************************
*///////////////////////////////////////////////////////*
*               This is CHAINER version 1.0             *
*                                                       *
*        Copyright by Mohamed I. Abouelhoda  (C) 2004   *
*                                                       *
*       Unauthorized commercial usage and distribution  *
*              of this program is prohibited.           *
*             Contact the author for a license.         *
*           Please report bugs and suggestions to       *
*             <mibrahim@informatik.uni-ulm.de>          *
*///////////////////////////////////////////////////////*
*********************************************************

INTRODUCTION

CHAINER implements chaining algorithms for comparing genomes.
The source code is based on the algorithms developed by
Abouelhoda-Ohlebusch (see references in the manual). The programming
languages is standard C++. The input/output functions are from
the "stdio" library of standard C. 
The distributed version of CHAINER includes also shell scripts
demonstrating its usage and also visualization modules written
in Perl by Kathrin Hockel. It also includes test dataset.


LICENCE AGREEMENT
There is no distribution fee for CHAINER for non commercial research
institutions. Commercial research institutions who want to use the program
should directly contact the author
		Mohamed I. Abouelhoda 
		<mibrahim@informatik.uni-ulm.de>. 
Unauthorized commercial usage and  distribution of this program is 
prohibited.  Please, report bugs and suggestions to the author.
Refer to the CHAINER publications whenever CHAINER used in your
research; see reference list in the manual or licence agreement. 
The complete licence agreement is in the same directoy.

DIRECTORY STRUCTURE

There are 5 directories with this version of chainer: /docs, /scripts, /src, 
/bin, and /testdata. 

The directory /docs contains this README file and a PS/PDF file 
containing CHAINER manual and tutorial.

The directory /bin contains the executable codes for CHAINER and the
format-transform program match2chainer. These binaries are generated
for Linux and for SUN-solaris. There is 64-bit version for 
Sun-solaris.

The directory /src contains the source code for chainer.

The directory /testdata contains a set of test data to demonstrate the usage
of CHAINER. This directory is compressed and can be downloaded from CHAINER website.

DOCUMENTATION

Please read the user manual to see how to run CHAINER on the data set.
This documentation is not only a manual presenting the options of CHAINER,
but it is also a tutorial that discusses how to efficiently compare genomes and
other related tasks such as mapping cDNA/EST database to a genomic sequence.


BINARY FILES

The binary files are compiled for the Linux and Sun-solaris platforms.
There is also a 64-bit version for Sun-solaris. If you have another machine, 
then you can go to the source code directory, run a Makefile with the correct 
(check the libraries and comiler paths against your system) options, and 
generate the binaries. For running the 64-bit version, be sure that the global 
settings on your machines are correct.

This binaries are restricted to work on 128 genomes. Any line in the fragment
file must not exceed 4998 characters. The maximum length of a file name must 
not exceed 1498 characters. The maximum fragment position in a genome must be 
less than 2x10^9 for 32-bits machines.  

This directory contains perl scripts for visualizing the output of chainer
using gnuplot.



SOURCE CODE

The source code is in the /src directory. Please do not make any changes, read 
the license agreement. For suggestions and bugs, please contact the author.
In this directory there is a Makefile that generates the binary files
for CHAINER and for the program match2chainer. If you want to generate 
binaries for 64-bit version, run the make file with the argument 64bits, 
i.e., write "make 64bits".  Be sure that you compiler is configured 
appropriately. Also check that the paths in the Makefile is the same as in your
system. Please also see the global definitions file "myglobaldef.h", 
as it contains setting for parsing files. It also contains the maximum and 
minimum values for fragment weights, maximum length of a file name, 
the maximum length of a line in a fragment file.


TEST DATA

The directory /testdata contain simple fragment files that are used
in the examples of the tutorial. It also contain 3 subdirectories:
\bacteria, \draftgenome, and \EST. These contain example genomes
and fragment files used in the examples of the manual 
that run the shell scripts.
The directory \bacteria contain the directories
\chlamydia_pneumoniae,  \chlamydia_trachomatis, and  \Staph4.
The first two directories contain fasta files
for the bacteria C. pneumoniae and C. trachomatis.
The third one contains strains of the bacteria Staphylococcus.
The directory /draftgenome contains the draft genome
M.bovis BCG Pasteur stored in the file BCG BCGdraft.dbs.
It contains the complete genome M. bovis MB.dbs.
The directory /EST contains a subsdirectory /mouse_test
that contains a fragment file generated between
a database of mouse ESTs and a part of the mouse genome.
This directory contains the description file "matches.seqinfo"
that is necessary for the chaining.



SHELL SCRIPTS

These shell scripts demonstrate how to combine chainer with other 
software tools and scripts to perform comparative genomic tasks.

We assume that CHAINER, the fragment-generation programs, and 
the visualization scripts are stored in the \bin directory. 
For running these scripts, we need the programs mkvtree, 
vmatch, vseqinfo, and multimat. These programs can be obtained 
from "www.vmatch.de". We also assume that gnuplot and perl are 
already installed on your system.

There are 6 scripts: 
comp2genomes2.sh, 
comp_multiple_genomes.sh,
global_comp_multiple_genomes.sh,
comp_draft_finished.sh,
map_genome_est.sh,
splicing_map_genome_est.sh. For details on these scripts, see the manual.






