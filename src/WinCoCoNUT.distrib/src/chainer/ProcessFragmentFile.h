// ProcessFragmentFile.h: interface for the ProcessFragmentFile class.
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/*
  Copyright by Mohamed I. Abouelhoda (C) 2003
  =====================================
  You may use, copy and distribute this file freely as long as you
   - do not change the file,
   - leave this copyright notice in the file,
   - do not make any profit with the distribution of this file
   - give credit where credit is due
  You are not allowed to copy or distribute this file otherwise
  The commercial usage and distribution of this file is prohibited
  Please report bugs and suggestions to <mibrahim@informatik.uni-ulm.de>
*/

#if !defined(AFX_PROCESSFRAGMENTFILE_H__84CD68B1_169D_4118_8470_D28E438BAC88__INCLUDED_)
#define AFX_PROCESSFRAGMENTFILE_H__84CD68B1_169D_4118_8470_D28E438BAC88__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "kdTree.h"
#include "myglobaldef.h"
#include <stdio.h>


/*
The class ProcessFragmentFile is responsible for parsing
the fragment file in the chainer format. It contains a function that parses  
the fragment file and computes the statistics file if not computed before. It
can also read the genomes sizes from extra log file.  
*/


class ProcessFragmentFile  
{
public:
	//////////////////// Attributes /////////////////////////////////
	char filename[500];
	long no_of_fragments;
	long no_of_genomes;
	long max_genome_size;
	long max_fragment_length;
	long min_fragment_length;
	long max_fragment_weight;
	long genome_size_array[MaxGenomes];
	long max_locations_array[MaxGenomes]; // max_f{f.x_r} maximum indicies
	
	//////////////////// Operations ////////////////////////////////
 public:
	ProcessFragmentFile();
	~ProcessFragmentFile();
public:	
	int GetFragmentFileStatisticsFormat2(char *infilename); 
	int  analyzefragmentline (char *linearr,struct Kdelem* elem);
	
	//int analyzefragmentline (FILE *fp, unsigned int mline,struct Kdelem* elem,long *max_dimensions);

};

#endif // !defined(AFX_PROCESSFRAGMENTFILE_H__84CD68B1_169D_4118_8470_D28E438BAC88__INCLUDED_)






