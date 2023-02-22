// ContigSorting.h: interface for the ContigSorting class.
//
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

#if !defined(AFX_CONTIGSORTING_H__6D0DF9A4_2EB4_45E4_8EC0_58ABE6A3C32D__INCLUDED_)
#define AFX_CONTIGSORTING_H__6D0DF9A4_2EB4_45E4_8EC0_58ABE6A3C32D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000




//#include "match.h"
#include "kdTree.h"

class chainblock{
public:
	chainblock();
	virtual ~chainblock();
public:
	long startx,endx;
	long starty,endy;
	long contig_id;
	long new_contig_id;
	double score;
	chainblock*last;
	char chainflag;
};


class kdchainblock{
public:
	kdchainblock();
	virtual ~kdchainblock();

public:
	Region* region;
	long contig_id;
	long new_contig_id;
	double score;
	double weight;
	chainblock* last;
	char chainflag;
};


class contigblock{
public:
	char contig_desc[300];
	char flag;
};

class ContigSorting  
{
public:
	ContigSorting();
	virtual ~ContigSorting();
public:
	///////////////// Attributes ///////////////////////
	long* sorted_list;
	char* point_type;
	long* recordseps;
	long* newrecordseps;
	long* mapingcontigs;
	contigblock* contig_block_array;
	long number_of_chains;
	long number_of_contigs;
//	Multimatch* multimatchtab;
	chainblock* chainblocks;
	long total_conitg_length;
	int perc_cover_flag;
	double MIN_VALUE;
	long newtotal_conitg_length;
	double one_dimensional_maximum_score;
	chainblock* one_dimensional_maxelement;
	long dimensions;
	long* offset_array;
	long* total_contig_length_array;
	long number_of_draft_genomes;
	long* number_of_contigs_array;
	
	///////////////// Implementation ///////////////////
public:
	
	int    qsort_chains();
	void   quicksortkdlistkeysOndim(unsigned char dim,long* a,long lo, long hi);
	double get_dimension_value(unsigned char dim, long index, long type_index);
	int    read_contigs(char* filename,int savingflag);
	long   kd_read_chain_file(char* filename,long dimensions);
	
	int    analyzechainline (FILE *fp, unsigned int mline, Kdelem*, long*);
	int    process_chains(char* filename);
	long   getrecordnum(long *recordseps,long numofrecords,long totalwidth,long position);
	long getcontigsize(long position,long* recordseps);
	void set_perc_flag(){perc_cover_flag=1;};
////////////////// Generalized

	long getrecordnumGeneralized(long* positions,long* contig_ids);
	int read_kd_contigs(char* filename, int savingflag);
	long TransformKdelemtoRefGeneralized(Kdelem* elem,long* contig_id);
	long kd_read_chain_file_kdraftgenome(char* filename,long in_dimensions);
	////////////////////////////////////////////////////////
	//////////////////One Dimesional Chaining Functions
	int TransformKdelemtoRef(Kdelem* elem);
	int specific_qsort_chains();
	void specific_quicksortkdlistkeysOndim(unsigned char dim,long* a,long lo, long hi);
	double specific_get_dimension_value(unsigned char dim, long index, long type_index);
	void one_dimensional_chaining();
	int transform2reference();
	int transform2absolute_sortedcontig();


};

#endif // !defined(AFX_CONTIGSORTING_H__6D0DF9A4_2EB4_45E4_8EC0_58ABE6A3C32D__INCLUDED_)
