// kdChainer.h: interface for the kdChainer class.
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

#if !defined(AFX_KDChainer_H__9BF172B0_4AB3_4589_AE5C_027196D93625__INCLUDED_)
#define AFX_KDChainer_H__9BF172B0_4AB3_4589_AE5C_027196D93625__INCLUDED_



//#include "chregion.h"

//#include "match.h"
#include "ProcessFragmentFile.h"
#include <stdio.h>
#include <stdlib.h>
#include "kdTree.h"
#include "ContigSorting.h"


#define MaxFragmentLine 4999

#ifndef KDDIMENSION
#define KDDIMENSION MAXSEQUENCENUMBER
#endif





class kdChainer  
{
public:
        kdChainer();
	kdChainer(long* );
        ~kdChainer();
/////////////////////////////////Data////////////////////////////
	long* genome_sizes_array;
	double MIN_VALUE;
	int chain_depth;
	double chain_score;
	double multiplyingfactor;
	int shiftvalue;
	long gap_threshold;
	char gap_threshold_flag;
	char cluster_flag;
	kdTree kdtreeobj;

	
	
	Kdelem_ptr *perm;    //stores the elements of the tree in order
	Kdelem_ptr nntarget, first, tail; 
	Kdnode_ptr *bucket_ptr;  //stores the bucket location of each element
	//Kdnode_ptr root;
	Kdelem_ptr blocks;   //permanent storage of the elements
	long* sorted_list;
	char* point_type;	
	ContigSorting* contigobj; // for processing EST and Contigs
	long cutoff;     //bucket size
	double nndist;
	long  nnptnum; 
	long dimensions;
	long numofblocks;
	double global_maxscore; // max score of all local chain or a global chain
	double chain_average_length;
	double maxscore_gapti;
	double max_fragment_weight;
	double min_connect;
	long min_match_length;
	struct child* chainroots;
	double percentagescore;
	int perc_cover_flag;
	int ref_pos_flag;
	int stdout_flag;
	int chainer_format_flag;
	int absolute_flag;
	int cluster_file_flag;
	long overlapping_value;
	int overlapping_side; // 1 =overlap in first sequence not in others, 0 overlap in all
	int detect_repeat_flag;
	/////////////////////// Reading Fragment file Functions/////////////

	long kd_read_fragment_file(char* filename);
	int analyzefragmentline (FILE *fp, unsigned int mline,Kdelem* elem,
				 long *max_dimensions);
	int transform_for_overlapping(Kdelem* elem);
	int inverse_transform_for_overlapping(Kdelem* elem);
	void gapped_kdelem(int gap_type);

        /////////////Sorting and Initialization Functions////////////////      			
	int quick_sort_list(long in_dimension);
	void recursive_subdim_sort(long in_dimension, long startkeys, long endkeys);
	long get_dimension_value(long dim, long index, long type_index);
	void quicksortkdlistkeysOndim(long dim,long* a,long lo, long hi);
	void initialize_blocks(char* filename);


 public:
        ///////////////////// Gnerating Global Chain ///////////////
	// IncludeDirction and includeWeight are obsolete
	int kdtree_mst();
	int KdtreebasedModgetchain(ProcessFragmentFile* InfoObj,int IncludeWeight,
				    int IncludeDirection); 
 public:
	/////////////// To handle Gaps in Global Chain ////////////
	
	int kdtree_gapped_chain();
	int GappedKdtreebasedModgetchain(ProcessFragmentFile* InfoObj,int IncludeWeight,
					  int IncludeDirection);


	//////////////////// Functions for Local Chaining //////////////////////
	int SW_chaining(double Threshold);
	int mod_SW_chaining(double Threshold);
	void storeMultiChain(ProcessFragmentFile* InfoObj,int IncludeDirection);
	Kdelem_ptr traverse_store_chain(child* root, FILE* fptr, double*);
	int KdtreeSW(ProcessFragmentFile* InfoObj, int IncludeWeight,
		      int IncludeDirection,double Threshold);
	void free_blocks_array();
	void free_chain_roots();
	int check_chain(Kdelem_ptr,int,double,double );
	/////////////////// Variation for handling Repeats /////////////////////
	int mod_SW_repeat(double Threshold);
	//////////////////// Functions for controlling the ouput chain///////////
	void set_overlapping_value(int in_overlapping_value){overlapping_value=in_overlapping_value;};
	void set_overlapping_side(int in_overlapping_side){overlapping_side=in_overlapping_side;};
	void set_detect_repeat_flag(){detect_repeat_flag=1;};
	void set_chain_depth(int in_depth){chain_depth=in_depth;};
	void set_chain_average_length(double in_avergae_length){chain_average_length=in_avergae_length;};
	void set_chain_score(double in_score){chain_score=in_score;};
	void set_multiplying_factor(double lw_value){multiplyingfactor=lw_value;};
	void set_gap_threshold(long gt_value){gap_threshold=gt_value;};
	void set_gap_threshold_flag(){gap_threshold_flag=1;};
	void set_cluster_flag(){cluster_flag=1;};
	void storeChain(ProcessFragmentFile* InfoObj,int IncludeDirection);
	void set_stdout_flag(int in_stdout_flag){stdout_flag=in_stdout_flag;};
	void set_ref_pos_flag(int in_ref_pos_flag){ref_pos_flag=in_ref_pos_flag;};
	void set_chainer_format_flag(int in_chainer_format_flag)
	    {chainer_format_flag=in_chainer_format_flag;};
	void set_absolute_flag(int in_absolute_flag)
	    {absolute_flag=in_absolute_flag;};
	void set_cluster_file_flag(int in_cluster_file_flag)
	    {cluster_file_flag=in_cluster_file_flag;};

	///////////////////  Functions for Contigs/EST Placement ///////////////////
 
	int ContigESTChainer(ProcessFragmentFile* InfoObj, int IncludeWeight,
			     int IncludeDirection,double Threshold, char* contigfilename,int); 
        
        ///////////////////  Functions for Contigs ///////////////////
	int mod_Contig_SW(double Threshold,char* contigfilename);
	///////////////////  Functions for EST Placement ///////////////////

	void EST_SW(double Threshold,char* contigfilename);
	void set_percentage_threshold(double perc_value)
	    {percentagescore=perc_value;perc_cover_flag=1;};	
	
	int  check_est(int depth, double th_score,double given_score,
		       ContigSorting* contigobj,double percentagescore);
	int mod_EST_SW(double Threshold,char* contigfilename);
	Kdelem_ptr traverse_store_EST(child* root, FILE* fptr, double* maxscoreptr);
	void storeMultiEST(ProcessFragmentFile* InfoObj,int IncludeDirection);
	//////////////////////Storing in database format///////////////7
	void kdChainer::storeChainDB(ProcessFragmentFile* InfoObj,int);
	void kdChainer::storeMultiChainDB(ProcessFragmentFile* InfoObj,int IncludeDirection);
	void kdChainer::storeMultiContigDB(ProcessFragmentFile* InfoObj,int IncludeDirection);
	void kdChainer::storeMultiEstDB(ProcessFragmentFile* InfoObj,int IncludeDirection);
	void kdChainer::silent_chain_traversal();
	void kdChainer::silent_traverse_cluster(child* root, double* maxscoreptr);

};

#endif // !defined(AFX_KDChainer_H__9BF172B0_4AB3_4589_AE5C_027196D93625__INCLUDED_)
