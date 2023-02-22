// kdTree.h: interface for the kdTree class.
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


#if !defined(AFX_KDTREE_H__B2911895_D084_4E82_96E5_347BEBB77207__INCLUDED_)
#define AFX_KDTREE_H__B2911895_D084_4E82_96E5_347BEBB77207__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


//#include "match.h"
#include "ProcessFragmentFile.h"
#include <stdio.h>
#include <stdlib.h>



#ifndef KDDIMENSION
#define MAXSEQUENCENUMBER 32
#define KDDIMENSION MAXSEQUENCENUMBER
#endif

typedef struct
{
  long start, end;
} Region;

typedef struct
{
  long length;                          // length of matching sequence
  long* positions;    // start positions of matches
} Multimatch;                           // \Typedef{Multimatch}


typedef struct Kdnode *Kdnode_ptr;

typedef struct Kdnode
{
  long bucket;      //1 = bucket, 0 = internal node
  long empty;       //1 if no blocks in the subtree
  long cutdim;
  long cutval;
  double max_score;   //of a chain ending in its subtree
  Kdnode_ptr loson, hison, father;
  long lopt, hipt;  //left and right index in the block array (perm)
} Kdnode;


typedef struct Kdelem *Kdelem_ptr;

//#define Weighttype double;

typedef struct Kdelem
{
  //Region region[KDDIMENSION];
  Region* region;
  long num;               //used for identification
//  char reversed;
  double globalscore;      //of a chain ending with that element
//  double offsetscore;
// Weighttype weight;     //length of the block
  double weight;           //length of the block
  struct Kdelem* last;
  struct Kdelem* next;     //next is used to loop over all blocks
  struct child* start_child_list;	
} Kdelem;                  //last is used to hold the chain


typedef struct child{
	struct Kdelem* child_elem;
	struct child* next;

}child;

class dimensions_flag
{
public:
	dimensions_flag();
	virtual ~dimensions_flag();
	int flag_array[KDDIMENSION];
};


class kdTree  
{
public:
	kdTree();
	virtual ~kdTree();
	kdTree(long dimensions);
/////////////////////////////////Data////////////////////////////
	long* genome_sizes_array;
	double MIN_VALUE;
	Kdelem_ptr *perm;    //stores the elements of the tree in order
	Kdelem_ptr nntarget, first, tail; 
	Kdnode_ptr *bucket_ptr;  //stores the bucket location of each element
	Kdnode_ptr root;
	Kdelem_ptr blocks;   //permanent storage of the elements
	long* sorted_list;
	char* point_type;
	//ArrayMultimatch multimatchtab;
	long cutoff;     //bucket size
	double nndist;
	long  nnptnum; 
	long dimensions;
	long numofblocks;
	double min_connect;
	long min_match_length;
	struct child* chainroots;
	int status;         // to indicate if the tree correctly built
	//////////////////////////////////////////////////////////////////
	/// Functions
	void showElems(long l, long r);
	void showKdTree(Kdnode_ptr start);
	void showScores(Kdnode_ptr start, long depth);
	void showChain();
	long randomized_partition(long l, long u, long dim);
	void randomized_select(long l, long u, long m, long dim);
	Kdnode_ptr build_tree(long l, long u, long cutdim);
	Kdnode_ptr build(long l, long u, long cutdim);
	void swapelems(long i, long j);
	long isroot(Kdnode_ptr p);
	void DeleteElem(Kdelem *elem);
	void DeleteAll();
	void undelete(Kdelem *elem);
	long inorder(Kdelem *from, Kdelem *to);
	double Connect(Kdelem *from, Kdelem *to);
	double dist(Kdelem_ptr from, Kdelem_ptr to);
	void recursiveNN(Kdnode_ptr p);
	long nearestNeighbor(Kdelem *j);
	void DeletekdTree(Kdnode* node);
	

	void gapped_recursiveNN(Kdnode_ptr p);
	double gapped_dist(Kdelem_ptr from, Kdelem_ptr to);
	long gapped_nearestNeighbor(Kdelem *j);
	void gapped_undelete(Kdelem *elem);

	
	// rectangular region based functions
	void region_gapped_recursiveNN(Kdnode_ptr p,Kdelem_ptr lower_bound);
	double region_gapped_dist(Kdelem_ptr, Kdelem_ptr p, Kdelem_ptr lower_bound);
	long region_gapped_nearestNeighbor(Kdelem *j, Kdelem_ptr lower_bound);
	void lowerbound_region_gapped_recursiveNN(Kdnode_ptr p, Kdelem_ptr lower_bound,
						  dimensions_flag on_dim);
	void DeactivateElem(Kdelem *elem);
	
	
};

#endif // !defined(AFX_KDTREE_H__B2911895_D084_4E82_96E5_347BEBB77207__INCLUDED_)
