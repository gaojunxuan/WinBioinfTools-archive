// kdTree.cpp: implementation of the kdTree class.
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



#include "kdTree.h"
#include "myglobaldef.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif



////////////////////My Declarations////////////////////////////////////
#define NUMOFCELLS 512

#define ERRNOTOPEN      "cannot open matchfile"
#define ERREMPTY        "empty matchfile"
#define ERRINCONSISTENT "inconsistent number of positions in line"
#define ERRTOOFEW       "too few positions in line"
#define FAIL 0

#ifndef max
#define max(a, b)  ((a)>=(b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b)  ((a)<=(b) ? (a) : (b))
#endif
#define px(m,d)  (perm[m]->region[d].end)

#define exchange(i,j)  { Kdelem_ptr buf = perm[i];\
                         perm[i] = perm[j]; perm[j] = buf; }

#define GETREGIONLENGTH(R)\
  ((R).end - (R).start + 1)

#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define NaturalLog(X) 2.3025*log(X)




//////////////////////////////////////////////////////////////////////


dimensions_flag ::dimensions_flag()
{

	for (int i=0;i<KDDIMENSION;i++){
		flag_array[i]=0;
	}

}

dimensions_flag ::~dimensions_flag()
{

}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

kdTree::kdTree()
{

		//genome_sizes_array=genome_sizes;
	cutoff = 10;     //bucket size 
	dimensions=dimensions;
	numofblocks = 0;
	min_connect =(double)0.00;
	sorted_list=NULL;
	point_type=NULL;
	root=NULL;
	MIN_VALUE=(double)-2000000000;
	chainroots=NULL;
}

kdTree::kdTree(long dimensions)
{

	//genome_sizes_array=genome_sizes;
	cutoff = 10;     //bucket size 
	dimensions=dimensions;
	numofblocks = 0;
	min_connect =(double)0.00;
	sorted_list=NULL;
	point_type=NULL;
	root=NULL;
	MIN_VALUE=(double)-2000000000;
	chainroots=NULL;

}

kdTree::~kdTree()
{
	if(root!=NULL)
	this->DeletekdTree(root);
}

//FUNCTIONALITY FOR PRINTING OUT INFORMATION

/*
  shows the elements in perm from perm[l] to perm[r]
*/
void kdTree::showElems(long l, long r){
  long i=l,j;
  for(;i<=r;i++){
    printf("%lu: ",i);
    for(j=0;j<dimensions;j++){
      printf("%lu  %lu  ||  ",perm[i]->region[j].start,perm[i]->region[j].end);
    }
    printf("\n");
  }
  printf("\n");
}

/*
  shows the tree (incl. subtrees) from node start
*/
void kdTree::showKdTree(Kdnode_ptr start){
  if(start->bucket == 1){
    showElems(start->lopt,start->hipt);
    return;
  }
  else{
    if(!start->empty){
      printf("internal node\n");
      showKdTree(start->loson);
      showKdTree(start->hison);
    }
  }
}

/*
  shows the scores of all nodes from node start on with their depth
*/
void kdTree::showScores(Kdnode_ptr start, long depth){
  if(start->bucket)
    {
      printf("Depth: %lu, Max_score: %f\n",depth,(float)start->max_score);
    }
  else if(!start->empty){
    printf("Depth: %lu, Max_score: %f\n",depth,(float)start->max_score);
    showScores(start->loson, depth+1);

    showScores(start->hison, depth+1);
  }
}
  
/*
  prints out the computed chain
*/
void kdTree::showChain()
{
  Kdelem_ptr elem = tail;  //last element of the chain
  long i;
  FILE* fptr;
  //Kdelem_ptr elem = tail;  //last element of the chain
  printf("Score:%f",tail->globalscore);
  fptr=FOPEN("showchain.txt","w");
  for(; elem; elem = elem->last)  //last references the previous elem
    {
      for(i=0;i<dimensions;i++)
	{
	  fprintf(fptr,"[%lu,%lu] ", elem->region[i].start,elem->region[i].end);
	}
      fprintf(fptr,"\n");
    }
  fclose(fptr);
}

//////////////////////////////////////////////////////////////////////
//-----FUNCTIONALITY RELATED TO CALCULATING THE MEDIAN VALUE--------//
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
//Function: randomized_partition
/*
  sorts perm in between l and u according to x = perm[(u+l)/2]
  returns index of the rightmost position of x
*/
//////////////////////////////////////////////////////////////////////

// The new version from Janina
long kdTree::randomized_partition(long l, long u, long dim)
{
  long lp = l;
  long up = u-1;
  long i = (u+l)/2;   //oder random-Wert?
  long x = perm[i]->region[dim].end;
  exchange(u,i);
  while(1)
    {
      while (perm[lp]->region[dim].end < x) lp++;
      while (perm[up]->region[dim].end > x && up > l) up--;
      if(lp < up){
	exchange(lp,up);
	lp++;
	up--;
      }
      else{
	exchange(lp,u);
	return lp;
      }
    }
}


//////////////////////////////////////////////////////////////////////
//Function: randomized_select
/*
  calls function randomized_partition until position m is returned
  perm[l] to perm[m] contains only elements <= perm[m]
  perm[m] to perm[u] contains only elements >= perm[m]
*/
//////////////////////////////////////////////////////////////////////

void kdTree::randomized_select(long l, long u, long m, long dim)
{
  long q;
  if(l == u) return;
  q = randomized_partition(l,u,dim);
  if(q == m) return;
  if(q > m) randomized_select(l, q-1, m, dim);
  else randomized_select(q+1, u, m, dim);
}

//////////////////////////////////////////////////////////////////////
//------FUNCTIONALITY FOR CONSTRUCTING AND USING A KD-TREE-----------
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//Function: build  
/*
  builds a kd-tree for elems of perm from perm[l] to perm[u] 
  using dimensin as the cut dimension
  using the median as the cut value
  function recursively calls itself with new dimension bounds
*/
//////////////////////////////////////////////////////////////////////
Kdnode_ptr kdTree::build(long l, long u, long cutdim){
    Kdnode_ptr tree_root;
    status=1;
    tree_root=build_tree(l, u,cutdim);
    if((status==FAIL)&&(tree_root!=NULL)){
	this->DeletekdTree(tree_root);
	tree_root=NULL;
    }
    return tree_root;
}

Kdnode_ptr kdTree::build_tree(long l, long u, long cutdim)
{
  Kdnode_ptr p;
  long m,i;
  long nextdim = (cutdim+1)%dimensions; //looping over all dimensions

  p = (Kdnode_ptr) malloc(sizeof (Kdnode)); //the new node in the tree
  if(p==NULL){
      printf("not enough memory available\n");
      status=FAIL;
      return(NULL); 
  }
  p->max_score=(double)MIN_VALUE-(double)20;
  p->hison=NULL;
  p->loson=NULL;
  if(p == NULL) 
  {
    printf("not enough memory available\n");
    exit(0);    
  }
  p->empty = 0;
  if (u-l + 1 <= cutoff) {  //all elements fit within p
    p->bucket = 1;
    p->lopt = l;
    p->hipt = u;
    for(i=l;i<=u;i++){
      bucket_ptr[i] = p;  //setting references for elements of p
    }                     //->direct access!
  }
  else{
    p->bucket = 0;
    p->cutdim = cutdim;
    m = (l+u)/2;         //postulate position of median value
    randomized_select(l,u,m,p->cutdim); //sort perm -> median at pos m
    p->cutval = px(m, p->cutdim);
    p->loson = build(l, m, nextdim);    //recursive call for elems <= cutval
    p->hison = build(m+1, u, nextdim);  //recursive call for elems >= cutval
    p->loson->father = p;
    p->hison->father = p;
   }
   return p;
}

//////////////////////////////////////////////////////////////////////
// Function: swapelems
/*
  enables exchanging of elements without their numbers being exchanged
  -> necessary for deleting and undeleting of elements
*/
//////////////////////////////////////////////////////////////////////

void kdTree::swapelems(long i, long j)
{
  Kdelem_ptr buf = perm[i];
  perm[i] = perm[j];
  perm[j] = buf;
  //perm[i]->num = i;
  //perm[j]->num = j;
}

//////////////////////////////////////////////////////////////////////
// Function: isroot
// checks if a node p is the tree root or not
//////////////////////////////////////////////////////////////////////

long kdTree::isroot(Kdnode_ptr p){
  if(p == root){
    return 1;
  }
  else return 0;
}



//////////////////////////////////////////////////////////////////////
// Function:DeleteElem
/*
  deletes an element from the tree by moving it to the end of its bucket
  and decreasing the right index of the bucket node into perm
*/
////////////////////////////////////////////////////////////////////

void kdTree::DeleteElem(Kdelem *elem)
{
  Kdnode_ptr p = bucket_ptr[elem->num];
  long j = p->lopt;
  while(perm[j]->num != elem->num) j++; //finding the element to be deleted
  swapelems(j, p->hipt);
  p->hipt--;
  if(p->lopt > p->hipt){
    p->empty = 1;
    while (!isroot(p) && (p = p->father) && p->loson->empty && p->hison->empty){
      p->empty = 1; //father nodes are empty as well
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Function: DeleteAll
// deletes all the tree elements (deactivate them)
//////////////////////////////////////////////////////////////////////

void kdTree::DeleteAll()
{
  long i = 0;
  for(;i <(numofblocks); i++)
    {
      DeleteElem(perm[i]);
    }
  return;
}


//////////////////////////////////////////////////////////////////////
// Function: undelete
// Activate an element in the tree
/*
  reverses the delete operation by increasing the index of the appropriate
  node into perm and placing the element on that new last position
*/
//////////////////////////////////////////////////////////////////////


void kdTree::undelete(Kdelem *elem)
{
  Kdnode_ptr p = bucket_ptr[elem->num];
  long j = p->lopt;
  double k;
  while(perm[j]->num != elem->num)
    {
      j++;
    }
  p->hipt++;
  swapelems(j, p->hipt);
  p->empty = 0;
  p->max_score = max(p->max_score, elem->globalscore);
  k = p->max_score;
  while(!isroot(p) && (p = p->father) && (p->empty || p->max_score < k)){
    p->empty = 0;
    if(p->max_score < k) p->max_score = k;
  }
}

//////////////////////////////////////////////////////////////////////
// Function: DeletekdTree
// frees the memory occupied by the kdtree
//////////////////////////////////////////////////////////////////////
void  kdTree::DeletekdTree(Kdnode* node){
	 if(node->loson!=NULL) DeletekdTree(node->loson);
	 if(node->hison!=NULL) DeletekdTree(node->hison);
	 free(node);
	 return;

}


//////////////////////////////////////////////////////////////////////
// Function: inorder
/*
  returns 1 when the blocks do not overlap, else 0
*/
//////////////////////////////////////////////////////////////////////

long kdTree::inorder(Kdelem *from, Kdelem *to){
  long i=0;
  for(;i<dimensions+1;i++)
    {
	//printf("%d\n",dimensions);
	//I added dimensions+1 to guarantee colineary even in the first dimension (done in 22.08.2004)
      if(to->region[i].start <= from->region[i].end) return 0;
    }
  return 1;
}


//////////////////////////////////////////////////////////////////////
// Function: Connect
//  returns cost for connecting elems from and to
//////////////////////////////////////////////////////////////////////
// obsolete

double kdTree::Connect(Kdelem *from, Kdelem *to)
{
  //  long i=0;
  //long maxconnect = 0;
  //for(;i<dimensions;i++)
  //  {
  //    if((to->start[i] - from->cval[i])>maxconnect)
  //	maxconnect = to->start[i] - from->cval[i];
  //  }
  //  return maxconnect;
  return 0;
}


//////////////////////////////////////////////////////////////////////
// Function: dist
/*
  returns the cost of connecting from and to minus the score
  of the best chain ending in from
  -> the smaller the result, the better
*/
//////////////////////////////////////////////////////////////////////

double kdTree::dist(Kdelem_ptr from, Kdelem_ptr to)
{
  if(!inorder(from,to)) return(0);
  else return (Connect(from,to) - from->globalscore);
}

//////////////////////////////////////////////////////////////////////
// Function: recursiveNN
/*
  performs a recursive nearest neighbor search for nntarget 
  (starting at node p)
  nearest neighbor = the elem that leads to the best chain
  ending with nntarget!
*/
//////////////////////////////////////////////////////////////////////

void kdTree::recursiveNN(Kdnode_ptr p)
{
  long i, val, thisx;
  double thisdist;
  if(p->empty) return;
  if(p->bucket)
    {
    for(i=p->lopt; i<=p->hipt; i++)
	{
	  thisdist = dist(perm[i], nntarget); //trying out all elems of p
	  if(thisdist < nndist)
	    {
	      nndist = thisdist;
	      nnptnum = i;
	    }
	}
    }
  else
    {
      val = p->cutval;
      thisx = nntarget->region[p->cutdim].start;
      if(thisx < val)  //only in left subtree are appropriate end points
	{
	  recursiveNN(p->loson);
	}
      else
	{
	  recursiveNN(p->hison);
	  //if p's max_score is smaller than the best chain so far
	  //there is no need to look any further
	  if((p->max_score - min_connect) > (0 - nndist))  
	    {
	      recursiveNN(p->loson);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////////////
// Function: nearestNeighbor
//  calls the recursive nearest neighbor function and sets nntarget to j
// based on the kdtree implementation by Bently
/////////////////////////////////////////////////////////////////////

long kdTree::nearestNeighbor(Kdelem *j)
{
  nntarget = j;
  nndist = 0;
  nnptnum = -1;
  recursiveNN(root);
  return nnptnum;
}




//////////////////////////////////////////////////////////////////////
// -------Functionality for Constructing Gapped Global Chain--------//
//////////////////////////////////////////////////////////////////////
// The following functionalities are the same as discribed before
// but considers the boundary conditions required when the gaps are
// included
//////////////////////////////////////////////////////////////////////
// Function: gapped_undelete
// Activate an element in the kdtree in case of using gap costs
//////////////////////////////////////////////////////////////////////

void kdTree::gapped_undelete(Kdelem *elem)
{
    Kdnode_ptr p = bucket_ptr[elem->num];
    long j = p->lopt;
    double k;
    while(perm[j]->num != elem->num)
    {
	j++;
    }
    p->hipt++;
    swapelems(j, p->hipt);
    p->empty = 0;
    p->max_score = max(p->max_score, elem->globalscore);
    k = p->max_score;
    while(!isroot(p) && (p = p->father) && (p->empty|| p->max_score < k)){
	p->empty = 0;
	if(p->max_score < k) p->max_score = k;
    }
}




//////////////////////////////////////////////////////////////////////
// Function: gapped_dist
// The same as the function dist but modified to return MIN_VALUE 
//////////////////////////////////////////////////////////////////////

double kdTree::gapped_dist(Kdelem_ptr from, Kdelem_ptr to)
{
  if(!inorder(from,to)) return(MIN_VALUE-(double)22.00);
  else return (from->globalscore);
}

//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// Function: gapped_recursiveNN
// returns the point with the hieghest score including the gap cost
//////////////////////////////////////////////////////////////////////
void kdTree::gapped_recursiveNN(Kdnode_ptr p)
{
  long i; 
  double thisdist;
  long val, thisx;
  if(p->empty) return;
  if(p->bucket)
    {
    for(i=p->lopt; i<=p->hipt; i++)
    {
	thisdist = gapped_dist(perm[i], nntarget); //trying out all elems of p
	if(thisdist > nndist)
	{
	    nndist = thisdist;
	    nnptnum = i;
	}
    }
    }
  else
  {
      val = p->cutval;
      thisx = nntarget->region[p->cutdim].start;
      if(thisx < val)  //only in left subtree are appropriate end points
      {
	  gapped_recursiveNN(p->loson);
      }
      else
      {
	  gapped_recursiveNN(p->hison);
	  //if p's max_score is smaller than the best chain so far
	  //there is no need to look any further
	  if((p->max_score-min_connect) > (nndist))  
	  {
	      gapped_recursiveNN(p->loson);
	  }
      }
  }
}

//////////////////////////////////////////////////////////////////////
// Function: The same as nearestNeighbor but considers the MIN_VALUE
//////////////////////////////////////////////////////////////////////


long kdTree::gapped_nearestNeighbor(Kdelem *j)
{
  nntarget = j;
  nndist = MIN_VALUE-(double)50;
  nnptnum = -1;
  gapped_recursiveNN(root);
  return nnptnum;
}

////////////////////////////RMQ in Rectangular Area////////////////////
// This RMQ takes as input a rectangle and fins the point with
// maximum score in this rectangle
///////////////////////////////////////////////////////////////////////

double kdTree::region_gapped_dist(Kdelem_ptr inpucket, Kdelem_ptr p, Kdelem_ptr lower_bound)
{
    long i=0;
    //printf("Dimsneions: %d \n",dimensions);
    for(i=0;i<dimensions;i++){
	if((p->region[i].start <= inpucket->region[i].end)||
	   (lower_bound->region[i].end > inpucket->region[i].end)){
	    return(MIN_VALUE-(double)100.00) ;
	}
    }
    return (inpucket->globalscore);
    
}
/////////////////////////////////////////////////////////////////////////////////////
// This function to locat the subtrees above the lower bounds
// Then for each of these subtrees it called the 
// function region_gapped_recursiveNN to get the RMQ in this region

void kdTree::lowerbound_region_gapped_recursiveNN(Kdnode_ptr p, Kdelem_ptr lower_bound,
dimensions_flag on_dim){

    long i; 
    //double thisdist;
    long val, thisx;
    if(p->empty) return;
    if(p->bucket)
    {
	region_gapped_recursiveNN(p,lower_bound);
    }
    else{
	val=p->cutval;
	thisx=lower_bound->region[p->cutdim].start;
	if(thisx > val)  //only in left subtree are appropriate end points
	{
	    lowerbound_region_gapped_recursiveNN(p->hison,lower_bound, on_dim);
	}
	else
	{
	    lowerbound_region_gapped_recursiveNN(p->loson,lower_bound, on_dim);
	    
	    //I need to check that all the cut values are in the range
	    on_dim.flag_array[p->cutdim]=1;
	    
	    // check if all dimensions are included
	    int flag=1;
	    for (i=0;i<dimensions;i++){
		if(on_dim.flag_array[p->cutdim]!=1)flag=0;
	    }
	    if(flag==1){
		region_gapped_recursiveNN(p->hison,lower_bound);
		
	    }
	    else{
		lowerbound_region_gapped_recursiveNN(p->hison,lower_bound, on_dim);
	    }
	}
    }
    
}

void kdTree::region_gapped_recursiveNN(Kdnode_ptr p, Kdelem_ptr lower_bound)
{
    long i; 
    long j;
    double thisdist;
    long val, thisx;
    long dist_i_to_target;
    long dist_nn_to_target=0;
    if(p->empty) return;
    if(p->bucket)
    {
	for(i=p->lopt; i<=p->hipt; i++)
	{
	    thisdist = region_gapped_dist(perm[i], nntarget,lower_bound); //trying out all elems of p
	    if(thisdist >= nndist){
		if(thisdist >nndist){
		    nndist = thisdist;
		    nnptnum = i;
		    // modification to improve EST
		    dist_nn_to_target=0;
		    for(j=0;j<dimensions;j++){
			dist_nn_to_target=dist_nn_to_target+nntarget->region[j].start-perm[i]->region[j].end;
		    }
		    // End of modification
		}
		else{
		    // modification added to improve EST
                    // if two points are with the same score choose the nearst one w.r.t. L1 dist
		    dist_i_to_target=0;
		    for(j=0;j<dimensions;j++){
			// distance to target
			dist_i_to_target=dist_i_to_target+nntarget->region[j].start-perm[i]->region[j].end;
		    }
		    if(dist_i_to_target<dist_nn_to_target){ // update with the nearst one
			nndist = thisdist;
			nnptnum = i;
		    }
		    // end of modification
		}
	    }
	}
    }
    else{
	val=p->cutval;
	thisx=nntarget->region[p->cutdim].start;
	if(thisx < val)  //only in left subtree are appropriate end points
	{
	    region_gapped_recursiveNN(p->loson,lower_bound);
	}
	else
	{
	    region_gapped_recursiveNN(p->hison,lower_bound);
	    //if p's max_score is smaller than the best chain so far
	    //there is no need to look any further
	    
	    // Moh.: Modification 1. check if the point whose score 
	    // annotates P is in the region or not.
	    // Moh.: Modification 2. check for the lower bound
	    if((p->max_score-min_connect) > (nndist))  
	    {
		// If the point annotates P is not in the region go further in the tree
		region_gapped_recursiveNN(p->loson,lower_bound);
		
	    }
	}
    }
}

long kdTree::region_gapped_nearestNeighbor(Kdelem *j, Kdelem_ptr lower_bound)
{
/*    
    Kdelem* j=new Kdelem;
    j->region=(Region*)malloc(sizeof(Region)*(dimensions+1));

    for(long i=0;i<=dimensions;i++){
	j->region[i].start=in_j->region[i].start-1;
	j->region[i].end=in_j->region[i].end-1;
    }
*/   
    nntarget = j;
    nndist = MIN_VALUE-(double)50;
    nnptnum = -1;
    dimensions_flag on_dim;
    lowerbound_region_gapped_recursiveNN(root,lower_bound,on_dim);
/*
    delete j->region;
    delete j;
*/
    return nnptnum;
}

///////////////////////////////////////////////////////////////////////////////////
/*
  deactivate an element from the tree by moving it to the end of its bucket
  and decreasing the right index of the bucket node into perm
  it also adjusts the scores on the ancestor internal nodes
*/
void kdTree::DeactivateElem(Kdelem *elem)
{
  Kdnode_ptr p = bucket_ptr[elem->num];
  long j = p->lopt;
  while(perm[j]->num != elem->num){
      //finding the element to be deleted
      j++;
  }
  swapelems(j, p->hipt);
  p->hipt--;
  if(p->lopt > p->hipt){
      p->empty = 1;
      //p->max_score=MIN_VALUE-(double)50;
      while (!isroot(p) && (p = p->father) && p->loson->empty && p->hison->empty){
	  p->empty = 1; //father nodes are empty as well
      }
  }
  else{
      // updating the maxscore in the tree
      p->max_score=MIN_VALUE-(double)50;
      for(j = p->lopt;j<=p->hipt;j++){
	  p->max_score=max(perm[j]->globalscore,p->max_score);  
      }
      while (!isroot(p) && (p->father!=NULL)){
	  if(p==p->father->loson){
	      p->father->max_score=max(p->father->hison->max_score,p->max_score); 
	  }
	  else{
	      p->father->max_score=max(p->father->loson->max_score,p->max_score); 
	  }
	  p=p->father;
      }
      
  }
}
