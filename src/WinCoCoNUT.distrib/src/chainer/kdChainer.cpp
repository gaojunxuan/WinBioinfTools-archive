// kdChainer.cpp: implementation of the kdChainer class.
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


//#include "stdafx.h"
//#include "ProjMSTree.h"
#include "kdChainer.h"
#include "math.h"
#include "ContigSorting.h"

#include "myglobaldef.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#include <sys/times.h>
#include <unistd.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "string.h"

#include "kdChainer.h"
#include "ProcessFragmentFile.h"


////////////////////My Declarations////////////////////////////////////
//#define NUMOFCELLS 512

#define ERRNOTOPEN      "cannot open matchfile"
#define ERREMPTY        "empty matchfile"
#define ERRINCONSISTENT "inconsistent number of positions in line"
#define ERRTOOFEW       "too few positions in line"

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

#ifndef moduls
#define moduls(a)  ((a)>=(-1*a) ? (a) : (-1*a))
#endif


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

kdChainer::kdChainer()
{
    //multimatchtab=NULL;
    cutoff = 10;     //bucket size   
    numofblocks = 0;
    min_connect = (double)0.00;
    sorted_list=NULL;
    point_type=NULL;
    //root=NULL;
    chainroots=NULL;
    chain_depth=1;
    chain_score=-50000;
    multiplyingfactor=(double)1;
    shiftvalue=3;
    gap_threshold_flag=0;
    gap_threshold=MyMAX_VALUE;
    perc_cover_flag=0;
    percentagescore=(double)0;
    cluster_file_flag=0;
    global_maxscore=MIN_VALUE;
    overlapping_side=0;
    overlapping_value=0;
    detect_repeat_flag=0;
    chain_average_length=0;
    gcfile=NULL;
    detect_repeat_cluster_flag=0;
    cluster_flag=0;
    maxscore_gapti=MIN_VALUE;
}
kdChainer::kdChainer(long * genome_sizes)
{
    //multimatchtab=NULL;
    genome_sizes_array=genome_sizes;    
    cutoff = 10;     //bucket size
    numofblocks = 0;
    min_connect =(double)0.00;
    sorted_list=NULL;
    point_type=NULL;
    //root=NULL;
    MIN_VALUE=(double)MyMIN_VALUE;
   
    chainroots=NULL;
    multiplyingfactor=(double)1;	
    //chain_depth=1;
    //chain_score=-50000;
    shiftvalue=3;
    gap_threshold=MyMAX_VALUE;
    
    gap_threshold_flag=0;
    perc_cover_flag=0;
    percentagescore=(double)0;
    cluster_file_flag=0;
    global_maxscore=MIN_VALUE;
    overlapping_side=0;
    overlapping_value=0;
    detect_repeat_flag=0;
    chain_average_length=0;
    gcfile=NULL;
    cluster_flag=0;
    detect_repeat_cluster_flag=0;
    maxscore_gapti=MIN_VALUE;
}

kdChainer::~kdChainer()
{
    //if(multimatchtab!=NULL)
    //free(multimatchtab);
    if(sorted_list!=NULL)free(sorted_list);
    if(point_type!=NULL)free(point_type);
    // recursive deletion of the tree

    if(bucket_ptr!=NULL)free(bucket_ptr);
    free_blocks_array();
    free_chain_roots();
    if(perm!=NULL)free(perm); 
    
}



////////Functionality For Reading Fragment File in CHAINER format////////

long kdChainer::kd_read_fragment_file(char* filename)
{
    FILE* fptr=NULL;
    
    char tempstr[300];
    long i=0;  
    long j=0;
    int inputchar;
    int flag=0;
    long line = 0;
    long number_of_fragments=0;
    int retval;
    char linearr[MaxFragmentLine+1];
    long max_dimensions;
    long given_dimensions=0;

    max_fragment_weight=MIN_VALUE;

    Kdelem elem;
    fptr=FOPEN(filename,"r"); 
    printf("# reading fragment file %s\n", filename);
    sprintf(tempstr,"%s.ctg",filename);
    
    flag=0;
    while(!feof(fptr)){
 	if(flag==0){
	    flag=1;
	    retval=fscanf(fptr,">CHA %lu\n",&given_dimensions);
	    if(retval<=0){
		printf("Non-Correct Format\n");
		return((long)-1);
	    }
	    continue;
	    
	}
	fgets(linearr,MaxFragmentLine,fptr);
	if(strcmp(linearr,"\n")==0)continue;
	if(strcmp(linearr,"")==0)continue;
	if(linearr[0]=='#')continue;
	if(linearr[0]==0)continue;
	number_of_fragments++;
	strcpy(linearr,"");
	
    }
    //fprintf(outfptr,"\n");
    flag=0;
    fclose(fptr);
    
    // storing the fragments   
    blocks = (Kdelem_ptr) malloc((number_of_fragments+1)* sizeof(Kdelem)+1);
    if(blocks==NULL){
	printf("Not Enough Memory");
	return((long)-1);
    }
    // initializing the ORIGIN point
    blocks[0].region=(Region*)malloc(sizeof(Region)*(dimensions+1));
    if(blocks[0].region==NULL){
	printf("Not Enough Memory \n");
	return((long)-3);
    }
    
    for(i=0;i<=dimensions;i++){
	blocks[0].region[i].start=(long)0;
	blocks[0].region[i].end=(long)0;
	blocks[0].weight=(double)(multiplyingfactor);
	
    }
    // storing the genomes
    elem.region=(Region*)malloc(sizeof(Region)*(given_dimensions));
    if(elem.region==NULL){
	printf("Not Enough Memory \n");
	return ((long)-3);
    }

    fptr=FOPEN(filename,"r"); 
    
    if (fptr == NULL)
    {
	printf("Error opening fragment file\n");
	return((long)-1);
    }
    inputchar = getc (fptr);
    if (inputchar == EOF)
    {
	printf("Error opening fragment file\n");
	return((long)-2);
    }
    ungetc (inputchar, fptr);
    
    
    i=1;
    j=0;
    while (!feof(fptr)){
	line++;
	inputchar = getc (fptr);
	if (inputchar == EOF){
	    break;
	}
	if (inputchar == '#'){	    
	    ungetc (inputchar, fptr);
	    float weight;
	    char strweight[100];
	    retval=fscanf(fptr,"#%s\n",strweight);
	    weight=atof(strweight);
	    if (retval<=0){
		return((long)-1);
	    }	
	    blocks[i].weight=(double)weight*(double)(multiplyingfactor);
	    max_fragment_weight=max(blocks[i].weight,max_fragment_weight);
	}
	else if(inputchar == '>')
	{
	    while ((inputchar = getc (fptr)) != '\n')
		; // Nothing, could store comment but skip it for now
	}
	else if(inputchar == '\n'){
	    continue;
	}
	else if (inputchar == '['){
	    
	    ungetc (inputchar, fptr);
	    max_dimensions=0;
	    retval = analyzefragmentline (fptr,line,&elem,&max_dimensions);
	    if(retval==-1)break;
	    if (retval < 0){
		return((long)-3);
	    }
	    blocks[i].region=(Region*)malloc(sizeof(Region)*(given_dimensions));
	    if(blocks[i].region==NULL){
		printf("Not Enough Memory \n");
		return((long)-3);
	    }
	    blocks[i].region[given_dimensions-1].start=elem.region[0].start+shiftvalue;	   
	    blocks[i].region[given_dimensions-1].end=elem.region[0].end+shiftvalue;
	    // NOTE
            /* the shift value is added to insure that all fragments do not
	       overlap with the ORIGIN fragment
	    */
	    // NOTE
            /*
	      Becuase the range maximum query works on k-1 dimensions.
	      We store the first dimension at the last position in the array region of every 
              fragment (Block) and let the datastructure be built on the remaining dimensions
              starting from 0 to k-1
	     */
	    for(j=0;j<given_dimensions-1;j++){
		// the first dimension was stored at the end
		blocks[i].region[j].start=elem.region[j+1].start+shiftvalue;
		blocks[i].region[j].end=elem.region[j+1].end+shiftvalue;
	    }
	    i++;
	}
	else{
	}
    }
    fclose(fptr);
    free(elem.region);
    if(overlapping_value>0){
	for(i=1;i<number_of_fragments+1;i++){
	    //if(blocks[i].num!=0){	   
	    transform_for_overlapping(&blocks[i]);
	    //}
	}
    }
    return(number_of_fragments+1);
} 
//////////////////////////////////////////////////////////////////////////
// Function: analyzefragmentline
// This function gets a fragment positions from a line in fragment
// file
int kdChainer::analyzefragmentline (FILE *fp, unsigned int mline,Kdelem* elem,long *max_dimensions)
{
    
    
    long i = 0; 
//  long length, position;
    int inputchar;
        
    while (1)
    {
	// bail out at end of line
	inputchar = getc (fp);
	if(inputchar==EOF){
	    ungetc (inputchar, fp);
	    return(-1);
	}
	if (inputchar == '\n'){
	    break;
	}
	
 	if (inputchar == '['){
	    ungetc (inputchar, fp);	    
	    if (fscanf(fp,"[%lu,%lu]",&(elem->region[i].start),&(elem->region[i].end)) != 2)
	    {
		return -1;
	    }	
	    i++;	
	    
	}
	
    }
    if(i>(*max_dimensions)){
	(*max_dimensions)=i;
    }

// make point transformation for certain applications   
//    if(overlapping_value>0)transform_for_overlapping(elem);
//    printf("%u\n",overlapping_side);
    return 0;
}

int kdChainer::transform_for_overlapping(Kdelem* elem){
    int i;
    if(overlapping_side==0){
	for(i=0; i<=dimensions; i++){
	    elem->region[i].start=elem->region[i].start+overlapping_value;
	}
    }
    else{
	if(overlapping_side==1){
	    
	    elem->region[dimensions].start=elem->region[dimensions].start+overlapping_value;	    
	}
	else{
	    elem->region[overlapping_side-2].start=elem->region[overlapping_side-2].start+overlapping_value;
	}
    }
    return (0);
}
int kdChainer::inverse_transform_for_overlapping(Kdelem* elem){
    int i;
    
    if(overlapping_side==0){
	for(i=0; i<=dimensions; i++){
	    elem->region[i].start=elem->region[i].start-overlapping_value;
	}
    }
    else{
	if(overlapping_side==1){
	    elem->region[dimensions].start=elem->region[dimensions].start-overlapping_value;
	}
	else{
	    elem->region[overlapping_side-2].start=elem->region[overlapping_side-2].start-overlapping_value;
	}	
    }
    return (0);
}


////////////////////////////////////////////////////////////////////////////////
// Function: gapped_kdelem
// The function intializes with score relative to the TERMINU point
// Input: gap_type indicates that the initialization is done for
// handling different kinds of gap costs.
// The current version considers gap costs in the rectlinear metric (L1).

void kdChainer::gapped_kdelem(int gap_type)
{
    long i, j;
    

    blocks[0].globalscore=(double)0;
    blocks[0].num=0;	
    for(j=0; j<dimensions; j++)
    {	
	blocks[0].region[j].start = (long)0;
	blocks[0].region[j].end = (blocks[0].region[j].start);
	blocks[0].globalscore=blocks[0].globalscore+(double)(genome_sizes_array[j+1]-blocks[0].region[j].end);
    }
    blocks[0].globalscore=blocks[0].globalscore+(double)(genome_sizes_array[0]-blocks[0].region[dimensions].end);	
    blocks[0].globalscore=(double)blocks[0].weight-blocks[0].globalscore;
    
    
    perm[0] = &blocks[0];
    first = &blocks[0];
    MIN_VALUE=blocks[0].globalscore-(double)10;
    
    for(i=1; i < numofblocks ;i++){
	
	blocks[i].globalscore=(double)0;
	for(j=0; j<dimensions; j++)
	{
	    // initializing gap costs in L1 metric
	    blocks[i].globalscore=blocks[i].globalscore+(double)(genome_sizes_array[j+1]-blocks[i].region[j].end);
	    
	}
	blocks[i].globalscore=blocks[i].globalscore+(double)(genome_sizes_array[0]-blocks[i].region[dimensions].end);	
	blocks[i].globalscore=(double)blocks[i].weight-blocks[i].globalscore;
	blocks[i-1].next = &blocks[i];
	perm[i] = &blocks[i];
    }
}

////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function: quick sort list
// 1- creats two arrays sorted_list and point_type
// the array sorted_list store indices for the sorted start/end points of 
// the fragments and the array point_type store the corresponding 
// point type (start/end).
// 2- the functions sorts the start/end points w.r.t. the first dimension
// which is stored in blocks[i].region[dimension]. If two points have the
// same first dimension value, then they are recursively sorted w.r.t. the
// second dimension and so on.
// the sorted indices are in the two arrays sorted_list with their corresponding
// point type
//////////////////////////////////////////////////////////////////////

int kdChainer::quick_sort_list(long in_dimension){
    
    long i,j;
    if(sorted_list==NULL){
	sorted_list=(long*)malloc(sizeof(long)*(numofblocks*2));
	if(sorted_list==NULL)return(-1);
	point_type=(char*)malloc(sizeof(char)*(numofblocks*2));
	if(point_type==NULL)return(-1);
	if(sorted_list==NULL){
	    return(-1);
	}
	j=0;
	for(i=0;i<(numofblocks*2);i=i+2){
	    sorted_list[i]=j;
	    sorted_list[i+1]=j;
	    point_type[i]=1; // left end
	    point_type[i+1]=2; // right end
	    j++;
	}
    }
    // Sort w.r.t. the first dimension which is stored at the end of
    // array region

    quicksortkdlistkeysOndim(in_dimension,sorted_list,1,(2*numofblocks)-1);
    
    /// sort according to other directions recursive
    recursive_subdim_sort(in_dimension,1,(2*numofblocks)-1);
    return(1);
    
}
////////////////////////////////////////////////////////////////////////////////
// Recursive sorting
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function recursive_subdim_sort
// sorts the indices w.r.t. the next dimension if they are equal w.r.t. the current one
// ////////////////////////////////////////////////////////////////////////////

void kdChainer::recursive_subdim_sort(long in_dimension, long startkeys, long endkeys){
    long i;
    long end=0;
    long start=0;
    unsigned char rec_dimension;
    
    if(in_dimension>dimensions) rec_dimension=1;
    else rec_dimension=in_dimension+1;
    //in_dimension=in_dimension%(dimensions+1)+1;
    char flag=0;
    for(i=startkeys;i<endkeys;i++){	
	if(get_dimension_value(in_dimension,sorted_list[i],i)==get_dimension_value(in_dimension, sorted_list[i+1],(i+1))){
	    if(flag==0){
		start=i;
		flag=1;
	    }
	    end=i+1;
	    if(end<endkeys)continue;
	}
	if((flag==1)&&(end>start)&&(rec_dimension<(dimensions+1))){
	    quicksortkdlistkeysOndim(rec_dimension,sorted_list,start,end);
	    recursive_subdim_sort(rec_dimension,start,end);
	}
	flag=0;
    }
}
////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function quicksortkdlistkeysOndim
// starts sorting w.r.t. the given dimension (dim) between the two array
// boundaries lo and hi
//////////////////////////////////////////////////////////////////////

void kdChainer::quicksortkdlistkeysOndim(long dim,long* a,long lo, long hi)
{
    //  lo ist der unterste Index, hi ist der oberste Index des
    //  zu sortierenden (Teil-)Feldes a
    long i=lo, j=hi;
    //kdMasterPoint temp;
    long temp;
    char tempchar;
    long x=this->get_dimension_value(dim,a[(long)(((double)lo+(double)hi)/2)],(long)(((double)lo+(double)hi)/2));  
    
    //  Aufteilung
    do
    {    
        while (this->get_dimension_value(dim,a[(long)i],i)<x) i++;
        while (this->get_dimension_value(dim,a[(long)j],j)>x) j--;
        
	if (i<=j){
            temp=a[(unsigned int)i]; a[(unsigned int)i]=a[(unsigned int)j]; a[(unsigned int)j]=temp;
	    tempchar=point_type[(unsigned int)i]; 
	    point_type[(unsigned int)i]=point_type[(unsigned int)j]; 
	    point_type[(unsigned int)j]=tempchar;
	    i++; 
	    j--;
        }
    } while (i<=j);
    
    // Rekursion
    if (lo<j) quicksortkdlistkeysOndim(dim,a,lo, j);
    if (i<hi) quicksortkdlistkeysOndim(dim,a,i, hi);    
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function get_dimension_value
//  returns the coordinate of a point referenced by both 
//  the index (in sorted_list array) and its type (in point_type array) according to the
//  given dimension value.
//////////////////////////////////////////////////////////////////////

long kdChainer::get_dimension_value(long dim, long index, long type_index){
    
    // modified to be as blocks data structure
    // Chek if the index is 
    if(point_type[type_index]==1){
	return(blocks[index].region[dim-1].start);
    }
    else{
	return(blocks[index].region[dim-1].end);
    }
    
}

//////////////////////////////////////////////////////////////////////
// Function kdtree_mst
// the function constructs a global chain without considering the gap cost
// using  Abouelhoda-Ohlebusch algorithm.
//////////////////////////////////////////////////////////////////////

int kdChainer::kdtree_mst()
{
  Kdelem_ptr elem;
  double maxscore = 0;
  long permindex;
  long i;

  elem=&blocks[sorted_list[0]];
  elem->globalscore=elem->weight;
  elem->last=NULL;
  kdtreeobj.undelete(elem); //reinsert element into the tree
  for(i = 2; i<(numofblocks*2); i++){
      
      //elem->weight = GETREGIONLENGTH(elem->region[0]);
      elem=&blocks[sorted_list[i]];
      if(point_type[i]==1){
	  // get its right end
	  // add weight to its right end
	  permindex = kdtreeobj.nearestNeighbor(elem);
	  nndist=kdtreeobj.nndist;
	  
	  if(permindex < 0) //no such element
	  {
	      elem->last = NULL;
	      elem->globalscore = elem->weight;
	      if(elem->globalscore > maxscore) //new best chain was found
	      {
		  maxscore = elem->globalscore; 
		  tail = elem; //current last element of the chain
	      }
	  }
	  else{
	      elem->last = perm[permindex]; //setting chain reference
	      elem->globalscore = elem->weight - nndist;
	      if(elem->globalscore > maxscore) //new best chain was found
	      {
		  maxscore = elem->globalscore; 
		  tail = elem; //current last element of the chain
	      }
	  }
      }
      else{
	  // insert it in the tree
	  kdtreeobj.undelete(elem); //reinsert element into the tree
	  if(elem->globalscore > maxscore) //new best chain was found
	  {
	      maxscore = elem->globalscore; 
	      tail = elem; //current last element of the chain
	  }
	  
      }
  }
  global_maxscore=maxscore;
  return(1);
}

//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Function KdtreebasedModgetchain
// the function constructs a global chain without considering the gap cost
// by calling the function kd_mst which uses Abouelhoda-Ohlebusch algorithm.
// The function processes the input match file by 
// 1- calling the parsing functions
// 2- build the kdtree in k-1 dimensions. 
// 3- runs the kdtree_mst function
//  Activate if you want : measuring the time required for the chaining procedure
// the function enables the inclusion of an additional weight file
// the input IncludeDirection has no use in this version.
////////////////////////////////////////////////////////////////////////

int kdChainer::KdtreebasedModgetchain(ProcessFragmentFile* InfoObj, int IncludeWeight,int IncludeDirection) 
{
    // TODO: Add your command handler code here
    //////////timing in Unix//////////////
    float ClockTicksPerSecond; 
    struct tms StartTime; 
    struct tms EndTime; 
    float StartTimeSeconds; 
    float EndTimeSeconds; 
    ClockTicksPerSecond = (double)sysconf(_SC_CLK_TCK); 
    
    

    long i;
    //   char temp[MaxFileNameSize];
    
    
    dimensions=InfoObj->no_of_genomes-1;
    numofblocks=kd_read_fragment_file(InfoObj->filename);
    if(numofblocks<0){
	if((numofblocks==-1)||(numofblocks==-2)){
	    printf("Not enough memory\n");
	    return(-1);
	}
	else{
	    printf("Error reading input file \n");
	    return(-1);
	}
    }
    
    
    bucket_ptr = (Kdnode_ptr *) malloc(numofblocks * sizeof(Kdnode_ptr));
    if(bucket_ptr==NULL)
    {
	printf("Not Enough Memory\n");
	return(-1);
    }    
    
    perm = (Kdelem_ptr *) malloc(numofblocks * sizeof(Kdelem_ptr));
    if(perm==NULL)
    {
	printf("Not Enough Memory\n");
	return(-1);
    } 
    for(i=0;i<numofblocks;i++){
	perm[i] = &blocks[i];
    }
    //My kdTree code
    kdtreeobj.dimensions=dimensions;
    kdtreeobj.blocks=blocks;
    kdtreeobj.bucket_ptr=bucket_ptr;
    kdtreeobj.perm=perm;
    kdtreeobj.numofblocks=numofblocks;
       
    int res=quick_sort_list(InfoObj->no_of_genomes);  // with the original file dim        
    if(res<0){
	printf("Not enough memory \n");
	return(-1);
    }
    kdtreeobj.root = kdtreeobj.build(0, numofblocks-1, 0);
    if(kdtreeobj.root==NULL){
	printf("Not enough memory for building a kdtree \n");
	return(-1);
    }
    for(i=0;i<numofblocks;i++)
    {
	perm[i]->num = i; //setting num according to their order
    }
    
    kdtreeobj.DeleteAll();
    
   // timing

    times(&StartTime);   
    
    kdtree_mst();
    
    times(&EndTime); 
    StartTimeSeconds = StartTime.tms_utime/ClockTicksPerSecond; 
    EndTimeSeconds = EndTime.tms_utime/ClockTicksPerSecond; 
    //printf("time for chaining %f\n",EndTimeSeconds - StartTimeSeconds); 
    //printf("Score:%f\n",tail->globalscore);
    //char tempstr[300];
    //sprintf(tempstr,"Score:%f\n",tail->globalscore);
   
    if(chainer_format_flag)
	storeChain(InfoObj,IncludeDirection);
    else storeChainDB(InfoObj,0);
    
    //delete_all();
    free(bucket_ptr);bucket_ptr=NULL;
    
    for(i=0;i<numofblocks;i++)
    {
	free(blocks[i].region); //setting num according to their order
	blocks[i].region=NULL;
    }    
    free(blocks);blocks=NULL;
    free(perm);perm=NULL;    
    return(1);
}

/*======================Handling Gaps Functions ===================*/
//////////////////////////////////////////////////////////////////////
//*===================    Including Gap Costs    ===================*//
//////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
// Function GappedKdtreebasedModgetchain
// the function constructs a global chain considering the gap cost
// by calling the function kdtree_gapped_chain which uses 
// Abouelhoda-Ohlebusch algorithm.
// The function processes the input fragment file by 
// 1- calling the parsing functions
// 2- build the kdtree in k-1 dimensions. 
// 3- runs the kdtree_mst function
// Optional(activate if you want): measure the time required for the chaining procedure.
// The gap costs between two fragments is the L_1 distance between 
// their start and end points
///////////////////////////////////////////////////////////////////////
int kdChainer::GappedKdtreebasedModgetchain(ProcessFragmentFile* InfoObj,\
					    int IncludeWeight,int IncludeDirection) 
{
    // TODO: Add your command handler code here

  long i;
//  char temp[MaxFileNameSize];

//  float ClockTicksPerSecond; 
//  struct tms StartTime; 
//  struct tms EndTime; 
//  float StartTimeSeconds; 
//  float EndTimeSeconds; 
//  ClockTicksPerSecond = (double)sysconf(_SC_CLK_TCK); 

  
  dimensions=InfoObj->no_of_genomes-1;
  numofblocks=kd_read_fragment_file(InfoObj->filename);
  if(numofblocks<0){
      if((numofblocks==-1)||(numofblocks==-2)){
	  printf("Not enough memory\n");
	  return(-1);
      }
      else{
	  printf("Error reading input file \n");
	  return(-1);
      }
  }
  bucket_ptr = (Kdnode_ptr *) malloc(numofblocks * sizeof(Kdnode_ptr));
  if(bucket_ptr==NULL)
  {
      printf("Not Enough Memory \n");
      return(-1);
  }

  perm = (Kdelem_ptr *) malloc(numofblocks * sizeof(Kdelem_ptr));  
  if(perm==NULL)
  {
      printf("Not Enough Memory\n");
      return(-1);
  }
  gapped_kdelem(1);
  
  
  int res=quick_sort_list(InfoObj->no_of_genomes);  // with the original file dim
  if(res<0){
      printf("Not enough memory \n");
      return(-1);
  }
  //My kdTree code
  kdtreeobj.dimensions=dimensions;
  kdtreeobj.blocks=blocks;
  kdtreeobj.bucket_ptr=bucket_ptr;
  kdtreeobj.perm=perm;
  kdtreeobj.numofblocks=numofblocks;
  
  kdtreeobj.root = kdtreeobj.build(0, numofblocks-1, 0);
  if(kdtreeobj.root==NULL)
  {
      printf("Not enough memory for building a kdtree\n");
      return(-1);
  }

  for(i=0;i<numofblocks;i++)
  {
      perm[i]->num = i; //setting num according to their order
      perm[i]->start_child_list=NULL;
  }
  
  kdtreeobj.DeleteAll();
// timing
  
//   times(&StartTime);   
  this->kdtree_gapped_chain();  
//   times(&EndTime); 
//   StartTimeSeconds = StartTime.tms_utime/ClockTicksPerSecond; 
//   EndTimeSeconds = EndTime.tms_utime/ClockTicksPerSecond; 
//   printf("time for chaining %f\n",EndTimeSeconds - StartTimeSeconds); 
//   changes 2.06.2004
//  sprintf(temp,"Score: %f",tail->globalscore);
//  sprintf(temp,"Score: %f",tail->globalscore+(double)(0*(dimensions+1))+(100*(dimensions+1)));
//   printf("%s\n",temp);//   AfxMessageBox(temp);

  if(chainer_format_flag)
      storeChain(InfoObj,IncludeDirection);
  else 
      storeChainDB(InfoObj,1);

  //delete_all();
  free(bucket_ptr);bucket_ptr=NULL;
  for(i=0;i<numofblocks;i++)
  {
      free(blocks[i].region); //setting num according to their order
      blocks[i].region=NULL;
  }
  free(blocks);blocks=NULL;
  free(perm);perm=NULL;
  //kdChainer::delete_kdTree();
  return(1);
}

//////////////////////////////////////////////////////////////////////
// ----------Functionality for constructig local chains-------------//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function: KdtreeSW
// The function constructs local chains by basically calling 
// the function SW_chaining which is an implementation of 
// Abouelhoda-Ohlebusch algorithm.
// The function calls different functions to:
// 1) constructs k-1 file
// 2) Intialize the blocks and construct the kdtree
// 2) build the local chains
// 3) store the local chains
//////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Local chaining of fragments
int kdChainer::KdtreeSW(ProcessFragmentFile* InfoObj, int IncludeWeight,int IncludeDirection,double Threshold) 
{
    // TODO: Add your command handler code here

    long i;
//    char temp[MaxFileNameSize];
    int res;

//  float ClockTicksPerSecond; 
//  struct tms StartTime; 
//  struct tms EndTime; 
//  float StartTimeSeconds; 
//  float EndTimeSeconds; 
//  ClockTicksPerSecond = (double)sysconf(_SC_CLK_TCK); 


    //printf("haloooooooooooo noooooooooceeeeeeer\n");
  // the new format begin 
    dimensions=InfoObj->no_of_genomes-1;
    numofblocks=kd_read_fragment_file(InfoObj->filename);
    if(numofblocks<0){
	if((numofblocks==-1)||(numofblocks==-2)){
	    printf("Not enough memory\n");
	    return(-1);
	}
	else{
	    printf("Error reading input file \n");
	    return(-1);
	}
    }
    bucket_ptr = (Kdnode_ptr *) malloc(numofblocks * sizeof(Kdnode_ptr));
    if(bucket_ptr==NULL)
    {
	printf("Not Enough Memory");
	return(-1);
    }    
    perm = (Kdelem_ptr *) malloc(numofblocks * sizeof(Kdelem_ptr));
    if(perm==NULL)
    {
	printf("Not Enough Memory\n");
	return(-1);
    } 

    gapped_kdelem(1);	
    
    
    
    //// my test after initializing blocks
    //  quick_sort_list(4);  // with the original file dim
    res=quick_sort_list(InfoObj->no_of_genomes);  // with the original file dim
    if(res<0){
	printf("Not enough memory \n");
	return(-1);
    }
    
    //My kdTree code
    kdtreeobj.dimensions=dimensions;
    kdtreeobj.blocks=blocks;
    kdtreeobj.bucket_ptr=bucket_ptr;
    kdtreeobj.perm=perm;
    kdtreeobj.numofblocks=numofblocks;
    kdtreeobj.root = kdtreeobj.build(0, numofblocks-1, 0);
    if(kdtreeobj.root==NULL){
	printf("Not enough memory \n");
	return(-1);
    }
    for(i=0;i<numofblocks;i++)
    {
	perm[i]->num = i; //setting num according to their order
	perm[i]->start_child_list=NULL;
    }
    
    kdtreeobj.DeleteAll();
// timing
//   times(&StartTime);   
    
    if((detect_repeat_flag==0)){
	// gap threshold flag is to indicate of a gap constraint option is on
	res=-1;
	
	if(gap_threshold_flag==0){
	    res=this->SW_chaining(Threshold);
	}
	else{	
	    if(cluster_flag==1){
		Threshold=MIN_VALUE-(double)500;
		res=this->mod_SW_chaining(Threshold);
	    }
	    else{
		res=this->mod_SW_chaining(Threshold);	 
	    }
	}
    }    
    else{    	    
	mod_SW_repeat(Threshold);	
    }
    
//   times(&EndTime); 
//   StartTimeSeconds = StartTime.tms_utime/ClockTicksPerSecond; 
//   EndTimeSeconds = EndTime.tms_utime/ClockTicksPerSecond; 
//   printf("time for chaining %f\n",EndTimeSeconds - StartTimeSeconds); 
//   sprintf(temp,"Score: %f",tail->globalscore);
//   printf("%s\n",temp);
    // inverse the transform before storing the points
    if(overlapping_value>0){
	for(i=1;i<numofblocks;i++){
	    if(blocks[i].num!=0){	   
		inverse_transform_for_overlapping(&blocks[i]);
	    }
	}
    }

    if(res>0){    
	// get the global_max_score of all local chains
	silent_chain_traversal();
	if(chainer_format_flag){
	    storeMultiChain(InfoObj,IncludeDirection);
	}
	else
	{
	    storeMultiChainDB(InfoObj,IncludeDirection);
	}
    }
    //delete_all();
    free(bucket_ptr);bucket_ptr=NULL;
    free_blocks_array();
    free_chain_roots();
    free(perm);perm=NULL;
    if(res>0)return(1);
    else return(-1);
}

////////////////////////////////////////////////////////////////////77
//////////////////////////////////////////////////////////////////////
//Function: SW_chaining
// constructs local chains using Abouelhoda-Ohlebusch Algorithm.
// The connected fragments are forming clusters (subtrees) which
// are stored using linked lists.
/////////////////////////////////////////////////////////////////////

int kdChainer::SW_chaining(double Threshold)
{
    Kdelem_ptr elem;
//  Kdelem_ptr temp_elem;
    double tempscore;
    double maxscore = MIN_VALUE;
    long permindex;
    long i,j;
    long gap_ij;
    long gap_ti;
    child* tempchild;
    
    // Activate the ORIGIN point
    elem=&blocks[sorted_list[0]];
    // elem->globalscore=elem->weight; it is already initilaized
    
    elem->last=NULL;
    if(chainroots==NULL){
	chainroots=(child*)malloc(sizeof(child));
	if(chainroots==NULL){
	    return(-1);
	}
	chainroots->child_elem=elem;
	chainroots->next=NULL;
    }
    else{
	tempchild=(child*)malloc(sizeof(child));
	tempchild->next=NULL; // initialization
	tempchild->child_elem=elem;
	tempchild->next=chainroots;
	chainroots=tempchild;
	
    }
    
    kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
    for(i = 2; i<(numofblocks*2); i++)
    {
	
	//elem->weight = GETREGIONLENGTH(elem->region[0]);
	elem=&blocks[sorted_list[i]];
	if(point_type[i]==1){  // it is start point of a fragment
	    // get its right end
	    // add weight to its right end
	    permindex = kdtreeobj.gapped_nearestNeighbor(elem);
	    nndist=kdtreeobj.nndist;
	    //permindex = gapped_nearestNeighbor(elem);
	    if((permindex < 0)||(perm[permindex]->num==0)) //no such element
	    {		
		// start a new chain
		elem->last = NULL;
		// the elem as a new root of a new chain
		// and it is added to the root list
		if(chainroots==NULL){
		    chainroots=(child*)malloc(sizeof(child));
		    chainroots->child_elem=elem;
		    chainroots->next=NULL;
		}
		else{
		    tempchild=(child*)malloc(sizeof(child));		    
		    if(tempchild==NULL)
		    {
			return(-1);
		    }
		    tempchild->child_elem=elem;
		    tempchild->next=chainroots;
		    chainroots=tempchild;
		}
		
		if(elem->globalscore > maxscore) //new best chain was found
		{
		    maxscore = elem->globalscore; 
		    tail = elem; //current last element of the chain
		}
	    }
	    else{
		
		/// geting the gap between elem and perm[permindex]
		gap_ij=-1*(dimensions+1);
		gap_ti=0;
		for(j=0;j<dimensions;j++){
		    gap_ij=gap_ij+elem->region[j].start-perm[permindex]->region[j].end;
		    gap_ti=gap_ti+genome_sizes_array[j+1]-perm[permindex]->region[j].end;
		    
		}
		gap_ij=gap_ij+elem->region[dimensions].start-perm[permindex]->region[dimensions].end;
		gap_ti=gap_ti+genome_sizes_array[0]-perm[permindex]->region[dimensions].end;
		
		//gap_ij=1+(gap_ij/(dimensions+1)); 	
		tempscore= elem->weight + nndist+gap_ti-gap_ij;
		maxscore_gapti=tempscore;
		
		///============= Deciding if start new chain
		if((tempscore>Threshold)&&(gap_ij<gap_threshold)){
		    
		    tempscore= elem->globalscore + nndist+(double)gap_ti-(double)gap_ij;
		    elem->last = perm[permindex]; //setting chain reference	
		    elem->globalscore =tempscore;
		    // Adding the fragment to the child list
		    if(perm[permindex]->start_child_list==NULL){
			perm[permindex]->start_child_list=(child*)malloc(sizeof(child));
			if(perm[permindex]->start_child_list==NULL)
			{
			    return(-1);
			}
			perm[permindex]->start_child_list->next=NULL;
			perm[permindex]->start_child_list->child_elem=elem;
			//perm[permindex]->start_child_list=elem;
		    }
		    else{
			tempchild=(child*)malloc(sizeof(child));
			if(tempchild==NULL){
			    return(-1);
			}
			tempchild->child_elem=elem;
			tempchild->next=perm[permindex]->start_child_list;
			perm[permindex]->start_child_list=tempchild;
			
		    }
		    
		}
		else{
		    elem->last=NULL; 
		    // the elem as a new root of a new chain
		    // and it is added to the root list
		    if(chainroots==NULL){
			chainroots=(child*)malloc(sizeof(child));
			chainroots->child_elem=elem;
			chainroots->next=NULL;
		    }
		    else{
			tempchild=(child*)malloc(sizeof(child));
			if(tempchild==NULL){
			    return(-1);
			}
			tempchild->child_elem=elem;
			tempchild->next=chainroots;
			chainroots=tempchild;

		    }
		    
		}
		if(maxscore_gapti > maxscore) //new best chain was found
		{
		    maxscore = maxscore_gapti; 
		    tail = elem; //current last element of the chain
		    //maxscore_gapti=gap_ti;
		}								
	    }
	}
	else{
	    // activate it in the tree
	    kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
	    if(elem->globalscore > maxscore) //new best chain was found
	    {
		maxscore = elem->globalscore; 
		tail = elem; //current last element of the chain
	    }
	    
	}
    }
    global_maxscore=maxscore;
    return(1);
}

///////////////////////////////////////////////////////////////////////////////
// this function is modified by making RMQ over a bounded region
// to allow a certain threshold on the gap cost.

int kdChainer::mod_SW_chaining(double Threshold)
{
  Kdelem_ptr elem;
//  Kdelem_ptr temp_elem;
  double tempscore;
  double maxscore = MIN_VALUE;
  long permindex;
  long i,j;
  long gap_ij;
  long gap_ti;
  child* tempchild;
  long lower_limit=2;
  long delta_x;
  int flag=0;

  Kdelem_ptr lower_bound=new Kdelem;
  lower_bound->region=(Region*)malloc(sizeof(Region)*(dimensions+1));
  if(lower_bound->region==NULL){
      return(-1);
  }
  elem=&blocks[sorted_list[0]];
  // elem->globalscore=elem->weight; it is already initilaized
  elem->last=NULL;
  if(chainroots==NULL){
      chainroots=(child*)malloc(sizeof(child));
      if(chainroots==NULL)return(-1);
      chainroots->child_elem=elem;
      chainroots->next=NULL;
  }
  else{
      tempchild=(child*)malloc(sizeof(child));
      tempchild->next=NULL; // initialization
      tempchild->child_elem=elem;
      tempchild->next=chainroots;
      chainroots=tempchild;
      
  }
  
  //debug 
  //printf("no blocks %d\n",numofblocks);
  kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
  for(i = 2; i<(numofblocks*2); i++)
  {
      //elem->weight = GETREGIONLENGTH(elem->region[0]);
      elem=&blocks[sorted_list[i]];
      if(point_type[i]==1){  // it is start point of a fragment
	  	  	  
	  // get its right end
	  // add weight to its right end
	  
	  //permindex = kdtreeobj.gapped_nearestNeighbor(elem);
	  // here is the modfication w.r.t. RMQ
	  // We check if the lower limit is within the gap region or not.	  
	  //delta_x=blocks[sorted_list[lower_limit]].region[dimensions].end;
	  j=lower_limit;
	  while(j<i){
	      if(point_type[j]!=1){
		  delta_x=elem->region[dimensions].start-blocks[sorted_list[j]].region[dimensions].end;
		  if(delta_x>gap_threshold){
		      kdtreeobj.DeactivateElem(&blocks[sorted_list[j]]);
		      flag=1;
		  }
		  else{
		      break;
		  }
		  j++;
	      }
	      else{
		  j++;
	      }
	      
	      
	  }
	  lower_limit=j;

	  // assiging the reion boundaries
	  lower_bound->region[dimensions].start=blocks[sorted_list[lower_limit]].region[dimensions].start;
	  lower_bound->region[dimensions].end=lower_bound->region[dimensions].start;  		
	  for(j=0;j<dimensions;j++){	      
	      lower_bound->region[j].start=(long)max(0,elem->region[j].start-gap_threshold);
	      lower_bound->region[j].end=lower_bound->region[j].start;  		
	  }	
	  permindex=kdtreeobj.region_gapped_nearestNeighbor(elem,lower_bound);	
	  nndist=kdtreeobj.nndist;
	  

	  if((permindex < 0)||((flag==1)&&(perm[permindex]->num==0))||(perm[permindex]->num==0)) //no such element
	  {
	      // start a new chain
	      elem->last = NULL;
	      // the elem as a new root of a new chain
	      // and it is added to the root list
	      if(chainroots==NULL){
		  chainroots=(child*)malloc(sizeof(child));
		  if(chainroots==NULL)return(-1);
		  chainroots->child_elem=elem;
		  chainroots->next=NULL;
	      }
	      else{
		  tempchild=(child*)malloc(sizeof(child));
		  if(tempchild==NULL)return(-1);
		  tempchild->child_elem=elem;
		  tempchild->next=chainroots;
		  chainroots=tempchild;
	      }
	      
	      if(elem->globalscore > maxscore) //new best chain was found
	      {
		  maxscore = elem->globalscore; 
		  tail = elem; //current last element of the chain
	      }
	  }
	  else{
	      
	      /// geting the gap between elem and perm[permindex]
	      gap_ij=-1*(dimensions+1);;
	      gap_ti=0;
	      for(j=0;j<dimensions;j++){
		  gap_ij=gap_ij+elem->region[j].start-perm[permindex]->region[j].end;
		  gap_ti=gap_ti+genome_sizes_array[j+1]-perm[permindex]->region[j].end;
		  
	      }
	      gap_ij=gap_ij+elem->region[dimensions].start-perm[permindex]->region[dimensions].end;
	      gap_ti=gap_ti+genome_sizes_array[0]-perm[permindex]->region[dimensions].end;
	      
	      //gap_ij=1+(gap_ij/(dimensions+1)); 	
	      tempscore= elem->weight + nndist+gap_ti-gap_ij;
	      maxscore_gapti=tempscore;
	      
	      ///============= Deciding if start new chain
	      if(tempscore>Threshold){		  		  
		  tempscore= elem->globalscore + nndist+(double)gap_ti-(double)gap_ij;
		  elem->last = perm[permindex]; //setting chain reference	
		  elem->globalscore =tempscore;
		  // Adding the fragment to the child list
		  if(perm[permindex]->start_child_list==NULL){
		      perm[permindex]->start_child_list=(child*)malloc(sizeof(child));
		      if(perm[permindex]->start_child_list==NULL)return(-1);
		      perm[permindex]->start_child_list->next=NULL;
		      perm[permindex]->start_child_list->child_elem=elem;
		      //perm[permindex]->start_child_list=elem;
		  }
		  else{
		      tempchild=(child*)malloc(sizeof(child));
		      if(tempchild==NULL)return(-1);
		      tempchild->child_elem=elem;
		      tempchild->next=perm[permindex]->start_child_list;
		      perm[permindex]->start_child_list=tempchild;
		      
		  }
		  
	      }
	      else{
		  elem->last=NULL; 
		  // the elem as a new root of a new chain
		  // and it is added to the root list
		  if(chainroots==NULL){
		      chainroots=(child*)malloc(sizeof(child));
		      if(chainroots==NULL)return(-1);
		      chainroots->child_elem=elem;
		      chainroots->next=NULL;
		  }
		  else{
		      tempchild=(child*)malloc(sizeof(child));
		      if(tempchild==NULL)return(-1);		      
		      tempchild->child_elem=elem;
		      tempchild->next=chainroots;
		      chainroots=tempchild;
		      
		  }
		  
	      }
	      if(maxscore_gapti > maxscore) //new best chain was found
	      {
		  //maxscore = elem->globalscore; 
		  maxscore=maxscore_gapti;
		  tail = elem; //current last element of the chain
	      
	      }	      	        	      
	  }
      }
      else{
	  // insert it in the tree
	  kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
	  if(maxscore_gapti > maxscore) //new best chain was found
	  {
	      maxscore = maxscore_gapti; 
	      //maxscore= elem->globalscore; 
	      tail = elem; //current last element of the chain
	  }
	  
      }
  }
  free(lower_bound->region);
  delete lower_bound;
  global_maxscore=maxscore;
  return(1);
}
// end mod_SW_chaining
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// this function is modified by making RMQ over a bounded region
// to allow a certain threshold on the gap cost.

int kdChainer::mod_SW_repeat(double Threshold)
{
  Kdelem_ptr elem;
//  Kdelem_ptr temp_elem;
  double tempscore;
  double maxscore = MIN_VALUE;
  long permindex;
  long i,j;
  long gap_ij;
  long gap_ti;
  child* tempchild;
  long lower_limit=2;
  long delta_x;
  int flag=0;
  int Rewindflag=0;
  long* rootarray;
  long index_of_current_chain_root=0;

  Kdelem_ptr lower_bound=new Kdelem;
  lower_bound->region=(Region*)malloc(sizeof(Region)*(dimensions+1));
  if(lower_bound->region==NULL){
      return(-1);
  }
  elem=&blocks[sorted_list[0]];
  // elem->globalscore=elem->weight; it is already initilaized
  elem->last=NULL;
  if(chainroots==NULL){
      chainroots=(child*)malloc(sizeof(child));
      if(chainroots==NULL)return(-1);
      chainroots->child_elem=elem;
      chainroots->next=NULL;
  }
  else{
      tempchild=(child*)malloc(sizeof(child));
      tempchild->next=NULL; // initialization
      tempchild->child_elem=elem;
      tempchild->next=chainroots;
      chainroots=tempchild;
      
  }

  // Array to store for each element the root of its chain
  rootarray=(long*)malloc(sizeof(long)*numofblocks);
  for(long x=0;x<numofblocks;x++){
      rootarray[x]=x;
      
  }
  


  //debug 
  //printf("no blocks %d\n",numofblocks);
  kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
  for(i = 2; i<(numofblocks*2); i++)
  {
      //elem->weight = GETREGIONLENGTH(elem->region[0]);
      elem=&blocks[sorted_list[i]];
      if(point_type[i]==1){  // it is start point of a fragment
	  	  	  
	  // get its right end
	  // add weight to its right end
	  
	  //permindex = kdtreeobj.gapped_nearestNeighbor(elem);
	  // here is the modfication w.r.t. RMQ
	  // We check if the lower limit is within the gap region or not.	  
	  //delta_x=blocks[sorted_list[lower_limit]].region[dimensions].end;
	  j=lower_limit;
	  while(j<i){
	      if(point_type[j]!=1){
		  delta_x=elem->region[dimensions].start-blocks[sorted_list[j]].region[dimensions].end;
		  if(delta_x>gap_threshold){
		      kdtreeobj.DeactivateElem(&blocks[sorted_list[j]]);
		      flag=1;
		  }
		  else{
		      break;
		  }
		  j++;
	      }
	      else{
		  j++;
	      }
	      
	      
	  }
	  lower_limit=j;

	  
	  // assiging the reion boundaries
	  lower_bound->region[dimensions].start=blocks[sorted_list[lower_limit]].region[dimensions].start;
	  lower_bound->region[dimensions].end=lower_bound->region[dimensions].start;  		
	  for(j=0;j<dimensions;j++){	      
	      
	      lower_bound->region[j].start=(long)max(0,elem->region[j].start-gap_threshold);
	      lower_bound->region[j].end=lower_bound->region[j].start;  		
	  }	
	  permindex=kdtreeobj.region_gapped_nearestNeighbor(elem,lower_bound);	
	  nndist=kdtreeobj.nndist;
	  
          /// The distnict part for repeats
          /*
	    we check the condition that the end of the first component of the fragment
	    is before the start of the second component of the first fragment of the chain
	    
	   */
	  Rewindflag=0;
	  //bookmark
	 
	  if(permindex >= 0){
	      //printf("HALLO %u\n",(long)(perm[permindex]-blocks));
	      index_of_current_chain_root=rootarray[(long)(perm[permindex]-blocks)];
	      //index_of_current_chain_root=rootarray[];
	      if((blocks[index_of_current_chain_root].region[0].start)<=(elem->region[dimensions].end)){
		  Rewindflag=1;
		  //printf("HALLO %d %u\n",permindex,rootarray[permindex]);
	      }
	  }

	  if((permindex < 0)||((flag==1)&&(perm[permindex]->num==0))||(perm[permindex]->num==0)||(Rewindflag==1)) //no such element
	  {
	      // start a new chain
	      elem->last = NULL;
	      // the elem as a new root of a new chain
	      // and it is added to the root list
	      if(chainroots==NULL){
		  chainroots=(child*)malloc(sizeof(child));
		  if(chainroots==NULL)return(-1);
		  chainroots->child_elem=elem;
		  chainroots->next=NULL;
	      }
	      else{
		  tempchild=(child*)malloc(sizeof(child));
		  if(tempchild==NULL)return(-1);
		  tempchild->child_elem=elem;
		  tempchild->next=chainroots;
		  chainroots=tempchild;
	      }
	      
	      if(elem->globalscore > maxscore) //new best chain was found
	      {
		  maxscore = elem->globalscore; 
		  tail = elem; //current last element of the chain
	      }
	      //rootarray[sorted_list[i]]=sorted_list[i];
	  }
	  else{
	      
	      /// geting the gap between elem and perm[permindex]
	      gap_ij=-1*(dimensions+1);;
	      gap_ti=0;
	      for(j=0;j<dimensions;j++){
		  gap_ij=gap_ij+elem->region[j].start-perm[permindex]->region[j].end;
		  gap_ti=gap_ti+genome_sizes_array[j+1]-perm[permindex]->region[j].end;
		  
	      }
	      gap_ij=gap_ij+elem->region[dimensions].start-perm[permindex]->region[dimensions].end;
	      gap_ti=gap_ti+genome_sizes_array[0]-perm[permindex]->region[dimensions].end;
	      
	      //gap_ij=1+(gap_ij/(dimensions+1)); 	
	      tempscore= elem->weight + nndist+gap_ti-gap_ij;
	      maxscore_gapti=tempscore;
	      
	      ///============= Deciding if start new chain
	      if(tempscore>Threshold){		  		  
		  tempscore= elem->globalscore + nndist+(double)gap_ti-(double)gap_ij;
		  elem->last = perm[permindex]; //setting chain reference	
		  elem->globalscore =tempscore;
		  // Adding the fragment to the child list
		  if(perm[permindex]->start_child_list==NULL){
		      perm[permindex]->start_child_list=(child*)malloc(sizeof(child));
		      if(perm[permindex]->start_child_list==NULL)return(-1);
		      perm[permindex]->start_child_list->next=NULL;
		      perm[permindex]->start_child_list->child_elem=elem;
		      //perm[permindex]->start_child_list=elem;
		  }
		  else{
		      tempchild=(child*)malloc(sizeof(child));
		      if(tempchild==NULL)return(-1);
		      tempchild->child_elem=elem;
		      tempchild->next=perm[permindex]->start_child_list;
		      perm[permindex]->start_child_list=tempchild;
		      
		  }

		  // update the start point of the chain		  
		  rootarray[(long)(elem-blocks)]=index_of_current_chain_root;		 
		 
		  
	      }
	      else{
		  elem->last=NULL; 
		  // the elem as a new root of a new chain
		  // and it is added to the root list
		  if(chainroots==NULL){
		      chainroots=(child*)malloc(sizeof(child));
		      if(chainroots==NULL)return(-1);
		      chainroots->child_elem=elem;
		      chainroots->next=NULL;
		  }
		  else{
		      tempchild=(child*)malloc(sizeof(child));
		      if(tempchild==NULL)return(-1);		      
		      tempchild->child_elem=elem;
		      tempchild->next=chainroots;
		      chainroots=tempchild;
		      
		  }
		  
	      }
	      if(maxscore_gapti > maxscore) //new best chain was found
	      {
		  //maxscore = elem->globalscore; 
		  maxscore=maxscore_gapti;
		  tail = elem; //current last element of the chain
	      
	      }	      	        	      
	  }
      }
      else{
	  // insert it in the tree
	  kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
	  if(maxscore_gapti > maxscore) //new best chain was found
	  {
	      maxscore = maxscore_gapti; 
	      //maxscore= elem->globalscore; 
	      tail = elem; //current last element of the chain
	  }
	  
      }
  }
  free(lower_bound->region);
  delete lower_bound;
  global_maxscore=maxscore;
  free(rootarray);
  return(1);
}
// end mod_SW_repeat
//
//////////////////////////////////////////////////////////////////////
//Function:storeMultiChain
// stores the clusters obtained from the SW_chaining in a file"matchfile.cst".
// where matchfile is the name of the match file.
// By a traversal of every cluster (subtree) a representative chain 
// (a chain with the highest score) is reported and stored in the file
// "fragmentfile.chn" and compact chain file "framentfile.ccn".
//////////////////////////////////////////////////////////////////////

void kdChainer::storeMultiChain(ProcessFragmentFile* InfoObj,int IncludeDirection)
{
    Kdelem_ptr elem = NULL;  //last element of the chain
    long i;
    FILE* fptr;
    FILE* clusterfptr;
//    FILE* fptr2;
    FILE* compactchainfptr;
    child* ptr;
    long j;
    char temp[MaxFileNameSize];
//    char tempstr[MaxFileNameSize];
    long* regionend=(long*)malloc(sizeof(long)*(dimensions+1));
    int flag=0;
    long chain_counter=0;
 
    sprintf(temp,"%s.chn",InfoObj->filename);
    //Kdelem_ptr local_tail=NULL;
    double maxscore = MIN_VALUE-10;
    
    if(chainroots==NULL)return;
    fptr=FOPEN(temp,"w");
    
    // opening cluster file
    if(cluster_file_flag)
    sprintf(temp,"%s.cst",InfoObj->filename);
    
    clusterfptr=FOPEN(temp,"w");
    
    // opening compact chain file
    sprintf(temp,"%s.ccn",InfoObj->filename);   
    compactchainfptr=FOPEN(temp,"w");

    fprintf(fptr,">CHA %d\n",(int)dimensions+1);
    fprintf(compactchainfptr,">CHA %d\n",(int)dimensions+1);
    
    for(ptr=chainroots;ptr!=NULL ; ptr = ptr->next)  //last references the previous elem
    {
	// putting the clusters in file
	// CHANGE: Allow Storing chains of one element
	//if(ptr->child_elem->start_child_list==NULL) continue;
	
	tail=NULL;
	maxscore = MIN_VALUE-10;
	if(cluster_file_flag){
	    if(ptr->child_elem->num!=0)fprintf(clusterfptr,"#cluster\n");
	}
	// the following function store clusters of chains and returns a pointer
        // to the opimal chain in this cluster
	
	if(cluster_file_flag) traverse_store_chain(ptr,clusterfptr,&maxscore);
	else silent_traverse_cluster(ptr,&maxscore);

	// putting the optimal chain of cluster in other file
	long gap_ti=0;
	for(j=0;j<dimensions;j++){
	    gap_ti=gap_ti+genome_sizes_array[j+1]-tail->region[j].end;
	    
	}
	gap_ti=gap_ti+genome_sizes_array[0]-tail->region[dimensions].end;
	
	double tmpcheckchain;
	// num =0 means that this fragment is the ORIGIN
	if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
	    tmpcheckchain=tail->globalscore+gap_ti;
	    
	}
	else{	    	    
	    if (tail->last!=NULL) // the chain is not ending with the origin
		tmpcheckchain=tail->globalscore+gap_ti+(shiftvalue*(dimensions+1));
	    else tmpcheckchain=tail->globalscore+gap_ti;	    
	}
	
	// check chain against the options: score and chain depth
	if(cluster_flag==0)
	    chain_score=max(chain_score,global_maxscore*(percentagescore/100)); 
	if(check_chain(ptr->child_elem,chain_depth,chain_score,tmpcheckchain)!=1)continue;       
	if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
	    fprintf(fptr, "#%lf\n",(tail->globalscore+gap_ti));
	    fprintf(compactchainfptr, "#%lf\n",tail->globalscore+gap_ti);
	    flag=0;
	}
	else{
	    if (tail->last!=NULL){ // tail->last is the ORIGIN, the shiftvalue should be considered
		fprintf(fptr, "#%lf\n",(tail->globalscore+gap_ti+(shiftvalue*(dimensions+1))));
		fprintf(compactchainfptr, "#%lf\n",tail->globalscore+gap_ti+(shiftvalue*(dimensions+1)));
	    }
	    else{
		//fprintf(fptr, "#%lf\n",tail->globalscore+gap_ti);
		//fprintf(compactchainfptr, "#%lf\n",tail->globalscore+gap_ti);
	    }
	    flag=0;
	}
	// store the fragments of the chain
	for(elem=tail;elem;elem=elem->last){
	    
	    if (elem->num!=0){
		fprintf(fptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,
			elem->region[dimensions].end-shiftvalue);
	    }
	    	    
	    for(i=0;i<dimensions;i++){
		if (elem->num!=0){
		    fprintf(fptr,"[%lu,%lu] ", elem->region[i].start-shiftvalue,elem->region[i].end-shiftvalue);
		}	       
	    }
	    if (elem->num!=0) {
		fprintf(fptr,"%lu",(long)(elem-blocks));
		fprintf(fptr,"\n");
	    }
	    //// storing compact chain
	    if(flag==0){
		regionend[0]=elem->region[dimensions].end;
		flag=1;
		for(j=0;j<dimensions;j++){
		    regionend[j+1]=elem->region[j].end;
		}	
		// Chain of 1 fragment
		
		if(elem->last==NULL){
		    if (elem->num!=0){
			fprintf(compactchainfptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,\
			  elem->region[dimensions].end-shiftvalue);
			chain_counter++;
		    }      
		    
		    for(i=0;i<dimensions;i++){
			if (elem->num!=0){
			    fprintf(compactchainfptr,"[%lu,%lu] ", elem->region[i].start-shiftvalue,
				    elem->region[i].end-shiftvalue);
			}		  
		    }
		    if (elem->num!=0){
			//fprintf(compactchainfptr,"%lu",(long)(elem-blocks));
			fprintf(compactchainfptr,"\n"); 
		    }
		}
	    }
	    else if(elem->last==NULL){
		if (elem->num!=0){
		    fprintf(compactchainfptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,
			    regionend[0]-shiftvalue);		  
		    chain_counter++;
		}
		
		for(j=0;j<dimensions;j++)
		{
		    if (elem->num!=0){
			fprintf(compactchainfptr,"[%lu,%lu] ", 
				elem->region[j].start-shiftvalue,regionend[j+1]-shiftvalue);		  
		    }
		}	
		if (elem->num!=0){
		    //   fprintf(compactchainfptr,"%lu",(long)(elem-blocks));
		    fprintf(compactchainfptr,"\n");
		}
		
	    }	   	    	    
	}
	
	
	
    }
    fclose(fptr);
    if(cluster_file_flag)fclose(clusterfptr);
    fclose(compactchainfptr);
// Append information to the statistics file .stc
      sprintf(temp,"%s.stc",InfoObj->filename);     
      FILE* stcptr=FOPEN(temp,"a");
      if(stcptr==NULL){
	  printf("Cannot append to statistics file\n");
	  return;
      }      
      printf("# no. of chains: %lu \n",chain_counter);
      fprintf(stcptr,"# no. of chains: %lu \n",chain_counter);
      fclose(stcptr);
//
    free(regionend);
}
///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//Function:traverse_store_chain
// traverses the subtree of every cluster to 
// 1) store the cluster itself.
// 2) the leaf (tail of chain) with the highest
// score in this cluster is stored in the tail.
//////////////////////////////////////////////////////////////////////

Kdelem_ptr kdChainer::traverse_store_chain(child* root, FILE* fptr, double* maxscoreptr){
	
	Kdelem_ptr elem;	
	child* child_ptr;
	int i;	
	if(root==NULL)return(NULL);
	if(root->child_elem==NULL)return(NULL);
	elem=root->child_elem;
	//tail=elem;
	// here I can process the element
	if(elem->num!=0){
	    
	    fprintf(fptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,
		    elem->region[dimensions].end-shiftvalue);
	}
	if(elem->num!=0){
	    for(i=0;i<dimensions;i++)
	    {	    
		fprintf(fptr,"[%lu,%lu] ", elem->region[i].start-shiftvalue,elem->region[i].end-shiftvalue);
	    }
	}
	if(elem->num!=0){
	    fprintf(fptr,"%lu",(long)(elem-blocks)); //mohamed
	    fprintf(fptr,"\n");	
	}
	if(elem->globalscore > *maxscoreptr) //new best chain was found
	{
	    *maxscoreptr = elem->globalscore; 
	    tail = elem; //current last element of the chain
	}
	child_ptr=elem->start_child_list;
	while(child_ptr!=NULL){
	    traverse_store_chain(child_ptr,fptr,maxscoreptr);
	    child_ptr=child_ptr->next;
	    
	}
	return(NULL);
}

//////////////////////////////////////////////////////////////////////
//------------ Functionality related to filtering out chains--------//
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//Function: check chain
// checks the chain depth (number of fragments on a chain)
// checks if the given chain has a score larger than the given score 
//////////////////////////////////////////////////////////////////////

int kdChainer::check_chain(Kdelem_ptr this_chain_root, int depth, double th_score,double given_score){
    
    double average_length=0;
    Kdelem_ptr depth_ptr=tail;
    // check score
    if(given_score<th_score)return(-1);
    // check average length
    for(int i=0;i<=dimensions;i++){
	average_length=average_length+tail->region[i].end-this_chain_root->region[i].start;
    }

    //average_length=average_length+tail->region[dimensions].end-this_chain_root->region[dimensions].end;
    average_length=average_length/(double)(dimensions+1);
    if(average_length<chain_average_length)return(-1);
    for(int i=0;i<depth;i++){
	if(depth_ptr->last==NULL)return(-1);
	depth_ptr=depth_ptr->last;
    }
    
    
    return(1);
}

//////////////////////////////////////////////////////////////////////
//------------ Functionality related to freeing Memory--------//
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Function: free_blocks_array
// The function free the block array
//////////////////////////////////////////////////////////////////////
void kdChainer::free_blocks_array(){
	
//	Kdelem_ptr elem;
    child* child_ptr;
    long j;
    if(blocks==NULL)return;
    for(j=0;j<numofblocks;j++){
	child_ptr=blocks[j].start_child_list;
	while(child_ptr!=NULL){
	    blocks[j].start_child_list=blocks[j].start_child_list->next;
	    free(child_ptr);
	    child_ptr=blocks[j].start_child_list;
	    
	}
	free(blocks[j].region); blocks[j].region=NULL;
	
    }
    free(blocks);blocks=NULL;
}


//////////////////////////////////////////////////////////////////////
//Function: free_chain_roots
// frees the chain roots created in SW_chaining
/////////////////////////////////////////////////////////////////////

void kdChainer::free_chain_roots(){
	
//	Kdelem_ptr elem;
	child* child_ptr;
//	long j;
 
	child_ptr=chainroots;
	while(child_ptr!=NULL){
		chainroots=chainroots->next;
		free(child_ptr);
		child_ptr=chainroots;

	} 
	chainroots=NULL;
 
}

//////////////////////////////////////////////////////////////////////
// Function: storeChain
// it stores the global chains either generated while including the
// gap costs or not in a file named "matchfile.chn", where matchfile 
// is the name of the match file
//////////////////////////////////////////////////////////////////////
void kdChainer::storeChain(ProcessFragmentFile* InfoObj,int IncludeDirection)
{
  Kdelem_ptr elem = tail;  //last element of the chain
  long i;
  FILE* fptr;
  char temp[MaxFileNameSize];

//  char tempstr[MaxFileNameSize];
 

  
  sprintf(temp,"%s.chn",InfoObj->filename);
  fptr=FOPEN(temp,"w");
  fprintf(fptr,">CHA %d\n",(int)dimensions+1);
  // adding the score of the terminus  
  //tail->globalscore=tail->globalscore+multiplyingfactor;

  fprintf(fptr,"#%lf\n",(float)(tail->globalscore+multiplyingfactor));
  // storing the terminus
  for(i=0;i<=dimensions;i++)
  {
      fprintf(fptr,"[%lu,%lu] ",genome_sizes_array[i],genome_sizes_array[i]);
  }
  fprintf(fptr,"%lu",numofblocks); //mohamed
  fprintf(fptr,"\n");
  // storing the rest of the chain
  for(; elem; elem = elem->last){  //last references the previous elem
      if(elem->num!=0){	   
	  if(overlapping_value>0)inverse_transform_for_overlapping(elem);
      }
      if (elem->last!=NULL){
	  fprintf(fptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,\
		  elem->region[dimensions].end-shiftvalue);
      }
      else{
	  fprintf(fptr,"[%lu,%lu] ", elem->region[dimensions].start,elem->region[dimensions].end);
      }
      for(i=0;i<dimensions;i++)
      {
	  if (elem->last!=NULL){
	      fprintf(fptr,"[%lu,%lu] ", elem->region[i].start-shiftvalue,\
		      elem->region[i].end-shiftvalue);
	  }
	  else{
	      fprintf(fptr,"[%lu,%lu] ", elem->region[i].start,elem->region[i].end);
	  }
      }
      //fprintf(fptr," %lf",(float)elem->globalscore);
      fprintf(fptr,"%lu",(long)(elem-blocks)); //mohamed
      fprintf(fptr,"\n");
    }
  fclose(fptr);
// Append information to the statistics file .stc
      sprintf(temp,"%s.stc",InfoObj->filename);     
      FILE* stcptr=FOPEN(temp,"a");
      if(stcptr==NULL){
	  printf("Cannot append to statistics file\n");
	  return;
      }      
      printf("# no. of chains: 1 \n");
      fprintf(stcptr,"# no. of chains: 1 \n");
      fclose(stcptr);
//
  
}


///////////////////////////////////////////////////////////////////////////
// 


int kdChainer::kdtree_gapped_chain()
{
  Kdelem_ptr elem;
//  Kdelem temp_elem;
  double maxscore = MIN_VALUE;
  long permindex;
  long i,j;
  long gap_ij;
  long gap_ti;

  elem=&blocks[sorted_list[0]];
  // elem->globalscore=elem->weight; it is already initilaized
  elem->last=NULL;
  kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
  //printf("%lf",(float)elem->globalscore);
  for(i = 2; i<(numofblocks*2); i++)
  {
      //elem->weight = GETREGIONLENGTH(elem->region[0]);
      elem=&blocks[sorted_list[i]];
      if(point_type[i]==1){  // it is start point of a fragment
	  // get its right end
	  // add weight to its right end
	  permindex = kdtreeobj.gapped_nearestNeighbor(elem);
	  nndist=kdtreeobj.nndist;
	  //permindex = nearestNeighbor(elem);
	  if(permindex < 0) //no such element
	  {
	      elem->last = NULL;
	      //elem->globalscore = elem->weight;
	      if(elem->globalscore > maxscore) //new best chain was found
	      {
		  maxscore = elem->globalscore; 
		  tail = elem; //current last element of the chain
	      }
	  }
	  else{
	      elem->last = perm[permindex]; //setting chain reference
	      //temp_elem= *perm[permindex];
	      /// geting the gap between elem and perm[permindex]
	      gap_ij=0;
	      gap_ti=0;
	      
	      for(j=0;j<dimensions;j++){
		  //gap_ij=gap_ij+elem->region[j].start-perm[permindex]->region[j].end;
		  gap_ij=gap_ij+elem->region[j].start-elem->last->region[j].end;
		  //gap_ti=gap_ti+genome_sizes_array[j+1]-perm[permindex]->region[j].end;
		  gap_ti=gap_ti+genome_sizes_array[j+1]-elem->last->region[j].end;
		  
	      }
	      //gap_ij=gap_ij+elem->region[dimensions].start-perm[permindex]->region[dimensions].end;
	      gap_ij=gap_ij+elem->region[dimensions].start-elem->last->region[dimensions].end;
	      //gap_ti=gap_ti+genome_sizes_array[0]-perm[permindex]->region[dimensions].end;
	      gap_ti=gap_ti+genome_sizes_array[0]-elem->last->region[dimensions].end;
	      
	      elem->globalscore = elem->globalscore + nndist+(double)gap_ti-(double)gap_ij;
	      
	      if(elem->globalscore > maxscore) //new best chain was found
	      {
		  maxscore = elem->globalscore; 
		  tail = elem; //current last element of the chain
	      }
	  }
      }
      else{
	  // insert it in the tree
	  kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
	  if(elem->globalscore > maxscore) //new best chain was found
	  {
	      maxscore = elem->globalscore; 
	      tail = elem; //current last element of the chain
	  }
	  
      }
  }
  global_maxscore=maxscore;
  return(1);

}

///////////////////////////////////////////////////////////////////////////
// ///////////////////      EST and Contigs Chaining     //////////////////////////
///////////////////////////////////////////////////////////////////////////
/*
If flag is 0,
then the function performs chaining of fragments of a finished genome and 
a draft genome composed of a set of contigs. 
If flag is 1,
then the function performs chaining of fragments of a genomic sequence and
a database (set) of cDNA/EST sequences.
The input IncludeWeight and IncludeDirection are obsolete
*/


int kdChainer::ContigESTChainer(ProcessFragmentFile* InfoObj, int IncludeWeight,int IncludeDirection,double Threshold, char* contigfilename, int flag) 
{
	// TODO: Add your command handler code here
	
  long i;
//  char temp[MaxFileNameSize];
  int res=-1;
//  float ClockTicksPerSecond; 
//  struct tms StartTime; 
//  struct tms EndTime; 
//  float StartTimeSeconds; 
//  float EndTimeSeconds; 
//  ClockTicksPerSecond = (double)sysconf(_SC_CLK_TCK); 

  

  dimensions=InfoObj->no_of_genomes-1;
  numofblocks=kd_read_fragment_file(InfoObj->filename);
  if(numofblocks<0){
      if((numofblocks==-1)||(numofblocks==-2)){
	  printf("Not enough memory\n");
	  return(-1);
      }
      else{
	  printf("Error reading input file \n");
	  return(-1);
      }
  }

  bucket_ptr = (Kdnode_ptr *) malloc(numofblocks * sizeof(Kdnode_ptr));
  if(bucket_ptr==NULL)
  {
      printf("Not Enough Memory \n ");
      return(-1);
  }
  
  perm = (Kdelem_ptr *) malloc(numofblocks * sizeof(Kdelem_ptr));
  
  if(perm==NULL)
  {
      printf("Not Enough Memory\n");
      return(-1);
  } 
    
  gapped_kdelem(1);
    
  res=quick_sort_list(InfoObj->no_of_genomes);  // with the original file dim
  if(res<0){
      return(-1);
  }
  
  //My kdTree code
  kdtreeobj.dimensions=dimensions;
  kdtreeobj.blocks=blocks;
  kdtreeobj.bucket_ptr=bucket_ptr;
  kdtreeobj.perm=perm;
  kdtreeobj.numofblocks=numofblocks;
  
  kdtreeobj.root = kdtreeobj.build(0, numofblocks-1, 0);
  for(i=0;i<numofblocks;i++)
  {
      perm[i]->num = i; //setting num according to their order
      perm[i]->start_child_list=NULL;
  }
  
  kdtreeobj.DeleteAll();
// timing
//   times(&StartTime);   
// timing
// CTime startTime = CTime::GetCurrentTime();
// DWORD starttick=::GetTickCount();
  contigobj=(ContigSorting*)malloc(sizeof(ContigSorting));
  if(contigobj==NULL)return(-1);
  //printf("GAP THRESHOLD %lu\n",gap_threshold);

  if(detect_repeat_cluster_flag==1){      
      repeat_cluster(Threshold);
      //printf("haloooooooooooo clusteeeeeeeeerd\n");
  }
  else{
      if(flag==1){
	  res=this->mod_EST_SW(Threshold,contigfilename);
	  
      }
      else{
	  res=this->mod_Contig_SW(Threshold,contigfilename);
	  
      }
  }
  
  
  //	DWORD endtick=::GetTickCount();
  //CTime endTime = CTime::GetCurrentTime();
  //CTimeSpan elapsedTime = endTime - startTime;
  //sprintf(temp,"Score %lf time %2f sec & %4f msec",(float)tail->globalscore,(float)elapsedTime.GetTotalSeconds(),(float)(endtick-starttick));
  //AfxMessageBox(temp);
  //printf("%s",temp);

  // inverse the transform before storing the points
  if(overlapping_value>0){
      for(i=1;i<numofblocks;i++){
	  // 
	  //
	  if(blocks[i].num!=0){	   
	      inverse_transform_for_overlapping(&blocks[i]);
	  }
      }
  }
  if(res>0){
      if((flag==1)||(detect_repeat_cluster_flag==1)){
	  
	  if(chainer_format_flag==1)
	  {	      
	      storeMultiEST(InfoObj,IncludeDirection);  
	  }
	  else{
	      
	      storeMultiEstDB(InfoObj,IncludeDirection); 
	  }
	  
      }
      else{
	  silent_chain_traversal();	  
	  if(chainer_format_flag==1){
	      storeMultiChain(InfoObj,IncludeDirection);
	  }
	  else{
	      storeMultiContigDB(InfoObj,IncludeDirection);
	  }
      }
  }
  
  free(contigobj->recordseps);contigobj->recordseps=NULL;
  free(contigobj);        
  contigobj=NULL;
 
  
  
  //	storeMultiChain(InfoObj,IncludeDirection);
  
//   times(&EndTime); 
//   StartTimeSeconds = StartTime.tms_utime/ClockTicksPerSecond; 
//   EndTimeSeconds = EndTime.tms_utime/ClockTicksPerSecond; 
//   printf("time for chaining %f\n",EndTimeSeconds - StartTimeSeconds); 
//   sprintf(temp,"Score: %f",tail->globalscore);
//   printf("%s\n",temp);

  
//  freeing memory
  if(res>0){
      free(bucket_ptr);bucket_ptr=NULL;
      free_blocks_array();
      free_chain_roots();
      free(perm);  perm=NULL;
      return(1);
  }
  else{
      return(0);
  }
}

//////////////////////////////////////////////////////////////////////
//////////////   Aligining Contigs to a finished Genome    /////////////
//////////////////////////////////////////////////////////////////////
// Function mod_Contig_SW: 
// It chains the fragments in the contigs taking the gap cost
// into account.In other words, a fragment is connected to the
// tail of a chain if it is in the same contig and the score > 0
////////////////////////////////////////////////////////////////////////////

int kdChainer::mod_Contig_SW(double Threshold,char* contigfilename)
{

    //ContigSorting contigobj;
    Kdelem_ptr elem;
//  Kdelem_ptr temp_elem;
    double tempscore;
    double maxscore = MIN_VALUE;
    long permindex;
    long i,j;
    long gap_ij;
    long gap_ti;
    child* tempchild;
    
    ///////Contigs-specific Functionality
    long contig1=0;
    long contig2=0;
    long current_contig=0;
    //char firsttimeflag=0;
    long lower_limit;
    int res;
    
    Kdelem_ptr lower_bound=new Kdelem;
    lower_bound->region=(Region*)malloc(sizeof(Region)*(dimensions+1));
    if(lower_bound->region==NULL)return(-1);

    res=contigobj->read_kd_contigs(contigfilename,0);

    if(res<0)return(-1);
    
    
    elem=&blocks[sorted_list[0]];
    // elem->globalscore=elem->weight; it is already initilaized
    elem->last=NULL;
    if(chainroots==NULL){
	chainroots=(child*)malloc(sizeof(child));
	if(chainroots==NULL)return(-1);
	chainroots->child_elem=elem;
	chainroots->next=NULL;
    }
    else{
	tempchild=(child*)malloc(sizeof(child));
	if(tempchild==NULL)return(-1);
	tempchild->next=NULL; // initialization
	tempchild->child_elem=elem;
	tempchild->next=chainroots;
	chainroots=tempchild;
	
    }

    long * contig_id_array=new long [dimensions+1];
    long * positions=new long [dimensions+1];
    
    kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
    lower_limit=2;
    for(i = 2; i<(numofblocks*2); i++)
    {	
	elem=&blocks[sorted_list[i]];
	
	if(point_type[i]==1){  // it is start point of a fragment
	    
	    // for gapthreshold
	    
	    j=lower_limit;
	    while(j<i){
		if(point_type[j]!=1){
		    long delta_x=elem->region[dimensions].start-blocks[sorted_list[j]].region[dimensions].end;
		    if(delta_x>gap_threshold){
			kdtreeobj.DeactivateElem(&blocks[sorted_list[j]]);
			//flag=1;
		    }
		    else{
			break;
		    }
		    j++;
		}
		else{
		    j++;
		}
		
		
	    }
	    lower_limit=j;
	    
	    //end for gap threshold

	    // get its right end
	    // add weight to its right end
	    /// We also have to check if the two fragments are in the same conitg
	    
//	    contig2=contigobj->getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
//					   contigobj->total_conitg_length,elem->region[dimensions].end-shiftvalue-1);
//			contig1=contigobj.getrecordnum(contigobj.recordseps,contigobj.number_of_contigs,
//				contigobj.total_conitg_length,perm[permindex]->region[dimensions].start);

	    //printf("HALOOOOOOOOOOOOOOOO %u\n",contigobj->offset_array[0]);
/*
	    contig2=contigobj->getrecordnum(contigobj->recordseps+contigobj->offset_array[0],contigobj->number_of_contigs_array[0],
					   contigobj->total_contig_length_array[0],elem->region[dimensions].start-shiftvalue-1);
	    
*/

	    positions[0]=elem->region[dimensions].start-shiftvalue+1;
	    long kk;
	    for(kk=0;kk<dimensions;kk++){
		contig_id_array[kk]=0;
		positions[kk+1]=elem->region[kk].start-shiftvalue+1;
	    }
	    contig_id_array[kk]=0;
	    contig2=contigobj->getrecordnumGeneralized(positions,contig_id_array);

	    

	    if(current_contig==contig2) // still in the same EST/contig
	    {
		// a modification using the region search function
		// get the lower_bound point	
						
		//lower_bound->region[dimensions].start=(long)max(0,elem->region[dimensions].start-gap_threshold);
		//lower_bound->region[dimensions].start=(long)max((contigobj->recordseps+contigobj->offset_array[0])[contig_id_array[0]],
		//						elem->region[dimensions].start-gap_threshold);
		//lower_bound->region[dimensions].end=lower_bound->region[dimensions].start; 
		if(contig_id_array[0]!=0){
		    lower_bound->region[dimensions].start=(long)max((contigobj->recordseps+contigobj->offset_array[0])[contig_id_array[0]-1]+1,
							   elem->region[dimensions].start-gap_threshold);		    
		}
		else{
		    lower_bound->region[dimensions].start=(long)max(0,elem->region[dimensions].start-gap_threshold);
		}
		lower_bound->region[dimensions].start=lower_bound->region[dimensions].start+shiftvalue;
		lower_bound->region[dimensions].end=lower_bound->region[dimensions].start;
		
                // mohamed debug
		//if(elem->region[dimensions].start==875784){
		//  printf("HALLOOOOOOO ElemStart %u,  lowerbound %u\n",elem->region[dimensions].start,lower_bound->region[dimensions].start);
		//}

		for(j=0;j<dimensions;j++){		    
		    //lower_bound->region[j].start=(long)max(0,elem->region[j].start-gap_threshold);
		    if(contig_id_array[j+1]!=0){
			lower_bound->region[j].start=(long)max((contigobj->recordseps+contigobj->offset_array[j+1])[contig_id_array[j+1]-1]+1,
							   elem->region[j].start-gap_threshold);
			//printf("Dim %u, ElemStart %u, GC %u, lowerbound %u\n",j+1,elem->region[j].start,gap_threshold,
			//     (contigobj->recordseps+contigobj->offset_array[j+1])[contig_id_array[j+1]-1]);
		    }
		    else{
			lower_bound->region[j].start=(long)max(0,elem->region[j].start-gap_threshold);
		    }
		    //if(elem->region[dimensions].start==875784){
		    //	printf("HALLOOOOOOOOO ElemStart %u,  lowerbound %u\n",elem->region[j].start,lower_bound->region[j].start);
		    //}
		    lower_bound->region[j].start=lower_bound->region[j].start+shiftvalue;
		    lower_bound->region[j].end=lower_bound->region[j].start;  				    
		}	
				
		permindex=kdtreeobj.region_gapped_nearestNeighbor(elem,lower_bound);
		
		//if(elem->region[dimensions].start==875784)
		//printf("Perm index %u, pos: st %u, end %u,---lbs %u, lbe %u\n",permindex, 
		//perm[permindex]->region[0].start,perm[permindex]->region[0].end,lower_bound->region[0].start,lower_bound->region[0].end);
		
		// permindex = kdtreeobj.gapped_nearestNeighbor(elem);
		nndist=kdtreeobj.nndist;
		
		if((permindex >= 0)) {
		    //contig1=contigobj->getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
		    //	    contigobj->total_conitg_length,perm[permindex]->region[dimensions].start-shiftvalue+1);
		    contig1=contigobj->getrecordnum(contigobj->recordseps+contigobj->offset_array[0],contigobj->number_of_contigs_array[0],
						    contigobj->total_contig_length_array[0],
						    perm[permindex]->region[dimensions].start-shiftvalue+1);
		    //permindex = gapped_nearestNeighbor(elem);
		}
		
		if((permindex < 0)||(perm[permindex]->num==0)) //no such element
		{
		    // start a new chain
		    elem->last = NULL;
		    // the elem as a new root of a new chain
		    // and it is added to the root list
		    if(chainroots==NULL){
			chainroots=(child*)malloc(sizeof(child));
			if(chainroots==NULL)return(-1);
			chainroots->child_elem=elem;
			chainroots->next=NULL;
		    }
		    else{
			tempchild=(child*)malloc(sizeof(child));
			if(tempchild==NULL)return(-1);
			tempchild->child_elem=elem;
			tempchild->next=chainroots;
			chainroots=tempchild;
			
		    }		    		    
		    if(elem->globalscore > maxscore) //new best chain was found
		    {
			maxscore = elem->globalscore; 
			tail = elem; //current last element of the chain
			//maxscore_gapti=gap_ti;
		    }
		}
		// here I have to check if it is in the same contig 
		else if((contig1!=contig2)){ // check for making sure
		    
		    // start a new chain
		    elem->last=NULL; 
		    // the elem as a new root of a new chain
		    // and it is added to the root list
		    if(chainroots==NULL){
			chainroots=(child*)malloc(sizeof(child));
			if(chainroots==NULL)return(-1);
			chainroots->child_elem=elem;
			chainroots->next=NULL;
		    }
		    else{
			tempchild=(child*)malloc(sizeof(child));
			if(tempchild==NULL)return(-1);
			tempchild->child_elem=elem;
			tempchild->next=chainroots;
			chainroots=tempchild;
			
		    }
		    
		    
		}
		else{
		    // both are in the same contig
		    /// geting the gap between elem and perm[permindex]
		    gap_ij=-1*(dimensions+1);;
		    gap_ti=0;
		    for(j=0;j<dimensions;j++){
			gap_ij=gap_ij+elem->region[j].start-perm[permindex]->region[j].end;
			gap_ti=gap_ti+genome_sizes_array[j+1]-perm[permindex]->region[j].end;
			
		    }
		    gap_ij=gap_ij+elem->region[dimensions].start-perm[permindex]->region[dimensions].end;
		    gap_ti=gap_ti+genome_sizes_array[0]-perm[permindex]->region[dimensions].end;
		    
		    //gap_ij=1+(gap_ij/(dimensions+1)); 	
		    tempscore= elem->weight + nndist+gap_ti-gap_ij;
		    maxscore_gapti=tempscore;
		    /// decide if we start a new chain
		    if((tempscore>Threshold)&&(gap_ij<gap_threshold)){
			
			tempscore= elem->globalscore + nndist+(double)gap_ti-(double)gap_ij;
			elem->last = perm[permindex]; //setting chain reference	
			elem->globalscore =tempscore;
			// Adding the fragment to the child list
			if(perm[permindex]->start_child_list==NULL){
			    perm[permindex]->start_child_list=(child*)malloc(sizeof(child));
			    if(perm[permindex]->start_child_list==NULL)return(-1);
			    perm[permindex]->start_child_list->next=NULL;
			    perm[permindex]->start_child_list->child_elem=elem;
			    //perm[permindex]->start_child_list=elem;
			}
			else{
			    tempchild=(child*)malloc(sizeof(child));
			    if(tempchild==NULL)return(-1);
			    tempchild->child_elem=elem;
			    tempchild->next=perm[permindex]->start_child_list;
			    perm[permindex]->start_child_list=tempchild;			    
			}
			
		    }
		    else{
			elem->last=NULL; 
			// the elem as a new root of a new chain
			// and it is added to the root list
			if(chainroots==NULL){
			    chainroots=(child*)malloc(sizeof(child));
			    if(chainroots==NULL)return(-1);
			    chainroots->child_elem=elem;
			    chainroots->next=NULL;
			}
			else{
			    tempchild=(child*)malloc(sizeof(child));
			    if(tempchild==NULL)return(-1);
			    tempchild->child_elem=elem;
			    tempchild->next=chainroots;
			    chainroots=tempchild;
			    
			}
			
		    }
		    if( maxscore_gapti> maxscore) //new best chain was found
		    {
			//maxscore = elem->globalscore; 
			maxscore =max(0,maxscore_gapti);
			tail = elem; //current last element of the chain
			//maxscore_gapti=elem->weight-elem->globalscore;
		    }		    		    
		}				
	    }
	    else{ // we are in a new contig
		
		// Deactivate the elements from lower_limit till i
		for(j=lower_limit;j<i;j++){
		    if(point_type[j]!=1)
			kdtreeobj.DeactivateElem(&blocks[sorted_list[j]]);
		}
		// set the new contig index
		lower_limit=i;
		current_contig=contig2;
		// start a new chain
		elem->last=NULL; 
		// the elem as a new root of a new chain
		// and it is added to the root list
		if(chainroots==NULL){
		    chainroots=(child*)malloc(sizeof(child));
		    if(chainroots==NULL)return(-1);
		    chainroots->child_elem=elem;
		    chainroots->next=NULL;
		}
		else{
		    tempchild=(child*)malloc(sizeof(child));
		    if(tempchild==NULL)return(-1);
		    tempchild->child_elem=elem;
		    tempchild->next=chainroots;
		    chainroots=tempchild;		    
		}
		
		if(maxscore_gapti> maxscore) //new best chain was found
		{
		    maxscore = maxscore_gapti; 
		    tail = elem; //current last element of the chain
		    
		}
	    }
	    
	}
	else{ // we are at the end point of a fragment
	    // insert it in the tree
	    kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
	    if(maxscore_gapti> maxscore) //new best chain was found
	    {
		maxscore = maxscore_gapti; 
		tail = elem; //current last element of the chain
	    }
	    
	}
    }
    delete contig_id_array;
    delete positions;
    free(lower_bound->region);
    delete lower_bound;
    global_maxscore=maxscore;
    return(1);
}

/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////


int kdChainer::check_est(int depth, double th_score,double given_score,
			 ContigSorting* contigobj,double percentagescore){
	
	
	//ContigSorting contigobj;
//	Kdelem_ptr elem = NULL;  //last element of the chain
//	long gap_ti=0;
	Kdelem_ptr depth_ptr=tail;
        //printf("Score %f\n",given_score);
	if((float)given_score<(float)th_score)return(-1);
	for(int i=0;i<depth;i++){
		if(depth_ptr->last==NULL)return(-1);
		depth_ptr=depth_ptr->last;
	}
	return(1);
}


/////////////////////////////////////////////////////////////////////////////
//                    Placement of ESTs on Complete Genomes                //    
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

int kdChainer::mod_EST_SW(double Threshold,char* contigfilename)
{
    //ContigSorting contigobj;
    Kdelem_ptr elem;
//  Kdelem_ptr temp_elem;
    double tempscore;
    double maxscore = MIN_VALUE;
    long permindex;
    long i,j;
//    long gap_ij;
//    long gap_ti;
    child* tempchild;
    
    ///////Contigs-specific Functionality
    long contig1=0;
    long contig2=0;
    long current_contig=0;
//    int contigflag=1;
//    int current_contig_flag=0;
//    double maxscoreincontig=MIN_VALUE-100;
//    double prevmaxscoreincontig=0;
//    long gapElem_ti=0;
//    char firsttimeflag=0;
    long lower_limit;
    int res;
    
    /// gcfile
    long *gc_vector=new long [dimensions+1];
    for(long j=0;j<=dimensions;j++){
	//gc_vector[j]=(MyMAX_VALUE);
	gc_vector[j]=gap_threshold;
    }
    //printf("HALOOOOOOOOOOOOO %lu gc_vec %lu\n",gap_threshold,gc_vector[0]);
    
    if(gcfile!=NULL){	
	char linearr[MaxFragmentLine+1];
	char* input_token=NULL;
	
	FILE* gcfptr=FOPEN(gcfile,"r");
	if(gcfptr==NULL){
	    printf("Error: cannot open file: %s\n",gcfile);
	    delete []gc_vector;
	    return(-1);
	}
	
	fgets(linearr,MaxFragmentLine,gcfptr);	
	//printf("line:      %s\n",linearr);
	if(linearr!=NULL){
	    //  printf("%s",linearr);
	    input_token=strtok(linearr," ");
	    if(input_token!=NULL){
		gc_vector[0]=atol(input_token);
		//printf("haloooooooooooooooooooo %lu\n",gc_vector[0]);
	    }
	    for(long j=1;j<=dimensions;j++){
		input_token=strtok(NULL," ");
		if(input_token==NULL)break;
		gc_vector[j]=atol(input_token);
		//printf("haloooooooooooooooooooo %lu\n",gc_vector[j]);
	    }
	}
	fclose(gcfptr);
	
    }
    gap_threshold=gc_vector[0];
    // gcfile


    Kdelem_ptr lower_bound=new Kdelem;
    lower_bound->region=(Region*)malloc(sizeof(Region)*(dimensions+1));
    if(lower_bound->region==NULL)return(-1);
    
    
    res=contigobj->read_contigs(contigfilename,0);
    if(res<0)return(-1);
    
    elem=&blocks[sorted_list[0]];
    
    elem->last=NULL;
    if(chainroots==NULL){
	chainroots=(child*)malloc(sizeof(child));
	if(chainroots==NULL)return(-1);
	chainroots->child_elem=elem;
	chainroots->next=NULL;
    }
    else{
	tempchild=(child*)malloc(sizeof(child));
	if(tempchild==NULL)return(-1);
	tempchild->next=NULL; // initialization
	tempchild->child_elem=elem;
	tempchild->next=chainroots;
	chainroots=tempchild;
	
    }
    elem->globalscore=elem->weight;
    
    kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
    lower_limit=2;

    for(i = 2; i<(numofblocks*2); i++)
    {	
	elem=&blocks[sorted_list[i]];
	
	if(point_type[i]==1){  // it is start point of a fragment
	    elem->globalscore=elem->weight;
	    // get its right end
	    // add weight to its right end
	    // for gapthreshold
	    j=lower_limit;
	    while(j<i){
		if(point_type[j]!=1){
		    long delta_x=elem->region[dimensions].start-blocks[sorted_list[j]].region[dimensions].end;
		    if(delta_x>gap_threshold){
			kdtreeobj.DeactivateElem(&blocks[sorted_list[j]]);
			//flag=1;
		    }
		    else{
			break;
		    }
		    j++;
		}
		else{
		    j++;
		}				
	    }
	    lower_limit=j;
	    //end for gap threshold


	    /// We also have to check if the two fragments are in the same conitg
	    
	    contig2=contigobj->getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
					    contigobj->total_conitg_length,
					    elem->region[dimensions].end-shiftvalue-1);
	
	    
	    if(current_contig==contig2) // still in the same EST/contig
	    {
		// a modification using the region search function
		// get the lower_bound point	
				
		//lower_bound->region[dimensions].start=(long)max(0,
		//		   blocks[sorted_list[lower_limit]].region[dimensions].start-gap_threshold);
		lower_bound->region[dimensions].start=(long)max(0,elem->region[dimensions].start-gap_threshold);
		lower_bound->region[dimensions].end=lower_bound->region[dimensions].start;  		
		for(j=0;j<dimensions;j++){		    
		    //lower_bound->region[j].start=(long)max(0,elem->region[j].start-gap_threshold);
		    lower_bound->region[j].start=(long)max(0,elem->region[j].start-gc_vector[j+1]);
		    lower_bound->region[j].end=lower_bound->region[j].start;  				    
		}	
		
		permindex=kdtreeobj.region_gapped_nearestNeighbor(elem,lower_bound);					
		
		nndist=kdtreeobj.nndist;
		
		if((permindex >= 0)) {
		    contig1=contigobj->getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
						    contigobj->total_conitg_length,
						    perm[permindex]->region[dimensions].start-shiftvalue+1);
		    //permindex = gapped_nearestNeighbor(elem);
		}
		
		if((permindex < 0)||(perm[permindex]->num==0)) //no such element
		{
		    // start a new chain
		    elem->last = NULL;
		    // the elem as a new root of a new chain
		    // and it is added to the root list
		    if(chainroots==NULL){
			chainroots=(child*)malloc(sizeof(child));
			if(chainroots==NULL)return(-1);
			chainroots->child_elem=elem;
			chainroots->next=NULL;
		    }
		    else{
			tempchild=(child*)malloc(sizeof(child));
			if(tempchild==NULL)return(-1);
			//tempchild->next=NULL; // initialization
			tempchild->child_elem=elem;
			tempchild->next=chainroots;
			chainroots=tempchild;
			
		    }
		    
		    
		    if(elem->globalscore > maxscore) //new best chain was found
		    {
			maxscore = elem->globalscore; 
			tail = elem; //current last element of the chain
		    }
		}
		// here I have to check if it is in the same contig 
		else if((contig1!=contig2)){ // check for making sure
		    
		    // start a new chain
		    elem->last=NULL; 
		    // the elem as a new root of a new chain
		    // and it is added to the root list
		    if(chainroots==NULL){
			chainroots=(child*)malloc(sizeof(child));
			if(chainroots==NULL)return(-1);
			chainroots->child_elem=elem;
			chainroots->next=NULL;
		    }
		    else{
			tempchild=(child*)malloc(sizeof(child));
			if(tempchild==NULL)return(-1);
			tempchild->child_elem=elem;
			tempchild->next=chainroots;
			chainroots=tempchild;
			
		    }
		    
		    
		}
		else{
		    // both are in the same contig
		    /// geting the gap between elem and perm[permindex]
		    /// decide if we start a new chain
		    tempscore= elem->weight + nndist;		    
		    // Debug Mohamed		    
		    if((tempscore>Threshold)){
			tempscore= elem->globalscore + nndist;
			if(overlapping_value>0){ // correct tempscore and remove overlap
			    if(elem->region[dimensions].start-overlapping_value<=perm[permindex]->region[dimensions].end){
				tempscore= tempscore-(perm[permindex]->region[dimensions].end-(elem->region[dimensions].start-overlapping_value)+1);
			    }
				
			}
			elem->last = perm[permindex]; //setting chain reference	
			elem->globalscore =tempscore;
			// Adding the fragment to the child list
			if(perm[permindex]->start_child_list==NULL){
			    perm[permindex]->start_child_list=(child*)malloc(sizeof(child));
			    if(perm[permindex]->start_child_list==NULL)return(-1);
			    perm[permindex]->start_child_list->next=NULL;
			    perm[permindex]->start_child_list->child_elem=elem;
			    //perm[permindex]->start_child_list=elem;
			}
			else{
			    tempchild=(child*)malloc(sizeof(child));
			    if(tempchild==NULL)return(-1);
			    tempchild->child_elem=elem;
			    tempchild->next=perm[permindex]->start_child_list;
			    perm[permindex]->start_child_list=tempchild;			    
			}
			
		    }
		    else{
			elem->last=NULL; 
			// the elem as a new root of a new chain
			// and it is added to the root list
			if(chainroots==NULL){
			    chainroots=(child*)malloc(sizeof(child));
			    if(chainroots==NULL)return(-1);
			    chainroots->child_elem=elem;
			    chainroots->next=NULL;
			}
			else{
			    tempchild=(child*)malloc(sizeof(child));
			    if(tempchild==NULL)return(-1);
			    tempchild->child_elem=elem;
			    tempchild->next=chainroots;
			    chainroots=tempchild;
			    
			}
			
		    }
		    if(elem->globalscore > maxscore) //new best chain was found
		    {
			maxscore = elem->globalscore; 
			tail = elem; //current last element of the chain
		    }		    		    
		}				
	    }
	    else{ // we are in a new contig
		
		// Deactivate the elements from lower_limit till i
		for(j=lower_limit;j<i;j++){
		    if(point_type[j]!=1)
			kdtreeobj.DeactivateElem(&blocks[sorted_list[j]]);
		}
		// set the new contig index
		lower_limit=i;
		current_contig=contig2;
		// start a new chain
		elem->last=NULL; 
		// the elem as a new root of a new chain
		// and it is added to the root list
		if(chainroots==NULL){
		    chainroots=(child*)malloc(sizeof(child));
		    if(chainroots==NULL)return(-1);
		    chainroots->child_elem=elem;
		    chainroots->next=NULL;
		}
		else{
		    tempchild=(child*)malloc(sizeof(child));
		    if(tempchild==NULL)return(-1);
		    tempchild->child_elem=elem;
		    tempchild->next=chainroots;
		    chainroots=tempchild;		    
		}
		
		if(elem->globalscore > maxscore) //new best chain was found
		{
		    maxscore = elem->globalscore; 
		    tail = elem; //current last element of the chain
		}
	    }
	    
	}
	else{ // we are at the end point of a fragment
	    // insert it in the tree
	    kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
	    if(elem->globalscore > maxscore) //new best chain was found
	    {
		maxscore = elem->globalscore; 
		tail = elem; //current last element of the chain
	    }
	    
	}
    }
    
    free(lower_bound->region);
    delete lower_bound;
    global_maxscore=maxscore;
    delete[] gc_vector;
    return(1);

}

/////////////////////////////////////////////////////////////////////////////
// Function: storeMultiEST
// The score of cDNA/ESTs does not include gap constraints
// If perc_coverage_flag==1 then the contigs
// are sorted with the percentage coverage score

void kdChainer::storeMultiEST(ProcessFragmentFile* InfoObj,int IncludeDirection)
{
    
  Kdelem_ptr elem = NULL;  //last element of the chain
  long i,j;
  FILE* fptr=NULL;
  FILE* clusterfptr=NULL;
  FILE* compactchainfptr=NULL;
  child* ptr;
  char temp[MaxFileNameSize];
//  char tempstr[MaxFileNameSize];
  long* regionend=(long*)malloc(sizeof(long)*(dimensions+1));
  int flag=0;
  long chain_counter=0;
 
  


  sprintf(temp,"%s.chn",InfoObj->filename);
  double maxscore = MIN_VALUE-10;

  if(chainroots==NULL)return;
  fptr=FOPEN(temp,"w");
  
  // opening cluster file
  if(cluster_file_flag){
      sprintf(temp,"%s.cst",InfoObj->filename);     
      clusterfptr=FOPEN(temp,"w");
  }
  sprintf(temp,"%s.ccn",InfoObj->filename);
  //Kdelem_ptr local_tail=NULL;
  compactchainfptr=FOPEN(temp,"w");
  fprintf(fptr,">CHA %d\n",(int)dimensions+1);
  fprintf(compactchainfptr,">CHA %d\n",(int)dimensions+1);
  
  for(ptr=chainroots;ptr!=NULL ; ptr = ptr->next)  //last references the previous elem
  {
      // putting the clusters in file
      
      tail=NULL;
      maxscore = MIN_VALUE-10;
      if(cluster_file_flag){
	  if(ptr->child_elem->num!=0) fprintf(clusterfptr,"#cluster\n");
      }
      if(cluster_file_flag)traverse_store_EST(ptr,clusterfptr,&maxscore);
      else silent_traverse_cluster(ptr,&maxscore);
      long gap_ti=0;
      for(j=0;j<dimensions;j++){
	  gap_ti=gap_ti+genome_sizes_array[j+1]-tail->region[j].end;
	  
      }
      gap_ti=gap_ti+genome_sizes_array[0]-tail->region[dimensions].end;
      
      long contig1=contigobj->getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
					   contigobj->total_conitg_length,
					   tail->region[dimensions].end-shiftvalue-1);
      long contiglength;
      if(contig1==0){
	  contiglength=(long)(contigobj->recordseps[contig1]);
	  
      }
      else{ 
	  contiglength=(long)(contigobj->recordseps[contig1]-contigobj->recordseps[contig1-1])-1;
      }
      
      if(check_est(chain_depth,chain_score,tail->globalscore,contigobj,percentagescore)!=1)continue;
      // mohamed 8.11.2006. checkchain
      //tmpcheckchain=tail->globalscore
      if(chain_average_length>0){
	  if(check_chain(ptr->child_elem,chain_depth,chain_score,tail->globalscore)!=1)continue;
      }

      double perc_coverage=tail->globalscore/(double)(contiglength*(double)multiplyingfactor);

      //      printf("coordinate= %d   \n",(tail->region[dimensions].end-shiftvalue)-1);
      // printf("contig number= %d   \n",contig1);
      // printf("perc flag= %f   \n",perc_coverage);
      // printf("contig length %d\n",contiglength);
      if((perc_coverage<percentagescore)&&(perc_cover_flag==1)){
	  continue;
      }
	//if(ptr->child_elem->num!=0){
	if(perc_cover_flag==1){
	    if(ptr->child_elem->num!=0)
		fprintf(fptr, "#%lf coverage %lf\n",((tail->globalscore)/(double)multiplyingfactor),perc_coverage);
	    if(ptr->child_elem->num!=0) 
		fprintf(compactchainfptr, "#%lf coverage %lf\n",((tail->globalscore)/(double)multiplyingfactor),perc_coverage);
	}
	else{
	    if(ptr->child_elem->num!=0)
	    fprintf(fptr, "#%lf\n",((tail->globalscore)/(double)multiplyingfactor));
	    if(ptr->child_elem->num!=0)
	    fprintf(compactchainfptr, "#%lf\n",((tail->globalscore)/(double)multiplyingfactor));
	}		
	flag=0;	    
	
	for(elem=tail;elem;elem=elem->last){
	    
	  if (elem->num!=0){
	      fprintf(fptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,\
		      elem->region[dimensions].end-shiftvalue);
	  }
	  
	  for(i=0;i<dimensions;i++){
	      if (elem->num!=0){
		  fprintf(fptr,"[%lu,%lu] ", elem->region[i].start-shiftvalue,elem->region[i].end-shiftvalue);
	      }
	    		  
	  }
	  if(elem->num!=0){
	      fprintf(fptr,"%lu",(long)(elem-blocks)); //mohamed
	      fprintf(fptr,"\n");
	  }
	  //// storing compact chain
	  if(flag==0){
	      
	      regionend[0]=elem->region[dimensions].end;
	      flag=1;
	      for(j=0;j<dimensions;j++){		  		  
		  regionend[j+1]=elem->region[j].end;
	      }	
	      // Chain of 1 fragment
	      if(elem->last==NULL){
		  if (elem->num!=0){
		      fprintf(compactchainfptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,\
			      elem->region[dimensions].end-shiftvalue);
		      chain_counter++;
		  }
		  
		  for(i=0;i<dimensions;i++){
		      if (elem->num!=0){
			  fprintf(compactchainfptr,"[%lu,%lu] ", elem->region[i].start-shiftvalue,
				  elem->region[i].end-shiftvalue);
		      }		  
		  }
		  if(elem->num!=0) fprintf(compactchainfptr,"\n");
	      
	      }
	      
	  }
	  else if(elem->last==NULL){
	      if(elem->num!=0){
		  fprintf(compactchainfptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,
			  regionend[0]-shiftvalue);		  
		  chain_counter++;
	      }
	      
	      for(j=0;j<dimensions;j++){
		  if(elem->num!=0){
		      fprintf(compactchainfptr,"[%lu,%lu] ", elem->region[j].start-shiftvalue,
			      regionend[j+1]-shiftvalue);		  
		  }
	      }	
	      if(elem->num!=0)fprintf(compactchainfptr,"\n");
	      
	  }
	  ////
	  
	  
	}
	
	
  }
  fclose(fptr);
  if(cluster_file_flag)fclose(clusterfptr);
  fclose(compactchainfptr);
  free(regionend);
  // Append information to the statistics file .stc
      sprintf(temp,"%s.stc",InfoObj->filename);     
      FILE* stcptr=FOPEN(temp,"a");
      if(stcptr==NULL){
	  printf("Cannot append to statistics file\n");
	  return;
      }
      printf("# no. of chains: %lu \n",chain_counter);
      fprintf(stcptr,"# no. of chains: %lu \n",chain_counter);
      fclose(stcptr);
//
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// In this function, the score of the chain will be
// reported without considering the gaps between the fragments

Kdelem_ptr kdChainer::traverse_store_EST(child* root, FILE* fptr, double* maxscoreptr){
	
	Kdelem_ptr elem;	
	child* child_ptr;
	int i;
	
	if(root==NULL)return(NULL);
	if(root->child_elem==NULL)return(NULL);
	elem=root->child_elem;
	//tail=elem;
	// here I can process the element
	if(elem->num!=0){	    
	    fprintf(fptr,"[%lu,%lu] ", elem->region[dimensions].start-shiftvalue,
		    elem->region[dimensions].end-shiftvalue);
	}
	if(elem->num!=0){
	    for(i=0;i<dimensions;i++)
	    {	    
		fprintf(fptr,"[%lu,%lu] ", elem->region[i].start-shiftvalue,elem->region[i].end-shiftvalue);
	    }
	  
	}
	if(elem->num!=0){
	    fprintf(fptr,"%lu",(long)(elem-blocks)); //mohamed
	    fprintf(fptr,"\n");
	}
	
	long gap_ti=0;
	for(int j=0;j<dimensions;j++){
	    gap_ti=gap_ti+genome_sizes_array[j+1]-elem->region[j].end;
	    
	}
	gap_ti=gap_ti+genome_sizes_array[0]-elem->region[dimensions].end;
	
	/*
	// Modified for EST
	
	gap_ti=0;
	if(elem->last==NULL) // it is a root
	{
	    elem->globalscore=elem->weight;
	}
	else{
	    elem->globalscore=elem->last->globalscore+elem->weight;
	}
	*/
	if(elem->globalscore> *maxscoreptr) //new best chain was found
	{
	    // Modification to correct EST score
	    //*maxscoreptr = elem->globalscore+gap_ti;
	    *maxscoreptr = elem->globalscore;
	    tail = elem; //current last element of the chain
	}
	child_ptr=elem->start_child_list;
	while(child_ptr!=NULL){
	    traverse_store_EST(child_ptr,fptr,maxscoreptr);
	    child_ptr=child_ptr->next;
	    
	}
	return(NULL);
}

//////////////////////////////////////////////////////////////////////
/////////////////////////// storing in database format ///////////////
//////////////////////////////////////////////////////////////////////
// global chains
/////////////////////////////////////////////////////////////////////

void kdChainer::storeChainDB(ProcessFragmentFile* InfoObj,int gappedflag)
{
  Kdelem_ptr elem = tail;  //last element of the chain
  long i;
  FILE* fptr=NULL;
  char temp[MaxFileNameSize];

//  char tempstr[MaxFileNameSize];
  double leastscore=0;
  int scorenumdigit;

  int* positiondigitarray=(int*)malloc(sizeof(long)*(dimensions+1));
  for(i=0;i<=dimensions;i++){
      positiondigitarray[i]=2+(long)log10((double)genome_sizes_array[i]);
      leastscore=leastscore+((double)genome_sizes_array[i]);
  }
  int fragmentnumdigits=2+(int)log10((double)numofblocks);
  
  leastscore=max(moduls(global_maxscore),leastscore);
  
  if(gappedflag)
      scorenumdigit=7+1+1+(long)log10(leastscore); // on digit for negative scores
  else
      scorenumdigit=7+1+(long)log10(global_maxscore); // on digit for negative scores
      

  sprintf(temp,"%s.dbf",InfoObj->filename);
  if(stdout_flag!=1){
      fptr=FOPEN(temp,"w");
      //fprintf(fptr,">CHA %d\n",(int)dimensions+1);
  }

  // adding the score of the terminus
  //tail->globalscore=tail->globalscore+multiplyingfactor;
  if(stdout_flag!=1){
      fprintf(fptr,"# format: id s1 e1 ... sk ek score\n");
      fprintf(fptr,"# id: fragment id \n");
      fprintf(fptr,"# s1: fragment's start position in genome 1 \n");
      fprintf(fptr,"# e1: fragment's end position in genome 1 \n");
      fprintf(fptr,"# sk: fragment's start position in genome k \n");
      fprintf(fptr,"# ek: fragment's end position in genome k \n");
      fprintf(fptr,"# score: score of the chain ending at this fragment \n");
      fprintf(fptr,"# num of genomes: %d\n",(int)dimensions+1);
      fprintf(fptr,"# total num. of fragments including origin and terminus: %lu\n",numofblocks+1);
      fprintf(fptr,"# global chain score: %lf\n",(tail->globalscore+multiplyingfactor));  
  }
  else{
      printf("# format: id s1 e1 ... sk ek score\n");
      printf("# id: fragment id \n");
      printf("# s1: fragment's start position in genome 1 \n");
      printf("# e1: fragment's end position in genome 1 \n");
      printf("# sk: fragment's start position in genome k \n");
      printf("# ek: fragment's end position in genome k \n");
      printf("# score:  score of the chain ending at this fragment \n");
      printf("# num. of genomes: %d\n",(int)dimensions+1);
      printf("# total num. of fragments including origin and terminus: %lu\n",numofblocks+1);
      printf("# global chain score: %lf\n",(tail->globalscore+multiplyingfactor));
  }
  // storing the terminus
  if(stdout_flag!=1)
  fprintf(fptr,"%*lu",fragmentnumdigits,numofblocks);
  else printf("%*lu",fragmentnumdigits,numofblocks);
//  long gti=0;
  for(i=0;i<=dimensions;i++)
  {
      if(stdout_flag!=1)fprintf(fptr,"  %*lu  %*lu",positiondigitarray[i],genome_sizes_array[i],
	      positiondigitarray[i],genome_sizes_array[i]);
      else printf("  %*lu  %*lu",positiondigitarray[i],genome_sizes_array[i],
	      positiondigitarray[i],genome_sizes_array[i]);
  }
  if(stdout_flag!=1)
      fprintf(fptr,"  %*.6lf\n",scorenumdigit,(tail->globalscore+multiplyingfactor));
  else  printf("  %*.6lf\n",scorenumdigit,(tail->globalscore+multiplyingfactor));
  //printf("%f\n",(float)gti);
  
  // storing the rest of the chain in a file
  if(stdout_flag!=1){
      for(elem=tail; elem; elem = elem->last)  {//last references the previous elem
	    // inverse the transform before storing the points
	  
	  if(elem->num!=0){	   
	      if(overlapping_value>0)inverse_transform_for_overlapping(elem);
	  }
	  
	  fprintf(fptr,"%*lu",fragmentnumdigits,(long)(elem-blocks));
	  if (elem->last!=NULL){
	      fprintf(fptr,"  %*lu  %*lu", positiondigitarray[0],elem->region[dimensions].start-shiftvalue,\
		      positiondigitarray[0],elem->region[dimensions].end-shiftvalue);
	  }
	  else{
	      fprintf(fptr,"  %*lu  %*lu", positiondigitarray[0],elem->region[dimensions].start,
		      positiondigitarray[0],elem->region[dimensions].end);
	  }
	  for(i=0;i<dimensions;i++)
	  {
	      if (elem->last!=NULL){
		  fprintf(fptr,"  %*lu  %*lu", positiondigitarray[i+1],elem->region[i].start-shiftvalue,\
			  positiondigitarray[i+1],elem->region[i].end-shiftvalue);
	      }
	      else{
		  fprintf(fptr,"  %*lu  %*lu", positiondigitarray[i+1],elem->region[i].start,\
			  positiondigitarray[i+1],elem->region[i].end);
	      }
	  }
	  
	  fprintf(fptr,"  %*.6lf\n",scorenumdigit,elem->globalscore);
      }
      fclose(fptr);
  }
  else{
      for(elem=tail; elem; elem = elem->last)  {//last references the previous elem
	  if(elem->num!=0){	   
	      if(overlapping_value>0)inverse_transform_for_overlapping(elem);
	  }
	  printf("%*lu",fragmentnumdigits,(long)(elem-blocks));
	  if (elem->last!=NULL){
	      printf("  %*lu  %*lu", positiondigitarray[0],elem->region[dimensions].start-shiftvalue,\
		      positiondigitarray[0],elem->region[dimensions].end-shiftvalue);
	  }
	  else{
	      printf("  %*lu  %*lu", positiondigitarray[0],elem->region[dimensions].start,
		      positiondigitarray[0],elem->region[dimensions].end);
	  }
	  for(i=0;i<dimensions;i++)
	  {
	      if (elem->last!=NULL){
		  printf("  %*lu  %*lu", positiondigitarray[i+1],elem->region[i].start-shiftvalue,\
			  positiondigitarray[i+1],elem->region[i].end-shiftvalue);
	      }
	      else{
		  printf("  %*lu  %*lu", positiondigitarray[i+1],elem->region[i].start,\
			  positiondigitarray[i+1],elem->region[i].end);
	      }
	  }
	  
	  printf("  %*.6lf\n",scorenumdigit,elem->globalscore);
      }
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////

void kdChainer::storeMultiChainDB(ProcessFragmentFile* InfoObj,int IncludeDirection){
    Kdelem_ptr elem = NULL;  //last element of the chain
    long i;
    FILE* fptr=NULL;
    FILE* clusterfptr=NULL;
    
    child* ptr;
    long j;
    char temp[MaxFileNameSize];
    //char tempstr[MaxFileNameSize];
    long* regionend=(long*)malloc(sizeof(long)*(dimensions+1));
    long counter=1;
    int flag=0;
    
    double maxscore = MIN_VALUE-10;
    double chscore=0;
    
    if(chainroots==NULL)return;
    
    if(stdout_flag!=1){
	sprintf(temp,"%s.dbf",InfoObj->filename);    
	fptr=FOPEN(temp,"w");
    }
    if(cluster_file_flag){
	// opening cluster file
	sprintf(temp,"%s.cst",InfoObj->filename);   
	clusterfptr=FOPEN(temp,"w");
    }
    // opening compact chain file
 
    if(stdout_flag!=1){
	//fprintf(fptr,">CHA %d\n",(int)dimensions+1); 
	fprintf(fptr,"# format: ch_id f_id s1 e1 ... sk ek score\n");
	fprintf(fptr,"# ch_id: chain's id \n");
	fprintf(fptr,"# f_id: fragment's id \n");
	fprintf(fptr,"# s1: fragment's start position in genome 1 \n");
	fprintf(fptr,"# e1: fragment's end position in genome 1 \n");
	fprintf(fptr,"# sk: fragment's start position in genome k \n");
	fprintf(fptr,"# ek: fragment's end position in genome k \n");
	fprintf(fptr,"# score:  score of the chain \n");
	fprintf(fptr,"# num. of genomes: %d\n",(int)dimensions+1);
	fprintf(fptr,"# total num. of fragments: %lu\n",numofblocks-1);
	if(global_maxscore>=0)
	    fprintf(fptr,"# max score of a local chain: %lf\n",global_maxscore);
    }
    else{
	
	printf("# format: ch_id f_id s1 e1 ... sk ek s\n");
	printf("# ch_id: chain's id \n");
	printf("# f_id: fragment's id \n");
	printf("# s1: fragment's start position in genome 1 \n");
	printf("# e1: fragment's end position in genome 1 \n");
	printf("# sk: fragment's start position in genome k \n");
	printf("# ek: fragment's end position in genome k \n");
	printf("# score:  score of the chain \n");
	printf("# num. of genomes: %d\n",(int)dimensions+1);
	printf("# total num. of fragments: %lu\n",numofblocks-1);
	if(global_maxscore>=0)
	    printf("# max score of a local chain: %lf\n",global_maxscore);
    }
    // computing number of digits to be displayed

    int* positiondigitarray=(int*)malloc(sizeof(long)*(dimensions+1));
    for(i=0;i<=dimensions;i++){
	positiondigitarray[i]=2+(long)log10((double)genome_sizes_array[i]);
    }
    int fragmentnumdigits=2+(long)log10((double)numofblocks);
    int scorenumdigit=7+1+(long)log10(max(maxscore_gapti,max_fragment_weight));
 
    for(ptr=chainroots;ptr!=NULL ; ptr = ptr->next)  //last references the previous elem
    {
	// putting the clusters in file
	// CHANGE: Allow Storing chains of one element
	//if(ptr->child_elem->start_child_list==NULL) continue;
	flag=0;
	tail=NULL;
	maxscore = MIN_VALUE-10;

	if(cluster_file_flag){
	    if(ptr->child_elem->num!=0)fprintf(clusterfptr,"#cluster\n");
	    // the following function store clusters of chains and returns a pointer
	    // to the opimal chain in this cluster
	    traverse_store_chain(ptr,clusterfptr,&maxscore);
	}
	else{
	    silent_traverse_cluster(ptr,&maxscore);
	}
	// putting the optimal chain of cluster in other file
	long gap_ti=0;
	for(j=0;j<dimensions;j++){
	    gap_ti=gap_ti+genome_sizes_array[j+1]-tail->region[j].end;	    
	}
	gap_ti=gap_ti+genome_sizes_array[0]-tail->region[dimensions].end;
	
	double tmpcheckchain;
	// num =0 means that this fragment is the ORIGIN
	if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
	    tmpcheckchain=tail->globalscore+gap_ti;	    
	}
	else{	    	    
	    if (tail->last!=NULL) // the chain is not ending with the origin
		tmpcheckchain=tail->globalscore+gap_ti+(shiftvalue*(dimensions+1));
	    else tmpcheckchain=tail->globalscore+gap_ti;	    
	}
	
	// check chain against the options: score and chain depth
	chain_score=max(chain_score,global_maxscore*(percentagescore/100)); 
	if(check_chain(ptr->child_elem,chain_depth,chain_score,tmpcheckchain)!=1)continue;       
	
	if(stdout_flag!=1){
	    // store the fragments of the chain	
	    for(elem=tail;elem;elem=elem->last){
		if(elem->num==0)continue;
		
		// storing chain id
		fprintf(fptr,"%*lu",fragmentnumdigits,counter);
		// storing fragment id id
// bookmark
		//	fprintf(fptr,"  %*lu",fragmentnumdigits,(int)(perm[elem->num]-&blocks[0]));
		fprintf(fptr,"  %*lu",fragmentnumdigits,(long)(elem-blocks));
		if (elem->num!=0){
		    fprintf(fptr,"  %*lu  %*lu", positiondigitarray[0],
			    elem->region[dimensions].start-shiftvalue,
			    positiondigitarray[0],elem->region[dimensions].end-shiftvalue);
		}
		
		for(i=0;i<dimensions;i++){
		    if (elem->num!=0){
			fprintf(fptr,"  %*lu  %*lu", positiondigitarray[i+1],elem->region[i].start-shiftvalue,
				positiondigitarray[i+1],elem->region[i].end-shiftvalue);
		    }      
		}	
		// storing the score    
		long gap_ti=0;
		for(j=0;j<dimensions;j++){
		    gap_ti=gap_ti+genome_sizes_array[j+1]-elem->region[j].end;
		    
		}
		gap_ti=gap_ti+genome_sizes_array[0]-elem->region[dimensions].end;
		
		
		if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
		    if(flag==0){
			fprintf(fptr, "  %*.6lf",scorenumdigit,(tail->globalscore+gap_ti));
			chscore =tail->globalscore+gap_ti;
			flag=1;
		    }
		    else{			
			fprintf(fptr,"  %*.6lf",scorenumdigit,chscore);
		    }
		}
		else{
		    if (elem->last!=NULL){ // tail->last is the ORIGIN, the shiftvalue should be considered
			fprintf(fptr, "  %*.6lf\n",scorenumdigit,
				(tail->globalscore+gap_ti+(shiftvalue*(dimensions+1))));		
		    }		 
		}
		if (elem->num!=0) fprintf(fptr,"\n");
		
	    }
	}
	else{ // store in stdout	    // store the fragments of the chain	
	    for(elem=tail;elem;elem=elem->last){
		if(elem->num==0)continue;
		
		// storing chain id
		printf("%*lu",fragmentnumdigits,counter);
		// storing fragment id id
		printf("  %*lu",fragmentnumdigits,(long)(elem-blocks));
		if (elem->num!=0){
		    printf("  %*lu  %*lu", positiondigitarray[0],elem->region[dimensions].start-shiftvalue,
			    positiondigitarray[0],elem->region[dimensions].end-shiftvalue);
		}
		
		for(i=0;i<dimensions;i++){
		    if (elem->num!=0){
			printf("  %*lu  %*lu", positiondigitarray[i+1],elem->region[i].start-shiftvalue,
				positiondigitarray[i+1],elem->region[i].end-shiftvalue);
		    }      
		}	
		// storing the score    
		long gap_ti=0;
		for(j=0;j<dimensions;j++){
		    gap_ti=gap_ti+genome_sizes_array[j+1]-elem->region[j].end;
		    
		}
		gap_ti=gap_ti+genome_sizes_array[0]-elem->region[dimensions].end;
		
		
		if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
		    if(flag==0){
			printf("  %*.6lf",scorenumdigit,
			       (tail->globalscore+gap_ti));
			chscore =tail->globalscore+gap_ti;
			flag=1;
		    }
		    else{
			
			printf("  %*.6lf",scorenumdigit,chscore);
		    }
		}
		else{
		    if (elem->last!=NULL){ // tail->last is the ORIGIN, the shiftvalue should be considered
			printf("  %*.6lf\n",scorenumdigit,
			       (tail->globalscore+gap_ti+(shiftvalue*(dimensions+1))));		
		    }		 
		}
		if (elem->num!=0) printf("\n");		
	    }	    
	}
	    
	counter++;
	    
	    
    }
    if(stdout_flag!=1)fclose(fptr);
    if(cluster_file_flag)fclose(clusterfptr);
    free(regionend);
}
//////////////////////////////////////////////////////////////////////
void kdChainer::storeMultiContigDB(ProcessFragmentFile* InfoObj,int IncludeDirection){
    Kdelem_ptr elem = NULL;  //last element of the chain
    long i;
    FILE* fptr=NULL;
    FILE* clusterfptr=NULL;
        
    child* ptr;
    long j;
    char temp[MaxFileNameSize];
    //char tempstr[MaxFileNameSize];
    
    long counter=1;
    int flag=0;
 
   //Kdelem_ptr local_tail=NULL;
    double maxscore = MIN_VALUE-10;
    double chscore= MIN_VALUE-10;

    if(chainroots==NULL)return;
   
    if(stdout_flag!=1){
	sprintf(temp,"%s.dbf",InfoObj->filename);
	fptr=FOPEN(temp,"w");
    }
    
    // opening cluster file
    if(cluster_file_flag){
	sprintf(temp,"%s.cst",InfoObj->filename);    
	clusterfptr=FOPEN(temp,"w");
    }
    // opening compact chain file
 
    if(stdout_flag!=1){
	//fprintf(fptr,">CHA %d\n",(int)dimensions+1);    
	if(!absolute_flag)
	    fprintf(fptr,"# format: ch_id f_id ctg_id s1 e1 s2 e2 score\n");
	else
	    fprintf(fptr,"# format: ch_id f_id  s1 e1 s2 e2 score\n");
	fprintf(fptr,"# ch_id: chain's id \n");
	fprintf(fptr,"# f_id: fragment's id \n");
	if(!absolute_flag)
	fprintf(fptr,"# ctg_id: contig's id \n");
	fprintf(fptr,"# s1: fragment's start position in draft genome \n");
	fprintf(fptr,"# e1: fragment's end position in draft genome  \n");
	fprintf(fptr,"# s2: fragment's start position in finished genome \n");
	fprintf(fptr,"# e2: fragment's end position in finished genome  \n");
	fprintf(fptr,"# score:  score of the chain\n");
	fprintf(fptr,"# num. of genomes: %d\n",(int)dimensions+1);
	fprintf(fptr,"# total num. of fragments: %lu\n",numofblocks-1);
	if(global_maxscore>=0)
	    fprintf(fptr,"# max score of a local chain: %lf\n",global_maxscore);
    }
    else{
	if(!absolute_flag)
	    printf("# format: ch_id f_id ctg_id s1 e1 s2 e2 score\n");
	else
	    printf("# format: ch_id f_id s1 e1 s2 e2 score\n");
	printf("# ch_id: chain's id \n");
	printf("# f_id: fragment's id \n");
	if(!absolute_flag)
	    printf("# ctg_id: contig's id \n");
	printf("# s1: fragment's start position in draft genome \n");
	printf("# e1: fragment's end position in draft genome  \n");
	printf("# s2: fragment's start position in finished genome \n");
	printf("# e2: fragment's end position in finished genome \n");
	printf("# score:  score of the chain\n");
	printf("# num. of genomes: %d\n",(int)dimensions+1);
	printf("# total num. of fragments: %lu\n",numofblocks-1);
	if(global_maxscore>=0)
	    printf("# max score of a local chain: %lf\n",global_maxscore);
	//printf("# max score of a chain: %lf\n",global_maxscore);		
    }
    // computing the number of digits

    int* positiondigitarray=(int*)malloc(sizeof(long)*(dimensions+1));
    for(i=0;i<=dimensions;i++){
	positiondigitarray[i]=1+(long)log10((double)genome_sizes_array[i]);
    }
    int fragmentnumdigits=1+(int)log10((double)numofblocks);
    int contignumdigits=1+(int)log10(contigobj->number_of_contigs);
    int scorenumdigit=7+1+(int)log10(max(maxscore_gapti,max_fragment_weight));
    

    Kdelem elem2;
    elem2.region=(Region*)malloc(sizeof(Region)*(dimensions+1));
    
    for(ptr=chainroots;ptr!=NULL ; ptr = ptr->next)  //last references the previous elem
    {
	// putting the clusters in file
	// CHANGE: Allow Storing chains of one element
	//if(ptr->child_elem->start_child_list==NULL) continue;
	flag=0;
	tail=NULL;
	maxscore = MIN_VALUE-10;
	if(cluster_file_flag){
	    if(ptr->child_elem->num!=0)fprintf(clusterfptr,"#cluster\n");
	    // the following function store clusters of chains and returns a pointer
	    // to the opimal chain in this cluster	
	    traverse_store_chain(ptr,clusterfptr,&maxscore);
	}
	else{
	    silent_traverse_cluster(ptr,&maxscore);
	}
	// putting the optimal chain of cluster in other file
	long gap_ti=0;
	for(j=0;j<dimensions;j++){
	    gap_ti=gap_ti+genome_sizes_array[j+1]-tail->region[j].end;
	    
	}
	gap_ti=gap_ti+genome_sizes_array[0]-tail->region[dimensions].end;
	
	double tmpcheckchain;
	// num =0 means that this fragment is the ORIGIN
	if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
	    tmpcheckchain=tail->globalscore+gap_ti;
	    
	}
	else{	    	    
	    if (tail->last!=NULL) // the chain is not ending with the origin
		tmpcheckchain=tail->globalscore+gap_ti+(shiftvalue*(dimensions+1));
	    else tmpcheckchain=tail->globalscore+gap_ti;	    
	}
	
	// check chain against the options: score and chain depth
	chain_score=max(chain_score,global_maxscore*(percentagescore/100)); 
	if(check_chain(ptr->child_elem,chain_depth,chain_score,tmpcheckchain)!=1)continue;       
	

	// store the fragments of the chain
	
	if(stdout_flag!=1){
	    for(elem=tail;elem;elem=elem->last){
		if(elem->num==0)continue;
		// storing chain id
		fprintf(fptr,"%*lu",fragmentnumdigits,counter);
		// storing fragment id id
		fprintf(fptr,"  %*lu",fragmentnumdigits,(long)(elem-blocks));
		for(i=0;i<dimensions;i++){
		    elem2.region[i+1].start=elem->region[i].start-shiftvalue;
		    elem2.region[i+1].end=elem->region[i].end-shiftvalue;		
		}
		elem2.region[0].start=elem->region[dimensions].start-shiftvalue;
		elem2.region[0].end=elem->region[dimensions].end-shiftvalue;
		
		if (elem->num!=0){
		    if(!absolute_flag){
			//long ContigId=contigobj->TransformKdelemtoRef(&elem2);
			long * contig_id_array=new long [dimensions+1];
			for(int kk=0;kk<=dimensions;kk++){
			    contig_id_array[kk]=0;
			}
			contigobj->TransformKdelemtoRefGeneralized(&elem2,contig_id_array);
			for(int kk=0;kk<=dimensions;kk++){
			    fprintf(fptr,"  %*lu",contignumdigits,contig_id_array[kk]+1);
			}
                        //fprintf(fptr,"  %*lu",contignumdigits,ContigId+1);
			delete contig_id_array;
		    }      		
		    for(i=0;i<=dimensions;i++){
		    
			fprintf(fptr,"  %*lu  %*lu", positiondigitarray[i],elem2.region[i].start,
				positiondigitarray[i],elem2.region[i].end);
		    }      
		}	
		// storing the score    
		long gap_ti=0;
		for(j=0;j<dimensions;j++){
		    gap_ti=gap_ti+genome_sizes_array[j+1]-elem->region[j].end;
		    
		}
		gap_ti=gap_ti+genome_sizes_array[0]-elem->region[dimensions].end;
		
		
		if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
		    if(flag==0){
			fprintf(fptr, "  %*.6lf",scorenumdigit,(tail->globalscore+gap_ti));
			chscore =tail->globalscore+gap_ti;
			flag=1;
		    }
		    else{
			
			fprintf(fptr, "  %*.6lf",scorenumdigit,chscore);
		    }
		}
		else{
		    if (elem->last!=NULL){ // tail->last is the ORIGIN, the shiftvalue should be considered
			fprintf(fptr, "%*.6lf",scorenumdigit,(tail->globalscore+gap_ti+(shiftvalue*(dimensions+1))));
		    }		 
		}
		if (elem->num!=0) fprintf(fptr,"\n");
		
	    }
	    
	}
	else{ // output in stdoutput
	    for(elem=tail;elem;elem=elem->last){
		if(elem->num==0)continue;
		// storing chain id
		printf("%*lu",fragmentnumdigits,counter);
		// storing fragment id id
		printf("  %*lu",fragmentnumdigits,(long)(elem-blocks));
		for(i=0;i<dimensions;i++){
		    elem2.region[i+1].start=elem->region[i].start-shiftvalue;
		    elem2.region[i+1].end=elem->region[i].end-shiftvalue;		
		}
		elem2.region[0].start=elem->region[dimensions].start-shiftvalue;
		elem2.region[0].end=elem->region[dimensions].end-shiftvalue;
		
		if (elem->num!=0){
		    if(!absolute_flag){
			//long ContigId=contigobj->TransformKdelemtoRef(&elem2);
			long * contig_id_array=new long [dimensions+1];
			for(int kk=0;kk<=dimensions;kk++){
			    contig_id_array[kk]=0;
			}
			contigobj->TransformKdelemtoRefGeneralized(&elem2,contig_id_array);
			for(int kk=0;kk<=dimensions;kk++){
			    printf("  %*lu",contignumdigits,contig_id_array[kk]+1);
			}
			//printf("  %*lu",contignumdigits,ContigId+1);
			delete contig_id_array;
		    }		    				
		    for(i=0;i<=dimensions;i++){
		   
			printf("  %*lu  %*lu", positiondigitarray[i],elem2.region[i].start,
			       positiondigitarray[i],elem2.region[i].end);
		    }      
		}	
		// storing the score    
		long gap_ti=0;
		for(j=0;j<dimensions;j++){
		    gap_ti=gap_ti+genome_sizes_array[j+1]-elem->region[j].end;
		    
		}
		gap_ti=gap_ti+genome_sizes_array[0]-elem->region[dimensions].end;
		
		
		if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
		    if(flag==0){
			printf("  %*.6lf",scorenumdigit,(tail->globalscore+gap_ti));
			chscore =tail->globalscore+gap_ti;
			flag=1;
		    }
		    else{		       		       
			printf("  %*.6lf",scorenumdigit,chscore);
		    }
		}
		else{
		    if (elem->last!=NULL){ // tail->last is the ORIGIN, the shiftvalue should be considered
			printf("  %*.6lf",scorenumdigit,(tail->globalscore+gap_ti+(shiftvalue*(dimensions+1))));
			
		    }		 
		}
		if (elem->num!=0) printf("\n");
		
	    }
	}
	counter++;
	
	
    }
    if(stdout_flag!=1)fclose(fptr);
    if(cluster_file_flag)fclose(clusterfptr);
    free(elem2.region);
}
//////////////////////////////////////////////////////////////////////
// Storing EST in dbformat
//////////////////////////////////////////////////////////////////////

void kdChainer::storeMultiEstDB(ProcessFragmentFile* InfoObj,int IncludeDirection){
    Kdelem_ptr elem = NULL;  //last element of the chain
    long i;
    FILE* fptr=NULL;
    FILE* clusterfptr=NULL;
    
    child* ptr;
    
    char temp[MaxFileNameSize];
    //char tempstr[MaxFileNameSize];
    
    long counter=1;
    
    
    
    double maxscore = MIN_VALUE-10;
    
    if(chainroots==NULL)return;
    if(stdout_flag!=1){
	sprintf(temp,"%s.dbf",InfoObj->filename);
	fptr=FOPEN(temp,"w");
    }
    // opening cluster file
    if(cluster_file_flag){
	sprintf(temp,"%s.cst",InfoObj->filename);   
	clusterfptr=FOPEN(temp,"w");
    }
    
    // opening compact chain file
    if(stdout_flag!=1){
	if(!absolute_flag)
	    fprintf(fptr,"# format: ch_id f_id est_id s1 e1 s2 e2 score cr\n");
	else
	    fprintf(fptr,"# format: ch_id f_id  s1 e1 s2 e2 score cr\n");
	fprintf(fptr,"# ch_id: chain's id \n");
	fprintf(fptr,"# f_id: fragment's id \n");
	if(!absolute_flag)
	    fprintf(fptr,"# est_id: cDNA/EST's id \n");
	fprintf(fptr,"# s1: fragment's start position in cDNA/EST \n");
	fprintf(fptr,"# e1: fragment's end position in cDNA/EST \n");
	fprintf(fptr,"# s2: fragment's start position in genomic seq.  \n");
	fprintf(fptr,"# e2: fragment's end position in genomic seq.  \n");
	fprintf(fptr,"# score:  score of the chain \n");
	fprintf(fptr,"# cr: coverage ratio of the chain w.r.t. its cDNA/EST\n");
	fprintf(fptr,"# num. of genomes: %d\n",(int)dimensions+1);
	fprintf(fptr,"# total num. of fragments: %lu\n",numofblocks-1);
	//fprintf(fptr,"# max score of a chain: %lf\n",global_maxscore);
    }
    else{ 
	if(!absolute_flag)
	    printf("# format: ch_id f_id est_id s1 e1 ... sk ek score cr\n");
	else
	     printf("# format: ch_id f_id s1 e1 ... sk ek score cr\n");
	printf("# ch_id: chain's id \n");
	printf("# f_id: fragment's id \n");
	if(!absolute_flag)
	    printf("# est_id: cDNA/EST's id \n");
	printf("# s1: fragment's start position in cDNA/EST \n");
	printf("# e1: fragment's end position in cDNA/EST \n");
	printf("# s2: fragment's start position in genomic seq. \n");
	printf("# e2: fragment's end position in genomic seq. \n");
	printf("# score:  score of the chain \n");
	printf("# cr: coverage ratio of the chain w.r.t. its cDNA/EST\n");
	printf("# num. of genomes: %d\n",(int)dimensions+1);
	printf("# total num. of fragments: %lu\n",numofblocks-1);
	//printf("# max score of a chain: %lf\n",global_maxscore);
    }
    // computing the number of digits
    
    int* positiondigitarray=(int*)malloc(sizeof(long)*(dimensions+1));
    for(i=0;i<=dimensions;i++){
	positiondigitarray[i]=1+(int)log10((double)genome_sizes_array[i]);
    }
    int fragmentnumdigits=1+(int)log10((double)numofblocks);
    int contignumdigits=1+(int)log10(contigobj->number_of_contigs);
    int scorenumdigit=7+1+(int)log10(global_maxscore);
    
    Kdelem elem2;
    elem2.region=(Region*)malloc(sizeof(Region)*(dimensions+1));
    
    for(ptr=chainroots;ptr!=NULL ; ptr = ptr->next)  //last references the previous elem
    {
	// putting the clusters in file
	
	tail=NULL;
	maxscore = MIN_VALUE-10;
	if(cluster_file_flag){
	    if(ptr->child_elem->num!=0) fprintf(clusterfptr,"#cluster\n");
	    traverse_store_EST(ptr,clusterfptr,&maxscore);
	}
	else{
	    silent_traverse_cluster(ptr,&maxscore);
	}
	long contig1=contigobj->getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
					     contigobj->total_conitg_length,
					     (tail->region[dimensions].end-shiftvalue)-1);
	long contiglength;
	if(contig1==0){
	    contiglength=(long)(contigobj->recordseps[contig1]);	   
	}
	else{ 
	    contiglength=(long)(contigobj->recordseps[contig1]-contigobj->recordseps[contig1-1])-1;
	}
	
	if(check_est(chain_depth,chain_score,tail->globalscore,contigobj,percentagescore)!=1)continue;
	
	double perc_coverage=(double)tail->globalscore/(double)(contiglength*(double)multiplyingfactor);
	
	if((perc_coverage<percentagescore)&&(perc_cover_flag==1)){	   
	    continue;
	}
	
	if(stdout_flag!=1){
	    // store the fragments of the chain  in .dbf file
	    for(elem=tail;elem;elem=elem->last){
		if(elem->num==0)continue;
		// storing chain id
		fprintf(fptr,"%*lu",fragmentnumdigits,counter);
		// storing fragment id id
		fprintf(fptr,"  %*lu",fragmentnumdigits,(long)(elem-blocks));
		for(i=0;i<dimensions;i++){
		      elem2.region[i+1].start=elem->region[i].start-shiftvalue;
		      elem2.region[i+1].end=elem->region[i].end-shiftvalue;		
		}
		elem2.region[0].start=elem->region[dimensions].start-shiftvalue;
		elem2.region[0].end=elem->region[dimensions].end-shiftvalue;
		
		if (elem->num!=0){
		    if(!absolute_flag){
			long ContigId=contigobj->TransformKdelemtoRef(&elem2);
			fprintf(fptr,"  %*lu",contignumdigits,ContigId+1);		     
		    }				
		    for(i=0;i<=dimensions;i++){
			
			fprintf(fptr,"  %*lu  %*lu", positiondigitarray[i],elem2.region[i].start,
				positiondigitarray[i],elem2.region[i].end);
		    }      
		}	
		
		// storing the score
		if(perc_cover_flag==1){
		    
		    if(ptr->child_elem->num!=0)
			fprintf(fptr, "  %*.6lf  %lf\n",scorenumdigit,
				((tail->globalscore)/(double)multiplyingfactor),perc_coverage);	    
		  }
		else{
		    if(ptr->child_elem->num!=0)
			fprintf(fptr, "  %*.6lf\n",scorenumdigit,
				((tail->globalscore)/(double)multiplyingfactor));
		  }
		// end storing the score	
	    }   
	     counter++;
	}
	else{ // store in stdout
	    
	    for(elem=tail;elem;elem=elem->last){
		if(elem->num==0)continue;
		// storing chain id
		printf("%*lu",fragmentnumdigits,counter);
		// storing fragment id id
		printf("  %*lu",fragmentnumdigits,(long)(elem-blocks));
		for(i=0;i<dimensions;i++){
		    elem2.region[i+1].start=elem->region[i].start-shiftvalue;
		    elem2.region[i+1].end=elem->region[i].end-shiftvalue;		
		}
		elem2.region[0].start=elem->region[dimensions].start-shiftvalue;
		elem2.region[0].end=elem->region[dimensions].end-shiftvalue;
		
		if (elem->num!=0){
		    if(!absolute_flag){
			long ContigId=contigobj->TransformKdelemtoRef(&elem2);
			printf("  %*lu",contignumdigits,ContigId+1);		     
		    }				
		    for(i=0;i<=dimensions;i++){
			
			printf("  %*lu  %*lu", positiondigitarray[i],elem2.region[i].start,
			       positiondigitarray[i],elem2.region[i].end);		          
		    }
		}	
		
		// storing the score
		if(perc_cover_flag==1){
		    
		    if(ptr->child_elem->num!=0)
			printf("  %*.6lf  %lf\n",scorenumdigit,
			       ((tail->globalscore)/(double)multiplyingfactor),perc_coverage);	    
		}
		else{
		    if(ptr->child_elem->num!=0)
			printf("  %*.6lf\n",scorenumdigit,
			       ((tail->globalscore)/(double)multiplyingfactor));
		}
		// end storing the score
	    }
	    counter++;
	}
    }
    if(stdout_flag!=1)fclose(fptr);
    if(cluster_file_flag)fclose(clusterfptr);
    free(elem2.region);
}
//////////////////////////////////////////////////////////////////////
// Function: silent_chain_traversal
// examines the scores of local chains and computes
// some statistics from it
// the results are stored in data structure 
/////////////////////////////////////////////////////////////////////
void kdChainer::silent_chain_traversal()
{
   

    // FILE* fptr;
    // FILE* clusterfptr;
    // FILE* fptr2;
    // FILE* compactchainfptr;
    child* ptr;
    child* prev_ptr=NULL;
    long j;
    
    
    long* regionend=(long*)malloc(sizeof(long)*(dimensions+1));
    int flag=0;
    double maxscore = MIN_VALUE-10;
    global_maxscore=MIN_VALUE-10;
      
    for(ptr=chainroots;ptr!=NULL ; ptr = ptr->next)  //last references the previous elem
    {

	tail=NULL;
	maxscore = MIN_VALUE-10;
	// the following function returns a pointer
        // to the opimal chain in this cluster
	silent_traverse_cluster(ptr,&maxscore);

	// putting the optimal chain of cluster in other file
	long gap_ti=0;
	for(j=0;j<dimensions;j++){
	    gap_ti=gap_ti+genome_sizes_array[j+1]-tail->region[j].end;
	    
	}
	gap_ti=gap_ti+genome_sizes_array[0]-tail->region[dimensions].end;
	
	double tmpcheckchain;
	// num =0 means that this fragment is the ORIGIN
	if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
	    tmpcheckchain=tail->globalscore+gap_ti;
	    
	}
	else{	    	    
	    if (tail->last!=NULL) // the chain is not ending with the origin
		tmpcheckchain=tail->globalscore+gap_ti+(shiftvalue*(dimensions+1));
	    else tmpcheckchain=tail->globalscore+gap_ti;	    
	}
	
	// check chain against the options: score and chain depth
	//chain_score=max(chain_score,global_maxscore*(percentagescore/100)); 
	if(check_chain(ptr->child_elem,chain_depth,chain_score,tmpcheckchain)!=1){
	    // we can remove this chain from the chain roots
	    if(prev_ptr!=NULL){
		prev_ptr->next=ptr->next;
		free(ptr);
		ptr=prev_ptr;
	    }	    			
	    continue;       
	}
	if(ptr->child_elem->num!=0){ // the chain is not starting with the ORIGIN
	    //fprintf(fptr, "#%lf\n",(tail->globalscore+gap_ti));
	    //fprintf(compactchainfptr, "#%lf\n",tail->globalscore+gap_ti);
	    flag=0;
	}
	else{
	    if (tail->last!=NULL){ // tail->last is the ORIGIN, the shiftvalue should be considered
		//fprintf(fptr, "#%lf\n",(tail->globalscore+gap_ti+(shiftvalue*(dimensions+1))));
		//fprintf(compactchainfptr, "#%lf\n",tail->globalscore+gap_ti+(shiftvalue*(dimensions+1)));
	    }

	    flag=0;
	}
	global_maxscore=max(global_maxscore,tmpcheckchain);	
	
	prev_ptr=ptr;	
    }
    free(regionend);
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//Function: silent_traverse_cluster
// silent_traverses the subtree of every cluster to 
// the leaf (tail of chain) with the highest
// score in this cluster is stored in the tail.
//////////////////////////////////////////////////////////////////////

void kdChainer::silent_traverse_cluster(child* root, double* maxscoreptr){
	
	Kdelem_ptr elem;	
	child* child_ptr;
	
	if(root==NULL)return;
	if(root->child_elem==NULL)return;
	elem=root->child_elem;
	//tail=elem;
	// here I can process the element

	// the element will not be printed
	if(elem->globalscore > *maxscoreptr) //new best chain was found
	{
	    *maxscoreptr = elem->globalscore; 
	    tail = elem; //current last element of the chain
	}
	child_ptr=elem->start_child_list;
	while(child_ptr!=NULL){
	    silent_traverse_cluster(child_ptr,maxscoreptr);
	    child_ptr=child_ptr->next;	    
	}
	return;
}
///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////// Repeat Clustering ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

int kdChainer::repeat_cluster(double Threshold)
{

    
    //ContigSorting contigobj;
    char* contigfilename=NULL;
    Kdelem_ptr elem;
//  Kdelem_ptr temp_elem;
    double tempscore;
    double maxscore = MIN_VALUE;
    long permindex;
    long i,j;
//    long gap_ij;
//    long gap_ti;
    child* tempchild;
    
    ///////Contigs-specific Functionality
    long contig1=0;
    long contig2=0;
    long current_contig=0;
//    int contigflag=1;
//    int current_contig_flag=0;
//    double maxscoreincontig=MIN_VALUE-100;
//    double prevmaxscoreincontig=0;
//    long gapElem_ti=0;
//    char firsttimeflag=0;
    long lower_limit;
    int res;

    int Rewindflag=0;
    long* rootarray;
    long index_of_current_chain_root=0;

    /// gcfile
    long *gc_vector=new long [dimensions+1];
    for(long j=0;j<=dimensions;j++){
	//gc_vector[j]=(MyMAX_VALUE);
	gc_vector[j]=gap_threshold;
    }
    //printf("HALOOOOOOOOOOOOO %lu gc_vec %lu\n",gap_threshold,gc_vector[0]);
    
    if(gcfile!=NULL){	
	char linearr[MaxFragmentLine+1];
	char* input_token=NULL;
	
	FILE* gcfptr=FOPEN(gcfile,"r");
	if(gcfptr==NULL){
	    printf("Error: cannot open file: %s\n",gcfile);
	    delete []gc_vector;
	    return(-1);
	}
	
	fgets(linearr,MaxFragmentLine,gcfptr);	
	//printf("line:      %s\n",linearr);
	if(linearr!=NULL){
	    //  printf("%s",linearr);
	    input_token=strtok(linearr," ");
	    if(input_token!=NULL){
		gc_vector[0]=atol(input_token);
		//printf("haloooooooooooooooooooo %lu\n",gc_vector[0]);
	    }
	    for(long j=1;j<=dimensions;j++){
		input_token=strtok(NULL," ");
		if(input_token==NULL)break;
		gc_vector[j]=atol(input_token);
		//printf("haloooooooooooooooooooo %lu\n",gc_vector[j]);
	    }
	}
	fclose(gcfptr);
	
    }
    gap_threshold=gc_vector[0];
    // gcfile
    
    Kdelem_ptr lower_bound=new Kdelem;
    lower_bound->region=(Region*)malloc(sizeof(Region)*(dimensions+1));
    if(lower_bound->region==NULL)return(-1);
    
    
    contigobj=(ContigSorting*)malloc(sizeof(ContigSorting));
    res=contigobj->read_contigs(contigfilename,0);   
    if(res<0)return(-1);
     
    elem=&blocks[sorted_list[0]];
    
    elem->last=NULL;
    if(chainroots==NULL){
	chainroots=(child*)malloc(sizeof(child));
	if(chainroots==NULL)return(-1);
	chainroots->child_elem=elem;
	chainroots->next=NULL;
    }
    else{
	tempchild=(child*)malloc(sizeof(child));
	if(tempchild==NULL)return(-1);
	tempchild->next=NULL; // initialization
	tempchild->child_elem=elem;
	tempchild->next=chainroots;
	chainroots=tempchild;
	
    }
    elem->globalscore=elem->weight;
    
    kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
    lower_limit=3;

    // Array to store for each element the root of its chain
    rootarray=(long*)malloc(sizeof(long)*numofblocks);
    for(long x=0;x<numofblocks;x++){
	rootarray[x]=x;
	
    }
    blocks[0].globalscore=blocks[0].weight;

    //printf("HALOOOOOOOOOOOOO %lu gc_vec %lu\n",gap_threshold,gc_vector[0]);
    for(i = 2; i<(numofblocks*2); i++)
    {	
	elem=&blocks[sorted_list[i]];
	
	if(point_type[i]==1){  // it is start point of a fragment
	    elem->globalscore=elem->weight;
	    // get its right end
	    // add weight to its right end
	    // for gapthreshold
	    j=lower_limit;
	    while(j<i){
		if(point_type[j]!=1){
		    long delta_x=elem->region[dimensions].start-blocks[sorted_list[j]].region[dimensions].end;
		    if(delta_x>gap_threshold){
			kdtreeobj.DeactivateElem(&blocks[sorted_list[j]]);
			//flag=1;
		    }
		    else{
			break;
		    }
		    j++;
		}
		else{
		    j++;
		}				
	    }
	    lower_limit=j;
	    //end for gap threshold


	    /// We also have to check if the two fragments are in the same conitg
	    
	    contig2=contigobj->getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
					    contigobj->total_conitg_length,
					    elem->region[dimensions].end-shiftvalue-1);
	
	    
	    if(current_contig==contig2) // still in the same EST/contig
	    {
		// a modification using the region search function
		// get the lower_bound point	
				
		//lower_bound->region[dimensions].start=(long)max(0,
		//		   blocks[sorted_list[lower_limit]].region[dimensions].start-gap_threshold);
		lower_bound->region[dimensions].start=(long)max(0,elem->region[dimensions].start-gap_threshold);
		lower_bound->region[dimensions].end=lower_bound->region[dimensions].start;  		
		for(j=0;j<dimensions;j++){		    
		    lower_bound->region[j].start=(long)max(0,elem->region[j].start-gc_vector[j+1]);
		    //lower_bound->region[j].start=(long)max(0,elem->region[j].start-gap_threshold);
		    lower_bound->region[j].end=lower_bound->region[j].start;  				    
		}	
		
		permindex=kdtreeobj.region_gapped_nearestNeighbor(elem,lower_bound);					
		
		nndist=kdtreeobj.nndist;
		

		Rewindflag=0;

		if((permindex >= 0)) {
		    contig1=contigobj->getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
						    contigobj->total_conitg_length,
						    perm[permindex]->region[dimensions].start-shiftvalue+1);
		    //check for repeat constraint
		    index_of_current_chain_root=rootarray[(long)(perm[permindex]-blocks)];		    
		    if((blocks[index_of_current_chain_root].region[0].start)<=(elem->region[dimensions].end)){
			Rewindflag=1;			
		    }
		    

		}
	
		if((permindex < 0)||(perm[permindex]->num==0)||(Rewindflag==1)) //no such element
		{
		    // start a new chain
		    elem->last = NULL;
		    // the elem as a new root of a new chain
		    // and it is added to the root list
		    if(chainroots==NULL){
			chainroots=(child*)malloc(sizeof(child));
			if(chainroots==NULL)return(-1);
			chainroots->child_elem=elem;
			chainroots->next=NULL;
		    }
		    else{
			tempchild=(child*)malloc(sizeof(child));
			if(tempchild==NULL)return(-1);
			//tempchild->next=NULL; // initialization
			tempchild->child_elem=elem;
			tempchild->next=chainroots;
			chainroots=tempchild;
			
		    }
		    
		    
		    if(elem->globalscore > maxscore) //new best chain was found
		    {
			maxscore = elem->globalscore; 
			tail = elem; //current last element of the chain
		    }
		   
		}
		// here I have to check if it is in the same contig 
		else if((contig1!=contig2)){ // check for making sure
		    
		    // start a new chain
		    elem->last=NULL; 
		    // the elem as a new root of a new chain
		    // and it is added to the root list
		    if(chainroots==NULL){
			chainroots=(child*)malloc(sizeof(child));
			if(chainroots==NULL)return(-1);
			chainroots->child_elem=elem;
			chainroots->next=NULL;
		    }
		    else{
			tempchild=(child*)malloc(sizeof(child));
			if(tempchild==NULL)return(-1);
			tempchild->child_elem=elem;
			tempchild->next=chainroots;
			chainroots=tempchild;
			
		    }
		    
		    
		}
		else{
		    // both are in the same contig
		    
		    /// decide if we start a new chain
		


		    tempscore= elem->weight + nndist;		    
		    // Debug Mohamed		    
		    if((tempscore>Threshold)){
			tempscore= elem->globalscore + nndist;
			if(overlapping_value>0){ // correct tempscore and remove overlap
			    if(elem->region[dimensions].start-overlapping_value<=perm[permindex]->region[dimensions].end){
				tempscore= tempscore-(perm[permindex]->region[dimensions].end-(elem->region[dimensions].start-overlapping_value)+1);
			    }
				
			}
			elem->last = perm[permindex]; //setting chain reference	
			elem->globalscore =tempscore;
			// Adding the fragment to the child list
			if(perm[permindex]->start_child_list==NULL){
			    perm[permindex]->start_child_list=(child*)malloc(sizeof(child));
			    if(perm[permindex]->start_child_list==NULL)return(-1);
			    perm[permindex]->start_child_list->next=NULL;
			    perm[permindex]->start_child_list->child_elem=elem;
			    //perm[permindex]->start_child_list=elem;
			}
			else{
			    tempchild=(child*)malloc(sizeof(child));
			    if(tempchild==NULL)return(-1);
			    tempchild->child_elem=elem;
			    tempchild->next=perm[permindex]->start_child_list;
			    perm[permindex]->start_child_list=tempchild;			    
			}
			 // update the start point of the chain		  
			rootarray[(long)(elem-blocks)]=index_of_current_chain_root;
		    }
		    else{
			elem->last=NULL; 
			// the elem as a new root of a new chain
			// and it is added to the root list
			if(chainroots==NULL){
			    chainroots=(child*)malloc(sizeof(child));
			    if(chainroots==NULL)return(-1);
			    chainroots->child_elem=elem;
			    chainroots->next=NULL;
			}
			else{
			    tempchild=(child*)malloc(sizeof(child));
			    if(tempchild==NULL)return(-1);
			    tempchild->child_elem=elem;
			    tempchild->next=chainroots;
			    chainroots=tempchild;
			    
			}
			
		    }
		    if(elem->globalscore > maxscore) //new best chain was found
		    {
			maxscore = elem->globalscore; 
			tail = elem; //current last element of the chain
		    }		    		    
		}				
	    }
	    else{ // we are in a new contig
		
		// Deactivate the elements from lower_limit till i
		for(j=lower_limit;j<i;j++){
		    if(point_type[j]!=1)
			kdtreeobj.DeactivateElem(&blocks[sorted_list[j]]);
		}
		// set the new contig index
		lower_limit=i;
		current_contig=contig2;
		// start a new chain
		elem->last=NULL; 
		// the elem as a new root of a new chain
		// and it is added to the root list
		if(chainroots==NULL){
		    chainroots=(child*)malloc(sizeof(child));
		    if(chainroots==NULL)return(-1);
		    chainroots->child_elem=elem;
		    chainroots->next=NULL;
		}
		else{
		    tempchild=(child*)malloc(sizeof(child));
		    if(tempchild==NULL)return(-1);
		    tempchild->child_elem=elem;
		    tempchild->next=chainroots;
		    chainroots=tempchild;		    
		}
		
		if(elem->globalscore > maxscore) //new best chain was found
		{
		    maxscore = elem->globalscore; 
		    tail = elem; //current last element of the chain
		}
	    }
	    
	}
	else{ // we are at the end point of a fragment
	    // insert it in the tree
	    kdtreeobj.gapped_undelete(elem); //reinsert element into the tree
	    if(elem->globalscore > maxscore) //new best chain was found
	    {
		maxscore = elem->globalscore; 
		tail = elem; //current last element of the chain
	    }
	    
	}
    }
    //printf("HALOOOOOOOOOOOOO %lu gc_vec %lu\n",gap_threshold,gc_vector[0]);
    free(lower_bound->region);
    free(rootarray);
    delete lower_bound;
    global_maxscore=maxscore;
    delete[] gc_vector;
//    free(contigobj->recordseps);contigobj->recordseps=NULL;
//    free(contigobj);    
//    contigobj=NULL;
    return(1);

}

