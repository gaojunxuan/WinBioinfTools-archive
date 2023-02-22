// ContigSorting.cpp: implementation of the ContigSorting class.
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
//#include "Gcomp.h"
#include "ContigSorting.h"
#include "myglobaldef.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "string.h"



#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ContigSorting::ContigSorting()
{
	total_conitg_length=0;	
	newtotal_conitg_length=0;
	sorted_list=NULL;
	recordseps=NULL;
	chainblocks=NULL;
	perc_cover_flag=0;
	sorted_list=NULL;
	point_type=NULL;
	MIN_VALUE=(double)MyMIN_VALUE;
	recordseps=NULL;
	contig_block_array=NULL;
	newrecordseps=NULL;
	mapingcontigs=NULL;
	number_of_draft_genomes=0;
	offset_array=NULL;
	total_contig_length_array=NULL;
	number_of_contigs_array=NULL;
}

ContigSorting::~ContigSorting()
{
	if(this->chainblocks!=NULL){
		free(chainblocks);
		chainblocks=NULL;
	}
	if(sorted_list!=NULL)free(sorted_list);
	if(point_type!=NULL)free(point_type);
	if(recordseps!=NULL)free(recordseps);
	if(contig_block_array!=NULL)free(contig_block_array);
	if(newrecordseps!=NULL)free(newrecordseps);
	if(mapingcontigs!=NULL)free(mapingcontigs);
	if(offset_array!=NULL) free(offset_array);

	//offset_array= new long [number_of_draft_genomes];
	if(total_contig_length_array!=NULL)free(total_contig_length_array);
	if(number_of_contigs_array!=NULL)free(number_of_contigs_array);
	//printf("HALOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n");

}


chainblock::chainblock(){
	last=NULL;

}

chainblock::~chainblock(){
}

//////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////
// this function reads a chain file or a compact chain file
// it reports *.ctg file containing the positions transformed
// to refrence values 

long ContigSorting::kd_read_chain_file(char* filename,long in_dimensions)
{
    FILE* fptr;
    FILE* contig_fptr;
    char tempstr[300];    
    
    number_of_chains=0;
    int inputchar;
    
    long line = 0;
    
    long retval;
    char linearr[MaxFragmentLine+1];
    long max_dimensions=0;
    this->dimensions=in_dimensions;
//  float chainscore;
//  long startx, endx, starty,endy;
//  long sposition, eposition;
//  long startcontig,endcontig;

    Kdelem elem;
    elem.region=(Region*)malloc(sizeof(Region)*(dimensions+1));
    if(elem.region==NULL){
        printf(" Error: Not enough memory\n");
    }
    fptr=FOPEN(filename,"r"); 
    if(fptr==NULL){
	printf("Error: unable to open chain file. \n");
    }
    sprintf(tempstr,"%s.ctg",filename);   
    contig_fptr=FOPEN(tempstr,"w");
    while(!feof(fptr)){
	fgets(linearr,MaxFragmentLine,fptr);
	if(strcmp(linearr,"\n")==0)continue;
	if(strcmp(linearr,"")==0)continue;
	if(linearr[0]=='#')continue;
	number_of_chains++;
	strcpy(linearr,"");
	
    }
    fclose(fptr);    
    fptr=FOPEN(filename,"r");     
    if (fptr == NULL){
	printf("Unable to open chain file\n");
	return -1;
    }
    inputchar = getc (fptr);
    if (inputchar == EOF){
	printf("Error: Empty chain file\n");
	return -2;
    }
    ungetc (inputchar, fptr);    
    while (!feof(fptr))
    {
	line++;
	inputchar = getc (fptr);
	if (inputchar == EOF)
	{
	    break;
	}
	ungetc (inputchar, fptr);
	if (inputchar == '#')
	{
	    //fprintf(contig_fptr,"%c",(char)inputchar);
	    while ((inputchar = getc (fptr)) != '\n')fprintf(contig_fptr,"%c",(char)inputchar);
		 // Nothing, could process the score but skip it for now
	    fprintf(contig_fptr,"\n");
	}
	else if (inputchar == '>')
	{
	    //fprintf(contig_fptr,"%c",inputchar);	    
	    while ((inputchar = getc (fptr)) != '\n'){
		fprintf(contig_fptr,"%c",inputchar);
	    } 
	    fprintf(contig_fptr,"%c",inputchar);
                 // Nothing, could store comment but skip it for now
	}
	else
	{
	    //retval = analyzemultimatchline (fp, line-1);
	    max_dimensions=0;
	    
	    retval=analyzechainline (fptr,line,&elem,&max_dimensions);
	    if(retval==-1)break;
	    if (retval < 0)
	    {
		return -3;
	    }	    
	    // processing the chain element for the contigs by reporting reference coordinates      
	    long ContigId=TransformKdelemtoRef(&elem);
	    fprintf(contig_fptr,"%lu ",ContigId+1);
	    for(int j=0;j<max_dimensions;j++){
		fprintf(contig_fptr,"[%lu,%lu] ", elem.region[j].start,elem.region[j].end);
	    }
	    fprintf(contig_fptr,"\n");
	    // if(ret == 1) break;
	}
    }
 
 
  
    //fprintf(outfptr,"\n");
    sprintf(tempstr,"Parsing finished");
    //AfxMessageBox(tempstr);
    //printf("%s",tempstr);
    fclose(fptr);
    fclose(contig_fptr);
    free(elem.region);
    return --line;
    
    return(1);
} 
//////////////////////////////////////////////////////////////
int ContigSorting::analyzechainline (FILE *fp, unsigned int mline,Kdelem* elem,long *max_dimensions)
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
	  if (fscanf(fp,"[%lu,%lu]",&(elem->region[i].start),&(elem->region[i].end)) != 2){
	      return -1;
	  }	
	  i++;	
	  
      }
      
  }
  if(i>(*max_dimensions)){
      (*max_dimensions)=i;
  }
  
  return 0;
}

////////////////////////////////////////////////////////////
// The function TransformKdelemtoRef transforms
// an element to reference format using binary search
// over the recordseps array filled during read_contigs
// function. The reason to make this binary search instead of
// storing the reference coordinate is to reduce the space consumption

int ContigSorting::TransformKdelemtoRef(Kdelem* elem){
  //long sposition, eposition;
  
    long startcontig,endcontig;
    startcontig=getrecordnum(recordseps,number_of_contigs,total_conitg_length,elem->region[0].start+1);
    endcontig=getrecordnum(recordseps,number_of_contigs,total_conitg_length,elem->region[0].end-1);
    
    if(startcontig == 0){
	//sposition=elem->region[0].start ; // do nothing it is the first contig
    } else{
	elem->region[0].start= elem->region[0].start - recordseps[startcontig-1] - 1;
    }
    if(endcontig == 0){
	//eposition= elem->region[0].end;
    } else{
	elem->region[0].end= elem->region[0].end - recordseps[endcontig-1] - 1;
    }
    return(startcontig);
    
}

////////////////////////////////////////////////////////////
// The following function reads the contigs file and store the contigs in 
// an array in the recordseps format.
// The contig file given here is the description file format of the vmatch package
// The necessary information is in the second column of the file.
// saving flag is obsolete
int ContigSorting::read_contigs(char* filename, int savingflag)
{
	
  FILE* fptr;  
  long i=0;
  number_of_contigs=0;
  total_conitg_length=0;
//  long markpos=0;
  char linearr[MaxFragmentLine*2+1];
  char* input_token;
  long contig_length;
  
  
  //  long length;
 
  if(filename==NULL){
      recordseps=(long*)malloc((1)*(sizeof(long)));
      if(recordseps==NULL){
	  printf("Error: Not enough memory to process description file \n");
	  return(-1);
      }
      recordseps[0]=2000000000;
      number_of_contigs=1;
      total_conitg_length=recordseps[0];
      return(2); 
  }
  fptr=FOPEN(filename,"r");  
  while(!feof(fptr)){
      fgets(linearr,MaxFragmentLine*2,fptr);
      if(strcmp(linearr,"\n")==0)continue;
      if(strcmp(linearr,"")==0)continue;
      //if(linearr[0]=='#')continue;
      number_of_contigs++;
      strcpy(linearr,"");
  }
  strcpy(linearr,"");
  //fprintf(outfptr,"\n");
  fclose(fptr);
  recordseps=(long*)malloc((number_of_contigs)*(sizeof(long)));
  if(recordseps==NULL){
      printf("Error: Not enough memory to process description file \n");
      return(-1);
  }

  if(savingflag!=0){
      contig_block_array=(contigblock*)malloc((number_of_contigs)*(sizeof(contigblock)));      
      if(contig_block_array==NULL){
	  printf("Error: Not enough memory to open description file \n");
	  return(-1);
      }
      for(i=0;i<number_of_contigs;i++){
	  contig_block_array[i].flag=0;
      }
  }
  fptr=FOPEN(filename,"r"); 
  i=0;
  while((!feof(fptr))&&(i<number_of_contigs)){
	fgets(linearr,MaxFragmentLine*2,fptr);
	if(strcmp(linearr,"")==0)continue;
	if(strcmp(linearr,"\n")==0)continue;

	//skip the first column
	input_token=strtok(linearr," ");
	if(input_token==NULL)
	{	    
	    continue;
	}
	if(savingflag!=0){
	    strcpy(contig_block_array[i].contig_desc,input_token);
	}
	//number_of_chains++;
	input_token=strtok(NULL," ");
	if(i==0){
	    contig_length=atoi(input_token);
	    recordseps[i]=contig_length;
	    total_conitg_length=recordseps[i];
	}
	else{
	    contig_length=atoi(input_token);
	    recordseps[i]=recordseps[i-1]+contig_length+1;
	    total_conitg_length=total_conitg_length+contig_length+1;
	}
	i++;	
	strcpy(linearr,"");
  }
  //fprintf(outfptr,"\n");
  fclose(fptr);
  return(1);
}

///////////////////////////////////////////////////////////
long ContigSorting::getrecordnum(long *recordseps,long numofrecords,long totalwidth,long position)
{
  long *leftptr, *midptr = NULL, *rightptr, len;

  if(numofrecords == 1 || position < recordseps[0])
  {
    return 0;
  }
  if(position > recordseps[numofrecords-2])
  { 
    if(position < totalwidth)
    {
      return numofrecords - 1;
    }
//    ERROR1("cannot find position %u",position);
    return -1;
  }

 // DEBUG2(3,"getrecordnum for pos %u in [0..%u]\n",position,numofrecords-2);
  leftptr = recordseps;
  rightptr = recordseps + numofrecords - 2;
  while (leftptr<=rightptr)
  {
    len = (long) (rightptr-leftptr);
    midptr = leftptr + (len >> 1);
 //   DEBUG1(3,"left=%u,",(Uint) (leftptr - recordseps));
 //   DEBUG1(3,"mid=%u,",(Uint) (midptr - recordseps));
 //   DEBUG1(3,"right=%u,",(Uint) (rightptr - recordseps));
    if(*midptr < position)
    {
      if(position < *(midptr+1))
      {
        return (int) (midptr - recordseps + 1);
      } 
      leftptr = midptr + 1;
    } else
    {
      if(*(midptr-1) < position)
      {
        return (int) (midptr - recordseps);
      }
      rightptr = midptr-1;
    }
  }
//  ERROR1("cannot find position %u",position);
  return -1;
}

//////////////////////////////////Generalized functions for kContigs/////////////////////////////////////
////////////////////////////////////////////////////////////
// The following function reads the contigs file and store the contigs in 
// an array in the recordseps format.
// Here we might have k-draft genome
// The contig file given here is the concatenation of the description file format of the vmatch package
// each draft genome is seperated by line starting with # symbol
// The necessary information is in the second column of the file.
// saving flag is obsolete
int ContigSorting::read_kd_contigs(char* filename, int savingflag)
{
	
  FILE* fptr;  
  long i=0;
  number_of_contigs=0;
  total_conitg_length=0;
//  long markpos=0;
  char linearr[MaxFragmentLine*2+1];
  char* input_token;
  long contig_length;
  
  
  //  long length;
  
  if(filename==NULL){
      //printf("HALOOOOOOOOOOOOOOOOOOOOOOOOO \n");
      number_of_draft_genomes=1;
      offset_array=(long*) malloc(sizeof(long)*number_of_draft_genomes);
      total_contig_length_array=new long [number_of_draft_genomes];
      number_of_contigs_array=new long [number_of_draft_genomes];
      recordseps=(long*)malloc((1)*(sizeof(long)));
      if(recordseps==NULL){
	  printf("Error: Not enough memory to process description file \n");
	  return(-1);
      }
      recordseps[0]=2000000000;
      number_of_contigs=1;
      total_conitg_length=recordseps[0];
      total_contig_length_array[0]=recordseps[0];
      offset_array[0]=0;
      number_of_contigs_array[0]=1;
      return(2); 
  }
  fptr=FOPEN(filename,"r"); 
  
  while(!feof(fptr)){
      fgets(linearr,MaxFragmentLine*2,fptr);
      if(strcmp(linearr,"\n")==0)continue;
      if(strcmp(linearr,"")==0)continue;
      if(linearr[0]=='#'){
	  number_of_draft_genomes++;
	  //offset_array[number_of_draft_genomes]=number_of_contigs;
	  continue;
      }
      number_of_contigs++;
      strcpy(linearr,"");
  }
  strcpy(linearr,"");
  //fprintf(outfptr,"\n");
  fclose(fptr);

  number_of_draft_genomes++;
  offset_array= new long [number_of_draft_genomes];
  total_contig_length_array=new long [number_of_draft_genomes];
  number_of_contigs_array=new long [number_of_draft_genomes];

  for(i=0;i<number_of_draft_genomes;i++){
      total_contig_length_array[i]=0;
  }
  offset_array[0]=0;

  recordseps=(long*)malloc((number_of_contigs)*(sizeof(long)));
  if(recordseps==NULL){
      printf("Error: Not enough memory to process description file \n");
      return(-1);
  }

  
  fptr=FOPEN(filename,"r"); 
  i=0;
  long kk=0;
  long contig_id_for_offset=0;
  int new_contig_flag=0;
  long tmp_number_of_contigs=0;
  while((!feof(fptr))&&(i<number_of_contigs)){
	fgets(linearr,MaxFragmentLine*2,fptr);
	if(strcmp(linearr,"")==0)continue;
	if(strcmp(linearr,"\n")==0)continue;
	if(linearr[0]=='#'){	    	    
	    total_contig_length_array[kk]=total_conitg_length;
	    total_conitg_length=0;
	    number_of_contigs_array[kk]=tmp_number_of_contigs;
	    tmp_number_of_contigs=0;
	    kk++;
	    offset_array[kk]=contig_id_for_offset;
//	    printf("HALOOOOOOOOOCONTIGSORTI number of draft genomes %u, offset %u  \n",number_of_draft_genomes,offset_array[kk]);
	    new_contig_flag=0;
	    continue;
	}
	contig_id_for_offset++;
	tmp_number_of_contigs++;
	//skip the first column
	input_token=strtok(linearr," ");
	if(input_token==NULL)
	{	    
	    continue;
	}

	
	input_token=strtok(NULL," ");
	if(new_contig_flag==0){
	    contig_length=atol(input_token);
	    recordseps[i]=contig_length;
	    total_conitg_length=recordseps[i];
	    new_contig_flag=1;
	}
	else{
	    contig_length=atol(input_token);
	    recordseps[i]=recordseps[i-1]+contig_length+1;
	    total_conitg_length=total_conitg_length+contig_length+1;
	}
	i++;	
	strcpy(linearr,"");
  }
  //fprintf(outfptr,"\n");
  fclose(fptr);
  total_contig_length_array[kk]=total_conitg_length;
  number_of_contigs_array[kk]=tmp_number_of_contigs;
/*
  for(int d=0;d<i;d++){
      printf("recordseps %u \n", recordseps[d]);
  }
  for(int d=0;d<number_of_draft_genomes;d++){
      printf("total_contig_len_arry %u, offset_array %u, num_of_contgs %u \n", total_contig_length_array[d],offset_array[d],number_of_contigs_array[d]);
  }
*/
  return(1);
}
////////////////////////////////////////////////////////////////////////


long ContigSorting::getrecordnumGeneralized(long* positions,long* contig_ids)
{
    //long num;
    if((positions==NULL)||(contig_ids==NULL)){
	return -1;
    }
    //long startcontig,endcontig;
    //  contig_ids[0]=getrecordnum(recordseps,offset_array[0],total_conitg_length_array[0],positions[0]);
    for(long i=0;i<number_of_draft_genomes;i++){
	//num=offset_array[i]-offset_array[i-1];
	contig_ids[i]=getrecordnum(recordseps+offset_array[i],number_of_contigs_array[i],total_contig_length_array[i],positions[i]);
	
    }

//  ERROR1("cannot find position %u",position);
    return (contig_ids[0]);
}
/////////////////////////////////////////////////////////////////////////////////////
long ContigSorting::TransformKdelemtoRefGeneralized(Kdelem* elem,long* contig_id){
  //long sposition, eposition;
  
    long startcontig;
    long* positions=new long[number_of_draft_genomes];
//    long* contig_id=new long[number_of_draft_genomes];
    for(long i=0;i<number_of_draft_genomes;i++){
	positions[i]=elem->region[i].start+1;
    }
    startcontig=getrecordnumGeneralized(positions,contig_id);     
    for(long i=0;i<number_of_draft_genomes;i++){
	if(contig_id[i] == 0){
	    //sposition=elem->region[0].start ; // do nothing it is the first contig
	} else{
	    elem->region[i].start= elem->region[i].start - (recordseps+offset_array[i])[contig_id[i]-1] - 1;
	    elem->region[i].end= elem->region[i].end -(recordseps+offset_array[i])[contig_id[i]-1] - 1;
	}
    }
    delete positions;
    return(startcontig);
    
}



//////////////////////////////////////////////////////////////
long ContigSorting::getcontigsize(long position,long* recordseps){
	return(recordseps[position]);
}

//////////////////////////////////////////////////////////////////////
// General Sorting Functions
// These functions sorts the chains w.r.t. the other genome
// or w.r.t. their score
//////////////////////////////////////////////////////////////////////


/// ???? The following functions are under construction ??????????////

int ContigSorting::qsort_chains()
{
		
	long i;
	if(sorted_list==NULL){
		sorted_list=(long*)malloc(number_of_chains*sizeof(long));
		if(sorted_list==NULL){
			return(-1);
		}
		for(i=0;i<number_of_chains;i++){
			sorted_list[i]=i;
		}
	}
	else{
		return(-1);
	}

	//Qobj.quicksortkdlistkeys(kdlist,sorted_list,0,terminal_no-1);
	quicksortkdlistkeysOndim(0,sorted_list,0,(number_of_chains)-1);
	return(1);
}

//////////////////////////////////////////////////////////////////////
// this function reads a chain file or a compact chain file considering k' draft genomes
// it reports *.ctg file containing the positions transformed
// to refrence values 

long ContigSorting::kd_read_chain_file_kdraftgenome(char* filename,long in_dimensions)
{
    FILE* fptr;
    FILE* contig_fptr;
    char tempstr[300];    
    
    number_of_chains=0;
    int inputchar;
    
    long line = 0;
    
    long retval;
    char linearr[MaxFragmentLine+1];
    long max_dimensions=0;
    this->dimensions=in_dimensions;
//  float chainscore;
//  long startx, endx, starty,endy;
//  long sposition, eposition;
//  long startcontig,endcontig;

    Kdelem elem;
    elem.region=(Region*)malloc(sizeof(Region)*(dimensions+1));
    long *contig_id_array=new long[dimensions+1];
    if(elem.region==NULL){
        printf(" Error: Not enough memory\n");
    }
    fptr=FOPEN(filename,"r"); 
    if(fptr==NULL){
	printf("Error: unable to open chain file. \n");
    }
    sprintf(tempstr,"%s.ctg",filename);   
    contig_fptr=FOPEN(tempstr,"w");
    while(!feof(fptr)){
	fgets(linearr,MaxFragmentLine,fptr);
	if(strcmp(linearr,"\n")==0)continue;
	if(strcmp(linearr,"")==0)continue;
	if(linearr[0]=='#')continue;
	number_of_chains++;
	strcpy(linearr,"");
	
    }
    fclose(fptr);    
    fptr=FOPEN(filename,"r");     
    if (fptr == NULL){
	printf("Unable to open chain file\n");
	return -1;
    }
    inputchar = getc (fptr);
    if (inputchar == EOF){
	printf("Error: Empty chain file\n");
	return -2;
    }
    ungetc (inputchar, fptr);    
    while (!feof(fptr))
    {
	line++;
	inputchar = getc (fptr);
	if (inputchar == EOF)
	{
	    break;
	}
	ungetc (inputchar, fptr);
	if (inputchar == '#')
	{
	    //fprintf(contig_fptr,"%c",(char)inputchar);
	    while ((inputchar = getc (fptr)) != '\n')fprintf(contig_fptr,"%c",(char)inputchar);
		 // Nothing, could process the score but skip it for now
	    fprintf(contig_fptr,"\n");
	}
	else if (inputchar == '>')
	{
	    //fprintf(contig_fptr,"%c",inputchar);
	    
	    while ((inputchar = getc (fptr)) != '\n'){
		fprintf(contig_fptr,"%c",inputchar);
	    } 
	    fprintf(contig_fptr,"%c",inputchar);
	    //while ((inputchar = getc (fptr)) != '\n'); 
                 // Nothing, could store comment but skip it for now
	}
	else
	{
	    //retval = analyzemultimatchline (fp, line-1);
	    max_dimensions=0;
	    
	    retval=analyzechainline (fptr,line,&elem,&max_dimensions);
	    if(retval==-1)break;
	    if (retval < 0)
	    {
		return -3;
	    }	    
	    // processing the chain element for the contigs by reporting reference coordinates      
	    //long ContigId=TransformKdelemtoRef(&elem);
	    //fprintf(contig_fptr,"%lu ",ContigId+1);

	    for(int j=0;j<dimensions+1;j++){
		contig_id_array[j]=0;
	    }
	    TransformKdelemtoRefGeneralized(&elem,contig_id_array);	    	    
	    for(int j=0;j<max_dimensions;j++){
		fprintf(contig_fptr,"%lu [%lu,%lu]  ", contig_id_array[j]+1,elem.region[j].start,elem.region[j].end);
	    }
	    fprintf(contig_fptr,"\n");
	    // if(ret == 1) break;
	}
    }
 
 
  
    //fprintf(outfptr,"\n");
    sprintf(tempstr,"Parsing finished");
    //AfxMessageBox(tempstr);
    //printf("%s",tempstr);
    fclose(fptr);
    fclose(contig_fptr);
    free(elem.region);
    delete contig_id_array;
    return --line;
    
    return(1);
} 
//////////////////////////////////////////////////////////////////////
void ContigSorting::quicksortkdlistkeysOndim(unsigned char dim,long* a,long lo, long hi)
{
    //  lo ist der unterste Index, hi ist der oberste Index des
    //  zu sortierenden (Teil-)Feldes a
    long i=lo, j=hi;
	//kdMasterPoint temp;
	long temp;

    double x=this->get_dimension_value(dim,a[(long)(((double)lo+(double)hi)/2)],(long)(((double)lo+(double)hi)/2));//.x_cord;
	
    //unsigned int start,end;
	
    //  Aufteilung
    do
    {    
        while (this->get_dimension_value(dim,a[(long)i],i)<x) i++;
        while (this->get_dimension_value(dim,a[(long)j],j)>x) j--;
        
		if (i<=j)
        {
            temp=a[(unsigned int)i]; a[(unsigned int)i]=a[(unsigned int)j]; a[(unsigned int)j]=temp;
		//	tempchar=point_type[(unsigned int)i]; point_type[(unsigned int)i]=point_type[(unsigned int)j]; point_type[(unsigned int)j]=tempchar;
				i++; 
				j--;
        }
    } while (i<=j);

    // Rekursion
    if (lo<j) quicksortkdlistkeysOndim(dim,a,lo, j);
    if (i<hi) quicksortkdlistkeysOndim(dim,a,i, hi);
	
}
double ContigSorting::get_dimension_value(unsigned char dim, long index, long type_index){

	// modified to be as blocks data structure
	// Chek if the index is 	
	return(chainblocks[index].starty);
	
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Specific Sorting Functions For the one dimensional chaining
// These functions sorts the start and end points w.r.t. the other genome
//////////////////////////////////////////////////////////////////////

int ContigSorting::specific_qsort_chains()
{
		
	long i,j;
	if(sorted_list==NULL){
		sorted_list=(long*)malloc((2*number_of_chains)*sizeof(long));
		point_type=(char*)malloc(sizeof(char)*(2*number_of_chains));
		if(sorted_list==NULL){
			return(-1);
		}
		j=0;
		for(i=0;i<(number_of_chains*2);i=i+2){
			sorted_list[i]=j;
			sorted_list[i+1]=j;
			point_type[i]=1; // left end
			point_type[i+1]=2; // right end
			j++;
		}
	}
	else{
		return(-1);
	}	
	specific_quicksortkdlistkeysOndim(0,sorted_list,0,(number_of_chains*2)-1);
	return(1);
}

//////////////////////////////////////////////////////////////////////
void ContigSorting::specific_quicksortkdlistkeysOndim(unsigned char dim,long* a,long lo, long hi)
{
    //  lo ist der unterste Index, hi ist der oberste Index des
    //  zu sortierenden (Teil-)Feldes a
    long i=lo, j=hi;
	//kdMasterPoint temp;
	long temp;
	char tempchar;

    double x=this->specific_get_dimension_value(dim,a[(long)(((double)lo+(double)hi)/2)],(long)(((double)lo+(double)hi)/2));//.x_cord;
	
    //unsigned int start,end;
//	char flag=0;
    //  Aufteilung
    do
    {    
        while (this->specific_get_dimension_value(dim,a[(long)i],i)<x) i++;
        while (this->specific_get_dimension_value(dim,a[(long)j],j)>x) j--;
        
		if (i<=j)
        {
            temp=a[(unsigned int)i]; a[(unsigned int)i]=a[(unsigned int)j]; a[(unsigned int)j]=temp;
			tempchar=point_type[(unsigned int)i]; point_type[(unsigned int)i]=point_type[(unsigned int)j]; point_type[(unsigned int)j]=tempchar;
				i++; 
				j--;
        }
    } while (i<=j);

    // Rekursion
    if (lo<j) specific_quicksortkdlistkeysOndim(dim,a,lo, j);
    if (i<hi) specific_quicksortkdlistkeysOndim(dim,a,i, hi);
	
}
double ContigSorting::specific_get_dimension_value(unsigned char dim, long index, long type_index){

	// modified to be as blocks data structure
	// Chek if the index is 	
	if(point_type[type_index]==1){
		return(chainblocks[index].starty);
	}
	else{
		return(chainblocks[index].endy);
	}
	
	
}
////////////////////////////////////////////////////////////////////////////////////////
// Function: One dimensional Chaining
void ContigSorting::one_dimensional_chaining(){

	double maxscore=(double)0;
	chainblock*  maxelementptr=NULL;
	long i;
	chainblock* chainelemptr;
	
	//read_contigs("C://homes/mibrahim/projects/randomtest/BCGTub/BCG.DES",1);
	//read_chain_file("C://homes/mibrahim/projects/randomtest/BCGTub/BCGMtubH37MEM30.match.mm.sorted.ccn");	
    specific_qsort_chains();
	FILE* fptr;
	fptr=FOPEN("C://homes/mibrahim/projects/randomtest/BCGTub/OneDimChain.chn","w"); 
	
	//chainblocks[sorted_list[0]].last=0;
	for(i = 0; i<(number_of_chains*2); i++)
    {
	  chainelemptr=&chainblocks[sorted_list[i]];
      if(point_type[i]==1){ // it is start point of a chain

		  chainelemptr->last=maxelementptr;
		  //chainelemptr->lastflag=1;
		  chainelemptr->score=chainelemptr->score+maxscore;
		  fprintf(fptr,"TY %d  %lu \n", point_type[i], chainblocks[sorted_list[i]].starty);

	
	  }
	  else{ // it is an end point of a chain
		  if(maxscore<chainelemptr->score){
			  maxelementptr=chainelemptr;
			  maxscore=chainelemptr->score;
			  fprintf(fptr,"TY %d, %lu \n", point_type[i],chainblocks[sorted_list[i]].endy);

		  }
	  }
	//  fprintf(fptr,"TY %u  [%u, %u] [%u, %u] \n", 
	//	  point_type[i], 
	//	  chainblocks[sorted_list[i]].startx,chainblocks[sorted_list[i]].endx,
	//	  chainblocks[sorted_list[i]].starty,chainblocks[sorted_list[i]].endy);


	}
	one_dimensional_maximum_score=maxscore;
	one_dimensional_maxelement=maxelementptr;



	// store the 1D chains along with their contig number
	//FILE* fptr;
	//fptr=FOPEN("C://homes/mibrahim/projects/randomtest/BCGTub/OneDimChain.chn","w"); 
	chainblock* ptr=maxelementptr;
	long j=0;
	while((ptr!=NULL)&&(j<number_of_chains)){
		fprintf(fptr,"ID %lu  [%lu, %lu] [%lu, %lu] \n", ptr->contig_id,
			ptr->startx,ptr->endx,
		ptr->starty,ptr->endy);
		//i=chainblocks[sorted_list[i]].last;
		ptr=ptr->last;
		j++;
	}
	fclose(fptr);

	char tempstr[200];
	sprintf(tempstr,"max 1D score: %f", (float)one_dimensional_maximum_score);
	//AfxMessageBox(tempstr);
	printf("%s",tempstr);
	return;
}

//////////////////////////////////////////////////////////////////////
// Placement functions
//////////////////////////////////////////////////////////////////////

int ContigSorting::process_chains(char* filename)
{
  FILE* fptr;
  long i=0;


	qsort_chains(); // sorts the chians w.r.t. the second genome (the complete)
//	long referencepos;
	fptr=FOPEN("C://homes/mibrahim/projects/randomtest/BCGTub/sortedcontigs.txt","w");
	for(i=0;i<number_of_chains;i++){
		//fprintf(fptr,"ContigID: %d %s score=%d starty=%d\n",chainblocks[sorted_list[i]].contig_id,
		//	contig_block_array[chainblocks[sorted_list[i]].contig_id].contig_desc,(long)chainblocks[sorted_list[i]].score,chainblocks[sorted_list[i]].starty);
		fprintf(fptr,"ContigID: %lu  score=%lu starty=%lu\n",chainblocks[sorted_list[i]].contig_id,
			(long)chainblocks[sorted_list[i]].score,chainblocks[sorted_list[i]].starty);
	}
	fclose(fptr);

	//transform2absolute_sortedcontig();


	return(1);
}

//////////////////////////////////////////////////////////////////////
// Creating Sorted Contigs with referenced values
// Creating Sorted Contigs with absolute values
//////////////////////////////////////////////////////////////////////
int ContigSorting::transform2reference()
{
 
  long i=0;
  long ref_startx, ref_endx;
  long startcontig=0;

  //  long length;
 
  
  for(i=0;(i<number_of_chains);i++){

	startcontig=getrecordnum(recordseps,number_of_contigs,total_conitg_length,chainblocks[i].startx+1);
	//endcontig=getrecordnum(recordseps,number_of_contigs,total_conitg_length,chainblocks[i].endx-1);
	if(startcontig>0){
		ref_startx=chainblocks[i].startx-recordseps[startcontig-1]-1;
		chainblocks[i].startx=ref_startx;
		ref_endx=chainblocks[i].endx-recordseps[startcontig-1]-1;
		chainblocks[i].endx=ref_endx;

	}
	else{
	    // it is the first contig do nothing
	//	ref_startx=chainblocks[i].startx;
	//	ref_endx=chainblocks[i].endx;
	}

  }
  return(1);
}
 

//////////////////////////////////////////////////////////////////////
// Creating Sorted Contigs with referenced values
// Creating Sorted Contigs with absolute values
// The function reports only one copy of every contig
//////////////////////////////////////////////////////////////////////
int ContigSorting::transform2absolute_sortedcontig()
{
 
  long i=0;
//  long ref_startx, ref_endx;
//  long startcontig=0;
  long contig_length;
 
  transform2reference();
 
  //Creat new recordseps array
 
  newrecordseps=(long*)malloc((number_of_contigs)*(sizeof(long)));
//  mapingcontigs=(long*)malloc((number_of_contigs)*(sizeof(long)));

  for(i=0;(i<number_of_chains);i++){

	  if(contig_block_array[chainblocks[sorted_list[i]].contig_id].flag==0){ 
		  // to check if the contig placed previousely
		if(i==0){
			if(chainblocks[sorted_list[i]].contig_id!=0){
		    	contig_length=recordseps[chainblocks[sorted_list[i]].contig_id]-
				recordseps[chainblocks[sorted_list[i]].contig_id-1];
				chainblocks[sorted_list[i]].new_contig_id=i;

				
			}
			else{
				contig_length=recordseps[chainblocks[sorted_list[i]].contig_id];
			}
			newrecordseps[i]=contig_length;
			total_conitg_length=recordseps[i];
			contig_block_array[chainblocks[sorted_list[i]].contig_id].flag=1;
			chainblocks[sorted_list[i]].new_contig_id=i;
		}
		else{
			if(chainblocks[sorted_list[i]].contig_id!=0){
				contig_length=recordseps[chainblocks[sorted_list[i]].contig_id]-
				recordseps[chainblocks[sorted_list[i]].contig_id-1];		
				chainblocks[sorted_list[i]].new_contig_id=i;
			}
			else{
				contig_length=recordseps[chainblocks[sorted_list[i]].contig_id];		
				chainblocks[sorted_list[i]].new_contig_id=i;
			}
			newrecordseps[i]=newrecordseps[i-1]+contig_length+1;
			newtotal_conitg_length=newtotal_conitg_length+contig_length+1;
			contig_block_array[chainblocks[sorted_list[i]].contig_id].flag=1;
		}
	  }
	  else{
		  chainblocks[sorted_list[i]].chainflag=2;
		  newrecordseps[i]=newrecordseps[i-1];
	  }


  }

// create new absolute numbers
  for(i=0;i<number_of_contigs;i++){
		contig_block_array[i].flag=0;
  }
  
  FILE* fptr;
  fptr=FOPEN("C://homes/mibrahim/projects/randomtest/BCGTub/sortedblocks.txt","w");


  long coverage=0;
  for(i=0;(i<number_of_chains);i++){
	  if(chainblocks[sorted_list[i]].chainflag!=0){
		  continue;

	  }
	if(contig_block_array[chainblocks[sorted_list[i]].new_contig_id].flag==0){ 
		// to check if the contig placed previousely

		if(chainblocks[sorted_list[i]].new_contig_id!=0){
			chainblocks[sorted_list[i]].startx=chainblocks[sorted_list[i]].startx+
				newrecordseps[chainblocks[sorted_list[i]].new_contig_id-1];
			chainblocks[sorted_list[i]].endx=chainblocks[sorted_list[i]].endx+
				newrecordseps[chainblocks[sorted_list[i]].new_contig_id-1];
			contig_block_array[chainblocks[sorted_list[i]].new_contig_id].flag=1;
		}

		// save it on file
	coverage=coverage+(long)chainblocks[sorted_list[i]].score;
	fprintf(fptr,"#chain score %f\n",chainblocks[sorted_list[i]].score);
	fprintf(fptr,"[%lu,%lu] [%lu,%lu]\n",chainblocks[sorted_list[i]].startx,chainblocks[sorted_list[i]].endx,
		chainblocks[sorted_list[i]].starty,chainblocks[sorted_list[i]].endy);


	}
 }
  fclose(fptr);
  char tempmesschar[5000];
  //sprintf(tempmess,"no of chains: %lu, total cov %lu",number_of_chains,coverage);
  sprintf(tempmesschar,"no of chains: %lu",number_of_chains);
 // AfxMessageBox(tempmesschar);
  printf("%s",tempmesschar);
  return(1);
}
