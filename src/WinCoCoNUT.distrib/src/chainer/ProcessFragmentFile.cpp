// ProcessFragmentFile.cpp: implementation of the ProcessFragmentFile class.
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/*
  Copyright by Mohamed I. Abouelhoda (C) 2004
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


#include "ProcessFragmentFile.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>


#include "myglobaldef.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ProcessFragmentFile::ProcessFragmentFile()
{
	no_of_fragments=0;
	no_of_genomes=0;
	max_genome_size=0;
	max_fragment_length=0;
	min_fragment_length=10000000;
	max_fragment_weight=0;
	filename[0]=0;

	for (int i=0;i<MaxGenomes;i++){
		genome_size_array[i]=0;
		max_locations_array[i]=0;
	}
	printf("*//////////////////////////////////////////////////*\n");
	printf("****************************************************\n");
	printf("********                                    ********\n");
	printf("********              CHAINER               ********\n");
	printf("********                                    ********\n");
	printf("****************************************************\n");
	printf("*//////////////////////////////////////////////////*\n");
	printf("*          This is CHAINER version 3.0             *\n");
	printf("*                 %d-bits version                  *\n",VERSION);
	printf("*                                                  *\n");
 	printf("*   Copyright by Mohamed I. Abouelhoda  (C) 2006   *\n");
	printf("*                                                  *\n");
	printf("*  Unauthorized commercial usage and distribution  *\n");
        printf("*         of this program is prohibited.           *\n");
	printf("*        Contact the author for a license.         *\n");
	printf("*      Please report bugs and suggestions to       *\n");
	printf("*           <mohamed.ibrahim@uni-ulm.de>           *\n");      
	printf("*//////////////////////////////////////////////////*\n");
	printf("****************************************************\n");	  
	
}

ProcessFragmentFile::~ProcessFragmentFile()
{

}

//////////////////////////////////////////////////////////////////
//The function opens the fragment file and its log file in the CHAINER Format
//The function reads the accompanied statistics file if it exists or
//The function  creates it.
/////////////////////////////////////////////////////////////////
int ProcessFragmentFile::GetFragmentFileStatisticsFormat2(char *infilename)
{
    FILE* fptr;
    Kdelem elem;
    char linearr[MaxFragmentLine+1];
    char* input_token;
    long i=0;
    
    int retval;
    
    char logfilename[MaxFragmentLine];
    char statisticsfilename[MaxFragmentLine];
    
    // getting number of matches
    fptr=FOPEN(infilename,"r");
    if(fptr==NULL){
	printf("# non existing fragments file\n");
	return(-1);
    }
    fclose(fptr);
	
   
    strcpy(filename,infilename);
    sprintf(statisticsfilename,"%s.stc",filename);
    //fptr=FOPEN(statisticsfilename,"r");
    fptr=NULL;
    if(fptr==NULL){		
	fptr=FOPEN(filename,"r");
	int tmpflag=0;
	while(!feof(fptr)){
	    if(tmpflag==0){
		tmpflag=1;
		retval=fscanf(fptr,">CHA %lu\n",&no_of_genomes);
		if(retval<=0){
		    printf("# non-correct CHAINER file format\n");
		    return(-1);
		}
		elem.region=(Region*)malloc(sizeof(Region)*(no_of_genomes));
		continue;		
	    }
	    fgets(linearr,MaxFragmentLine,fptr);
	    if(strcmp(linearr,"\n")==0)continue;
	    if(strcmp(linearr,"")==0)continue;	
	    if(linearr[0]=='#')continue;
	    if(linearr[0]==' ')continue;
	    if(strlen(linearr)==0)continue;
	    
	    retval=analyzefragmentline (linearr,&elem);	
	    
	    if(retval<0){		
		break;
	    }
	    
	    no_of_fragments++;
	  
	    strcpy(linearr,"");	    
	}// end while
	
	fclose(fptr);
	free(elem.region);

	fptr=NULL;
	printf("# creating statistics file, it will be created\n");
	fptr=FOPEN(statisticsfilename,"w");
	if(fptr==NULL){
	    printf("Not enough memory for creating statistics file\n");
	    exit(-1);
	}
	fprintf(fptr,"%lu %lu %lu %lu \n",no_of_genomes,no_of_fragments,max_fragment_length,min_fragment_length);
	for(i=0;i<no_of_genomes;i++)fprintf(fptr,"%lu ",max_locations_array[i]);
	fprintf(fptr,"\n");
	fclose(fptr);	
    }
    else{ // load only the information in the statistics file	
	//
    }
  
    fptr=NULL;
    sprintf(logfilename,"%s.log",filename);
    fptr=FOPEN(logfilename,"r");
    if(fptr==NULL){
	for(i=0;i<no_of_genomes;i++){	
	    //genome_size_array[i]=0;
	    genome_size_array[i]=max_locations_array[i]+100;
	    max_genome_size=MAX(max_genome_size,genome_size_array[i]);	
	}
	//fclose(fptr);
	return(1);
    }
    linearr[0]=0;
    fgets(linearr,MaxFragmentLine,fptr);
    input_token=strtok(linearr," ");
    if(input_token!=NULL){		
	genome_size_array[0]=(long)atof(input_token)+100;
	max_genome_size=genome_size_array[0];
	for (i=1;i<MaxGenomes;i++){
	    input_token=strtok(NULL," ");
	    if(input_token==NULL){
		break;
	    }
	    genome_size_array[i]=(long)atof(input_token)+100;
	    max_genome_size=MAX(max_genome_size,genome_size_array[i]);
	}
    }
    fclose(fptr);
    // in case the genome sizes are not given
    if(genome_size_array[0]==0){
	for(i=0;i<no_of_genomes;i++){	
	    genome_size_array[i]=max_locations_array[i]+100;
	    max_genome_size=MAX(max_genome_size,genome_size_array[i]);	
	}
    }
    
    return(1); 
}


int ProcessFragmentFile::analyzefragmentline (char* linearr,Kdelem* elem)
{

  
 long i = 0; 
//  long length, position;
 int inputchar;
 int retval=0;
 int position=0; 
 char tmparr[1000];
 while (1)
 {
      // bail out at end of line
    inputchar = linearr[position];
    if(inputchar==EOF){
	//ungetc (inputchar, fp);
      return(-1);
    }
    if (inputchar == '\n')
    {
		position++;
		break;
    }
    
    if (inputchar == '[')
    {
		//ungetc (inputchar, fp);
		retval=sscanf(linearr+position,"[%lu,%lu]",&(elem->region[i].start),&(elem->region[i].end));
		sprintf(tmparr,"[%lu,%lu]",elem->region[i].start,elem->region[i].end);
		position=position+strlen(tmparr);
		if ( retval!= 2)
		{
			return -1;
		}
		else{
		  	max_fragment_length=MAX(max_fragment_length,(elem->region[i].end-elem->region[i].start+1));			  
			min_fragment_length=MIN(min_fragment_length,(elem->region[i].end-elem->region[i].start+1));
			max_locations_array[i]=MAX(max_locations_array[i],elem->region[i].end);
		
		}	
		i++;	

    }
    else{
	break;
    }
	position++;
    
  }
 
// if(i>(*max_dimensions)){
//	(*max_dimensions)=i;
// }

  return 0;
}

////////////////////////////////////////////////////////////
