// ProcessOptions.cpp: Implementation for the ProcessOptions class.
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

#include "ProcessOptions.h"
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

ProcessOptions::ProcessOptions()
{
  optionsarray=NULL;
  fileindexarray=NULL;
}

ProcessOptions::~ProcessOptions()
{
  if(optionsarray!=NULL)free(optionsarray);
  if(fileindexarray!=NULL)free(fileindexarray);
 
}
int ProcessOptions::ParseOptions(long argc, char **argv)
{
    long i,j;
    long value;
    double fvalue;
    int operation=0;
    long filepos_inargv=0;
    FILE* fptr=NULL;
    // checking for the basic options -l, -g, -gg

    optionsarray=(double*)malloc(MaxOptions*sizeof(double));
    fileindexarray=(long*)malloc(MaxFiles*sizeof(long));
    fileindexarray[1]=0;
  
    // checking for the input file

    for(j=1;j<argc;j++){
	fptr=FOPEN(argv[j],"r");
	if((fptr!=NULL)){
	    if(strcmp(argv[j-1],"-des")==0){
		fclose(fptr);
		fptr=NULL;
		continue;
	    }
	    else if((strcmp(argv[j-1],"-s")==0)||(strcmp(argv[j-1],"-d")==0)){
		fclose(fptr);
		fptr=NULL;
		continue;
	    }
	    else if((strcmp(argv[j-1],"-coverage")==0)||
		    (strcmp(argv[j-1],"-gc")==0)||
		    (strcmp(argv[j-1],"-gcfile")==0)||
		    (strcmp(argv[j-1],"-length")==0)){
		fclose(fptr);
		fptr=NULL;
		continue;
	    }
	    else{
		//fclose(fptr);
		//fptr=NULL;
		break;
	    }
	}
    }
    if(fptr==NULL){
	printf("Invalid options: The fragment  file is not given \n");
	return(-1);
    }
    else{
	fclose(fptr);
    }
    
    //printf("%s",argv[j]);
    filepos_inargv=j;
    fileindexarray[0]=j;
    
    
    for(i=0;i<argc;i++){
	
	// processing local chaining options 
	if(strcmp(argv[i],"-l")==0){
	    // initialize options array with default values
	    optionsarray[0]=(double)0; // -s minimum score value
	    optionsarray[1]=(double)0; // -perc percentage of max. score
	    optionsarray[2]=(double)0; // -d
	    optionsarray[3]=(double)0; // -w
	    optionsarray[4]=(double)0; // -c // Obsolete Option compact chain or not
	    optionsarray[5]=(double)10; // -lw value
	    optionsarray[6]=(double)-1; // -gc gap constraints
	    optionsarray[7]=(double)0; // stdout flag
	    optionsarray[8]=(double)0; // chainer format flag
	    optionsarray[9]=(double)0; // cluster flag
	    optionsarray[10]=(double)0; // average_chain_length
	    optionsarray[11]=(double)0; // overlappign_value
	    optionsarray[12]=(double)0; // overlappign_side


	    
	    // update options array with given values
	    operation=5;
	    
	    for(j=1;j<argc;j++){
		if(strcmp(argv[j],"-s")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[0]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-length")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[10]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-perc")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if((fvalue>=(float)0)&&(fvalue<=(float)100)){
			    optionsarray[1]=(double)fvalue;
			}
			else{
			    printf("The percentage has to be between 0 and 100, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The percentage has to be between 0 and 100, the default is taken.\n");
		    }
		}
		else if(strcmp(argv[j],"-stdout")==0){		    
			optionsarray[7]=(double)1;
		
		}
		else if(strcmp(argv[j],"-chainerformat")==0){		    
			optionsarray[8]=(double)1;
		
		}
		else if(strcmp(argv[j],"-cluster")==0){		    
			optionsarray[9]=(double)1;
		
		}
		else if(strcmp(argv[j],"-gc")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>2){
			    optionsarray[6]=(double)value;
			}
			else{
			    printf("The gap constraint value should be > 2, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The gap constraint value be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-overlap")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[11]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap value has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-side")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[12]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap side has to be integer, the default is assigned. \n");
		    }
		}

		else if(strcmp(argv[j],"-d")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[2]=(double)value-1;
			}
			else{
			    printf("The depth value should be >= 1, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The depth value should be integer, the default is assigned \n");
		    }
		}
/*
  //obsolete options
		else if(strcmp(argv[j],"-w")==0){
		    optionsarray[3]=1;
		    operation=6;
		}
		else if(strcmp(argv[j],"-c")==0){
		    optionsarray[4]=1;
		}
*/
		else if(strcmp(argv[j],"-lw")==0){
		    if(getnumfloat(argv[j+1],&fvalue)){
			if(fvalue>(double)0){
			    optionsarray[5]=(double)fvalue;
			}
			else{
			    printf("The -lw should be positive, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("lw value should be float or integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-l")==0){
		    continue;
		}
		else if(filepos_inargv==j){
		    continue;
		}
		else{
		    printf("Invalid option: %s, see the manual\n",argv[j]);
		    return(-1);
		}
		
	    }
	    
	    
	    break;
	}
	// process neighboring cluster
	if(strcmp(argv[i],"-neighbor")==0){
	    // initialize options array with default values
	    optionsarray[0]=(double)0; // -s minimum score value
	    optionsarray[1]=(double)0; // -perc percentage of max. score
	    optionsarray[2]=(double)0; // -d
	    optionsarray[3]=(double)0; // -w
	    optionsarray[4]=(double)0; // -c // Obsolete Option compact chain or not
	    optionsarray[5]=(double)10; // -lw value
	    optionsarray[6]=(double)-1; // -gc gap constraints
	    optionsarray[7]=(double)0; // stdout flag
	    optionsarray[8]=(double)0; // chainer format flag
	    optionsarray[9]=(double)0; // cluster flag
	    optionsarray[10]=(double)0; // average_chain_length
	    optionsarray[11]=(double)0; // overlappign_value
	    optionsarray[12]=(double)0; // overlappign_side
	    optionsarray[13]=(double)0; // gc_file

	    
	    // update options array with given values
	    operation=10;
	    
	    for(j=1;j<argc;j++){
		if(strcmp(argv[j],"-s")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[0]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-length")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[10]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-perc")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if((fvalue>=(float)0)&&(fvalue<=(float)100)){
			    optionsarray[1]=(double)fvalue;
			}
			else{
			    printf("The percentage has to be between 0 and 100, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The percentage has to be between 0 and 100, the default is taken.\n");
		    }
		}
		else if(strcmp(argv[j],"-stdout")==0){		    
			optionsarray[7]=(double)1;
		
		}
		else if(strcmp(argv[j],"-chainerformat")==0){		    
			optionsarray[8]=(double)1;
		
		}
		else if(strcmp(argv[j],"-cluster")==0){		    
			optionsarray[9]=(double)1;
		
		}
		else if(strcmp(argv[j],"-gc")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>2){
			    optionsarray[6]=(double)value;
			}
			else{
			    printf("The gap constraint value should be > 2, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The gap constraint value be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-overlap")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[11]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap value has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-side")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[12]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap side has to be integer, the default is assigned. \n");
		    }
		}

		else if(strcmp(argv[j],"-d")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[2]=(double)value-1;
			}
			else{
			    printf("The depth value should be >= 1, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The depth value should be integer, the default is assigned \n");
		    }
		}
/*
  //obsolete options
		else if(strcmp(argv[j],"-w")==0){
		    optionsarray[3]=1;
		    operation=6;
		}
		else if(strcmp(argv[j],"-c")==0){
		    optionsarray[4]=1;
		}
*/
		else if(strcmp(argv[j],"-lw")==0){
		    if(getnumfloat(argv[j+1],&fvalue)){
			if(fvalue>(double)0){
			    optionsarray[5]=(double)fvalue;
			}
			else{
			    printf("The -lw should be positive, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("lw value should be float or integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-gcfile")==0){
		    
		    fptr=NULL;
		    fptr=FOPEN(argv[j+1],"r");
		    if(fptr!=NULL){
			optionsarray[13]=1;
			fileindexarray[2]=j+1;			
			//operation=11;
			fclose(fptr);
			j++;
		    }
		    else{
			//printf("Invalid options: No description file is given \n");
			//printf("# Warning: No description file given, the database is considered as 1 sequence. \n");
			//return(-1);
		    }
		}
		else if(strcmp(argv[j],"-neighbor")==0){
		    continue;
		}
		else if((filepos_inargv==j)||(fileindexarray[1]==j)||(fileindexarray[2]==j)){
		    continue;
		}
		
		else{
		    printf("Invalid option: %s, see the manual\n",argv[j]);
		    return(-1);
		}
		
	    }
	    
	    
	    break;
	}
	// end process neighboring cluster
	// processing global chaining options 
	if(strcmp(argv[i],"-g")==0){
	    // initialize options array with default values
	    // update options array with given values
	    operation=1;
	    optionsarray[0]=(double)0; // -w value
	    optionsarray[1]=(double)1; // -lw value
	    optionsarray[2]=(double)0; // stdout flag
	    optionsarray[3]=(double)0; // chainer format flag
	    optionsarray[4]=(double)0; // overlap value
	    optionsarray[5]=(double)0; // overlap side
	    
	    for(j=1;j<argc;j++){

	       
                if(strcmp(argv[j],"-stdout")==0){
		    optionsarray[2]=(double)1;			
		    
		}
		else if(strcmp(argv[j],"-chainerformat")==0){
		    optionsarray[3]=(double)1;			
		    
		}
		else if(strcmp(argv[j],"-lw")==0){
		    if(getnumfloat(argv[j+1],&fvalue)){
			if(fvalue>(double)0){
			    optionsarray[1]=(double)fvalue;
			}
			else{
			    printf("The -lw should be positive, the default is taken.\n");
			}
			j++;
		    } 
		    else{
			printf("lw value should be float or integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-overlap")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[4]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap value has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-side")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[5]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap side has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-g")==0){
		}
		else if(filepos_inargv==j){
		    continue;
		}
		else{
		    printf("Invalid option: %s, see the manual\n",argv[j]);
		    return(-1);
		} 
	    }
	    
	    break;
	}

	// processing global chaining with gap cost options 
	if(strcmp(argv[i],"-gg")==0){
	    // initialize options array with default values
	    // update options array with given values
	    operation=3;
	    optionsarray[0]=(double)0; // -w value
	    optionsarray[1]=(double)1; // -lw value   
	    optionsarray[2]=(double)0;
	    optionsarray[3]=(double)0; // chainerformat flag
	    optionsarray[4]=(double)0; // overlap value
	    optionsarray[5]=(double)0; // overlap side

	    for(j=1;j<argc;j++){

		if(strcmp(argv[j],"-stdout")==0){
		    optionsarray[2]=(double)1;			
		    
		}
		else if(strcmp(argv[j],"-chainerformat")==0){
		    optionsarray[3]=(double)1;			
		    
		}
		else if(strcmp(argv[j],"-lw")==0){
		    if(getnumfloat(argv[j+1],&fvalue)){
			if(fvalue>(double)0){
			    optionsarray[1]=(double)fvalue;
			    j++;
			}
			else{
			    optionsarray[1]=(double)1;
			    j++;
			    printf("-lw value has to be larger than 0, the default is taken.\n");
			}
		    }
		    else{
			printf("lw value should be float or integer, the default is assigned \n");
		    }
		} 
		else if(strcmp(argv[j],"-overlap")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[4]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap value has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-side")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[5]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap side has to be integer, the default is assigned. \n");
		    }
		}

		else if(strcmp(argv[j],"-gg")==0){
		}
		else if(filepos_inargv==j){
		    continue;
		}
		else{
		    printf("Invalid option: %s, see the manual\n",argv[j]);
		    return(-1);
		}
	    }
	    break;
	}

	// processing draft chaining options 
	if(strcmp(argv[i],"-draft")==0){
	    // initialize options array with default values
	    optionsarray[0]=(double)0;  // s minimum score value
	    optionsarray[1]=(double)0;  // stdout
	    optionsarray[2]=(double)0;  // -d
	    optionsarray[3]=(double)0;  // -w
	    optionsarray[4]=(double)0;  // -perc
	    optionsarray[5]=(double)10; // -lw value
	    optionsarray[6]=(double)-1; // -gc gap constraints
	    optionsarray[7]=(double)-1; // - description file
	    optionsarray[8]=(double)0;  // chainerformat
	    optionsarray[9]=(double)0;  // absolute positions
	    optionsarray[10]=(double)0;  // cluster file
	    optionsarray[11]=(double)0;  // chain_average_length
	    optionsarray[12]=(double)0;  // overlap value
	    optionsarray[13]=(double)0;  // overlap side
	    
	    // update options array with given values
	    operation=7;
	    
	    for(j=1;j<argc;j++){
		if(strcmp(argv[j],"-s")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[0]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-length")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[11]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-stdout")==0){
		 
			optionsarray[1]=(double)1;
		}
		else if(strcmp(argv[j],"-chainerformat")==0){
		 
			optionsarray[8]=(double)1;
		}
		else if(strcmp(argv[j],"-absolute")==0){		    
		    optionsarray[9]=(double)1;
		    
		}
		else if(strcmp(argv[j],"-cluster")==0){		    
		    optionsarray[10]=(double)1;
		    
		}
		else if(strcmp(argv[j],"-gc")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>2){
			    optionsarray[6]=(double)value;
			}
			else{
			    printf("The gap constraint value should be > 2, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The gap constraint value be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-overlap")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[12]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap value has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-side")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[13]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap side has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-d")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[2]=(double)value-1;
			}
			else{
			    printf("The depth value should be >= 1, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The depth value should be integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-w")==0){
		    optionsarray[3]=1;
		    //operation=6;
		}
		else if(strcmp(argv[j],"-des")==0){
		    fptr=NULL;
		    fptr=FOPEN(argv[j+1],"r");
		    if(fptr!=NULL){
			optionsarray[7]=1;
			fileindexarray[1]=j+1;
			operation=7;
			fclose(fptr);
			j++;
		    }
		    else{
			//printf("Invalid options: No description file is given \n");
			printf("# Warning: No description file given, the database is considered as 1 sequence. \n");
			//return(-1);
		    }
		}
		else if(strcmp(argv[j],"-perc")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if((fvalue>=(float)0)&&(fvalue<=(float)100)){
			    optionsarray[4]=(double)fvalue;
			}
			else{
			    printf("The percentage has to be between 0 and 100, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The percentage has to be between 0 and 100, the default is taken.\n");
		    }
		}
		else if(strcmp(argv[j],"-c")==0){
		    optionsarray[4]=1;
		}
		else if(strcmp(argv[j],"-lw")==0){
		    if(getnumfloat(argv[j+1],&fvalue)){
			if(fvalue>(double)0){
			    optionsarray[5]=(double)fvalue;
			}
			else{
			    printf("The -lw should be positive, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("lw value should be float or integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-draft")==0){
		    continue;
		}
		else if((filepos_inargv==j)||(fileindexarray[1]==j)){
		    continue;
		}
		else{
		    printf("Invalid option: %s, see the manual\n",argv[j]);
		    return(-1);
		}
		
	    }
	    if( optionsarray[7]!=1){
		//printf("Invalid options: No description file is given\n");
		printf("# warning: No description file is given, the database is considered as 1 sequence. \n");
		//return(-1);
	    }
	    
	    break;
	}
	
	// Processing Options for EST
	
	if(strcmp(argv[i],"-est")==0){
	    // initialize options array with default values
	    optionsarray[0]=(double)0; // s minimum score value
	    optionsarray[1]=(double)0; // stdout flag
	    optionsarray[2]=(double)0; // -d
	    optionsarray[3]=(double)0; // -w
	    optionsarray[4]=(double)0; // -c // Obsolute Option compact chain or not
	    optionsarray[5]=(double)1; // -lw value
	    optionsarray[6]=(double)-1; // perecentage value
	    optionsarray[7]=(double)-1; // description file
	    optionsarray[8]=(double)-1; // gap constraint
	    optionsarray[9]=(double)0; // chainer format
	    optionsarray[10]=(double)0; // absolute
	    optionsarray[11]=(double)0; // cluster
	    optionsarray[12]=(double)0; // overlap value
	    optionsarray[13]=(double)0; // overlap side

	    // update options array with given values
	    operation=8;
	    
	    for(j=1;j<argc;j++){
		if(strcmp(argv[j],"-s")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[0]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-stdout")==0){
		   
			optionsarray[1]=(double)1;
		}
		else if(strcmp(argv[j],"-chainerformat")==0){		   
			optionsarray[9]=(double)1;
		}
		else if(strcmp(argv[j],"-absolute")==0){		    
		    optionsarray[10]=(double)1;
		    
		}
		else if(strcmp(argv[j],"-cluster")==0){		    
		    optionsarray[11]=(double)1;
		    
		}
		else if(strcmp(argv[j],"-coverage")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if((fvalue>=(float)0)&&(fvalue<=(float)1)){
			    optionsarray[6]=(double)fvalue;
			}
			else{
			    printf("The coverage ratio has to be between 0 and 1, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The coverage ratio has to be between 0 and 1, the default is taken.\n");
		    }
		}
		else if(strcmp(argv[j],"-gc")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>2){
			    optionsarray[8]=(double)value;
			}
			else{
			    printf("The gap constraint value should be > 2, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The gap constraint value be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-overlap")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[12]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap value has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-side")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[13]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap side  has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-d")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[2]=(double)value-1;
			    
			}
			else{
			    printf("The depth value should be >=1, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The depth value should be integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-w")==0){
		    //optionsarray[3]=1;
		    //operation=6;
		}
		else if(strcmp(argv[j],"-des")==0){
		    fptr=NULL;
		    fptr=FOPEN(argv[j+1],"r");
		    if(fptr!=NULL){
			optionsarray[7]=1;
			fileindexarray[1]=j+1;
			operation=8;
			fclose(fptr);
			j++;
		    }
		    else{
			//printf("Invalid options: No description file is given \n");
			//return(-1);
                        //modification to accept a database of 1 sequence.
			printf("# warning: No description file is given, the database is considered as 1 sequence. \n");
		    }
		}
		else if(strcmp(argv[j],"-c")==0){
		    optionsarray[4]=1;
		}
		else if(strcmp(argv[j],"-lw")==0){
		    if(getnumfloat(argv[j+1],&fvalue)){
			if(fvalue>(double)0){
			    optionsarray[5]=(double)fvalue;
			}
			else{
			    printf("The -lw should be positive, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("lw value should be float or integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-est")==0){
		    continue;
		}
		else if((filepos_inargv==j)||(fileindexarray[1]==j)){
		    continue;
		}
		else{
		    printf("Invalid option: %s, see the manual\n",argv[j]);
		    return(-1);
		}
		
	    }
	    if( optionsarray[7]!=1){
		//modification to accept a database of 1 sequence.
		printf("# warning: No description file is given, the database is considered as 1 sequence. \n");
		//printf("Invalid option: No description file is given\n");
		//return(-1);
	    }
	    
	    break;
	}		
	///End of Processing Options for EST

	// start processing repeats
	if(strcmp(argv[i],"-r")==0){
	    // initialize options array with default values
	    optionsarray[0]=(double)0; // -s minimum score value
	    optionsarray[1]=(double)0; // -perc percentage of max. score
	    optionsarray[2]=(double)0; // -d
	    optionsarray[3]=(double)0; // -w
	    optionsarray[4]=(double)0; // -c // Obsolete Option compact chain or not
	    optionsarray[5]=(double)10; // -lw value
	    optionsarray[6]=(double)-1; // -gc gap constraints
	    optionsarray[7]=(double)0; // stdout flag
	    optionsarray[8]=(double)0; // chainer format flag
	    optionsarray[9]=(double)0; // cluster flag
	    optionsarray[10]=(double)0; // average_chain_length
	    optionsarray[11]=(double)0; // overlappign_value
	    optionsarray[12]=(double)0; // overlappign_side
	    optionsarray[13]=(double)0; // gc_file
	    

	    
	    // update options array with given values
	    operation=9;
	    
	    for(j=1;j<argc;j++){
		if(strcmp(argv[j],"-s")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[0]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-length")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[10]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-perc")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if((fvalue>=(float)0)&&(fvalue<=(float)100)){
			    optionsarray[1]=(double)fvalue;
			}
			else{
			    printf("The percentage has to be between 0 and 100, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The percentage has to be between 0 and 100, the default is taken.\n");
		    }
		}
		
		else if(strcmp(argv[j],"-stdout")==0){		    
			optionsarray[7]=(double)1;
		
		}
		else if(strcmp(argv[j],"-chainerformat")==0){		    
			optionsarray[8]=(double)1;
		
		}
		else if(strcmp(argv[j],"-cluster")==0){		    
			optionsarray[9]=(double)1;
		
		}
		else if(strcmp(argv[j],"-gc")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>2){
			    optionsarray[6]=(double)value;
			}
			else{
			    printf("The gap constraint value should be > 2, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The gap constraint value be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-overlap")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[11]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap value has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-side")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[12]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap side has to be integer, the default is assigned. \n");
		    }
		}

		else if(strcmp(argv[j],"-d")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[2]=(double)value-1;
			}
			else{
			    printf("The depth value should be >= 1, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The depth value should be integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-gcfile")==0){
		    fptr=NULL;
		    //printf("HALOOOOOOOOO %s\n",argv[j]);
		    
		    fptr=FOPEN(argv[j+1],"r");
		    if(fptr!=NULL){
			//optionsarray[13]=1;
			fileindexarray[2]=j+1;
			//printf("HALOOOOOOOOO %s",argv[j+1]);
			operation=7;
			fclose(fptr);
			j++;
		    }
		    else{
			//printf("Invalid options: No description file is given \n");
			printf("# Warning: No gap file given, no limit will be considered \n");
			//return(-1);
		    }
		}
		else if(strcmp(argv[j],"-lw")==0){
		    if(getnumfloat(argv[j+1],&fvalue)){
			if(fvalue>(double)0){
			    optionsarray[5]=(double)fvalue;
			}
			else{
			    printf("The -lw should be positive, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("lw value should be float or integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-r")==0){
		    continue;
		}
		else if((filepos_inargv==j)||(fileindexarray[1]==j)||(fileindexarray[2]==j)){
		    continue;
		}
		else{
		    printf("Invalid option: %s, see the manual\n",argv[j]);
		    return(-1);
		}
		
	    }
	    
	    
	    break;
	}
	// end processing repeats
	// start repeat cluster 
	if(strcmp(argv[i],"-rc")==0){
	    // initialize options array with default values
	    optionsarray[0]=(double)0; // -s minimum score value
	    optionsarray[1]=(double)0; // -perc percentage of max. score
	    optionsarray[2]=(double)0; // -d
	    optionsarray[3]=(double)0; // -w
	    optionsarray[4]=(double)0; // -c // Obsolete Option compact chain or not
	    optionsarray[5]=(double)10; // -lw value
	    optionsarray[6]=(double)-1; // -gc gap constraints
	    optionsarray[7]=(double)0; // stdout flag
	    optionsarray[8]=(double)0; // chainer format flag
	    optionsarray[9]=(double)0; // cluster flag
	    optionsarray[10]=(double)0; // average_chain_length
	    optionsarray[11]=(double)0; // overlappign_value
	    optionsarray[12]=(double)0; // overlappign_side
	    optionsarray[13]=(double)0; // gc_file
	    

	    
	    // update options array with given values
	    operation=11;
	    
	    for(j=1;j<argc;j++){
		
		if(strcmp(argv[j],"-s")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[0]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-length")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if(fvalue>(double)0){ 
			    optionsarray[10]=(double)fvalue;
			}
			else{
			    printf("s value should be positive. The default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("s value should be float or integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-perc")==0){
		    if(getnumfloat(argv[j+1],&fvalue)==1){
			if((fvalue>=(float)0)&&(fvalue<=(float)100)){
			    optionsarray[1]=(double)fvalue;
			}
			else{
			    printf("The percentage has to be between 0 and 100, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The percentage has to be between 0 and 100, the default is taken.\n");
		    }
		}
		
		else if(strcmp(argv[j],"-stdout")==0){		    
			optionsarray[7]=(double)1;
		
		}
		else if(strcmp(argv[j],"-chainerformat")==0){		    
			optionsarray[8]=(double)1;
		
		}
		else if(strcmp(argv[j],"-cluster")==0){		    
			optionsarray[9]=(double)1;
		
		}
		else if(strcmp(argv[j],"-gc")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>2){
			    optionsarray[6]=(double)value;
			}
			else{
			    printf("The gap constraint value should be > 2, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The gap constraint value be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-overlap")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[11]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap value has to be integer, the default is assigned. \n");
		    }
		}
		else if(strcmp(argv[j],"-side")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[12]=(double)value;
			}
			j++;
		    }
		    else{
			j++;
			printf("The overlap side has to be integer, the default is assigned. \n");
		    }
		}

		else if(strcmp(argv[j],"-d")==0){
		    if(getnumint(argv[j+1],&value)==1){
			if(value>0){
			    optionsarray[2]=(double)value-1;
			}
			else{
			    printf("The depth value should be >= 1, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("The depth value should be integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-gcfile")==0){
		    
		    fptr=NULL;
		    fptr=FOPEN(argv[j+1],"r");
		    if(fptr!=NULL){
			optionsarray[13]=1;
			fileindexarray[2]=j+1;			
			//operation=11;
			fclose(fptr);
			j++;
		    }
		    else{
			//printf("Invalid options: No description file is given \n");
			//printf("# Warning: No description file given, the database is considered as 1 sequence. \n");
			//return(-1);
		    }
		}
		else if(strcmp(argv[j],"-lw")==0){
		    if(getnumfloat(argv[j+1],&fvalue)){
			if(fvalue>(double)0){
			    optionsarray[5]=(double)fvalue;
			}
			else{
			    printf("The -lw should be positive, the default is taken.\n");
			}
			j++;
		    }
		    else{
			j++;
			printf("lw value should be float or integer, the default is assigned \n");
		    }
		}
		else if(strcmp(argv[j],"-rc")==0){
		    continue;
		}
		else if((filepos_inargv==j)||(fileindexarray[1]==j)||(fileindexarray[2]==j)){
		    continue;
		}
		else{
		    printf("Invalid option: %s, see the manual\n",argv[j]);
		    return(-1);
		}
		
	    }
	    
	    
	    break;
	}
	// end repeat cluster

	//else{
	    //printf("Invalid options: 'l', 'g' or 'gg' not specified");
	//}
	
    } // end for loop 
    return(operation);
}
////////////////////////////////////////////////////////////////
// Function getnumfloat checks if the given input is numerical 
//float value or not. The value is stored in fvalue

int ProcessOptions::getnumfloat(char* token, double* fvalue)
{
  if(cmpdecimal(token,1)==1){
    *fvalue=(double)atof(token);
    return (1);
  }
  else{
    return(0);
  }

}
////////////////////////////////////////////////////////////////
// Function getnumint checks if the given input is  
// integer value or not. The value is stored in value

int ProcessOptions::getnumint(char* token, long* value){
   if(cmpdecimal(token,0)==1){
     *value=atoi(token);
     return (1);
   }
   else{
     return(0);
   }
}

////////////////////////////////////////////////////////////////
// Function cmpdecimal checks if the given input is 
// integer or float. returns 1 on success. 

int ProcessOptions::cmpdecimal(char* token, int type){
  
  long length=strlen(token);
  for(long i=0;i<length;i++){
    if(token[i]=='1'){
      continue;
    }
    else if(token[i]=='2'){
      continue;
    }    
    else if(token[i]=='3'){
      continue;
    }
    else if(token[i]=='4'){
      continue;
    }
    else if(token[i]=='5'){
      continue;
    }
    else if(token[i]=='6'){
      continue;
    }
    else if(token[i]=='7'){
      continue;
    }
    else if(token[i]=='8'){
      continue;
    }
    else if(token[i]=='9'){
      continue;
    }
    else if(token[i]=='0'){
      continue;
    }
    else{
      if(type==1){
	if(token[i]=='.'){
	  continue;
	}
	else{
	  return(-1);
	}
      }
      else{
	return(-1);
      }
    }
  }
  return(1);
}
