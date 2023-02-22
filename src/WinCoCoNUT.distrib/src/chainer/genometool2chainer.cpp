#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include "myglobaldef.h"

#define max(a, b)  ((a)>=(b) ? (a) : (b))

#define FAIL 1
#define SUCCESS 0
int main(int argc,char *argv[]){
  
  char filename[1000];
  char line[MaxFragmentLine+50];
  char tmp_line[MaxFragmentLine+50];
  // char* line2;
  long length1=0;
  // long  seqno;
  long pos1=0;
  // unsigned char teampchar;
  long length2=0;
  // long seqno2;
  long pos2=0;
  int type=0; // match type
  int filetype=0;
  char  count=0;

  char RC_FLAG=0;

  FILE* fptr=NULL;
  FILE* fptr2=NULL;
  FILE* fptr3=NULL;
  FILE* logfptr=NULL;


  FILE* fptrdb=NULL;
  FILE* fptrq=NULL;

//  int i;
  int k=0;
//  double slope;
  char* token;
  char tempstr[5000];
  char tempstr1[5000];
  long  no_of_genomes=0;
  int flag=0;
//  long kk=0;
  long offset=0;
  long genomesizearray[128];
  long maxfragmentlength=0;
  long maxfragmentlength1=0;
  long maxfragmentlength2=0;
  int open_forward_file_flag=1;
  int open_reverse_file_flag=1;
  int j;
  long db_length=0;
  long query_length=0;
  if(argc!=3)
    {
    printf("Invalid options\n");
    printf("Usage: genometool2chainer -dir [p|n] -m macth_filename -db databasefile -q queryfile\n");
    printf("p=input file is forward/forward starnt\n");
    printf("n=input file is forward/reverse complement \n");
    printf("Output:\n ");
    printf("  -One file: either 'filename.pp' or  'filename.pm' depend on options\n");
    return(FAIL);
  } 
  
  // checking options
  for(j=1;j<argc;j++){
    if (strcmp("-dir",argv[j])==0){
	j++;
	if (strcmp("p",argv[j])==0){
	    RC_FLAG=0;
	}
	else if(strcmp("m",argv[j])==0){
	    RC_FLAG=1;
	}
	else{
	    printf("Error, unknown option\n");
	    break;
	}
	
      continue;
    }
    else if (strcmp("-m",argv[j])==0){
	j++;
	strcpy(filename,argv[j]);
	filetype=1;
	fptr=FOPEN(filename,"r");
	if(fptr==NULL){
	    printf("Invalid option, or file not found \n");
	    return(FAIL);
	}
	if(fptr!=NULL)fclose(fptr);	
	if(fptrdb!=NULL)fclose(fptrdb);	
	if(fptrq!=NULL)fclose(fptrq);	
	continue;
    }
    else if (strcmp("-db",argv[j])==0){
	j++;
	strcpy(filename,argv[j]);
	filetype=1;
	fptrdb=FOPEN(filename,"r");
	if(fptr==NULL){
	    printf("Invalid option, or file not found \n");
	    return(FAIL);
	}
	if(fptrdb!=NULL)fclose(fptrdb);	
	if(fptrq!=NULL)fclose(fptrq);	
	if(fptr!=NULL)fclose(fptr);	
	continue;
    }
    else if (strcmp("-q",argv[j])==0){
	j++;
	strcpy(filename,argv[j]);
	filetype=1;
	fptrq=FOPEN(filename,"r");
	if(fptr==NULL){
	    printf("Invalid option, or file not found \n");
	    return(FAIL);
	}
	if(fptrdb!=NULL)fclose(fptrdb);	
	if(fptrq!=NULL)fclose(fptrq);	
	if(fptr!=NULL)fclose(fptr);	
	continue;
    }
    
    else{
	printf("Error: unknown option: %s",argv[j]);
	exit(2);
    }
  }

  
  strcpy(filename,argv[1]);
  fptr=FOPEN(filename,"r");
  if(fptr==NULL){
      printf("Invalid Options or FILE not found");
      return(FAIL);
  }

  if(1){ 
      // filetype = 0 a multimat file
      // parsing the multimat file
      
      sprintf(tempstr,"%s.mat",filename);
      fptr2=FOPEN(tempstr,"w");
      if(fptr2==NULL){
	  printf("unable to open file \n"); return(FAIL);
      } 	  
      
      // Getting the number of genomes and priniting them
      while(feof(fptr)==0){
	  strcpy(line,"");
	  line[0]=0;
	  fgets(line,MaxFragmentLine,fptr);     
	  token = strtok (line, " " );  
	  if(token==NULL){
	      printf("Corrupted File\n");
	      return(FAIL);
	  }
	  if (strcmp(token,"#")==0) {
	      continue;
	  }
	  count=0;
	  do{
	      token=strtok ( NULL, " ");
	      if(token==NULL)break;
	      if(strcmp(token,"\n")==0)continue;
	      if(strcmp(token,"\n")==0)continue;
	      count++;
	  }while(token);
	  no_of_genomes=count;
	  break;
      }
      fprintf(fptr2,">CHA %d\n",(int)no_of_genomes);  
      fclose(fptr);
      fptr=FOPEN(filename,"r");
      if(fptr==NULL){
	  printf("FILE not found \n");
	  return(FAIL);
      }
      
      while(feof(fptr)==0){
      strcpy(line,"");
      line[0]=0;
      fgets(line,MaxFragmentLine,fptr);
      token = strtok (line, " " ); 
      if(token==NULL)break;
      if (strcmp(token,"#")==0)continue;
      if (strcmp(token,"\n")==0)continue;
      if (strcmp(token,"")==0)continue;
      length1=atoi(token);
      fprintf(fptr2,"#%.2f\n",(float)length1);
      count=0;     
      do{
	  token = strtok (NULL," " );  
	  if(token==NULL)break;
	  if (strcmp(token,"")==0) {
	      break;
	  }
	  if (strcmp(token,"#")==0) {
	  break;
	  }
	  if(strcmp(token,"\n")==0){
	      flag=1; 
	      break;
	  }
	  count++;	
	  // print the data with respect to the length
	  pos1=atoi(token);
	  fprintf(fptr2,"[%lu,%lu] ",pos1,pos1+length1-1);
      }while(token);      
      //flag=1;
      fprintf(fptr2,"\n"); 
      if(count!=no_of_genomes){
	  printf("Error in input file\n");
	  break;
      }
      }
      if(fptr2!=NULL)fclose(fptr2);
  }
  if(fptr2!=NULL)fclose(fptr);    
  return(SUCCESS);
}
