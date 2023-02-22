#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include "myglobaldef.h"

#define max(a, b)  ((a)>=(b) ? (a) : (b))

#define FAIL 1
#define SUCCESS 0
int main(int argc,char *argv[]){
  
  char* filename;
  char* filename2;
  char line[MaxFragmentLine+50];
  long length1=0;
  // long  seqno;
  long pos1=0;
  // unsigned char teampchar;
  long length2=0;

  char  count=0;

  char RC_FLAG=0;

  FILE* fptr=NULL;
  FILE* fptr2=NULL;
  FILE* logfptr=NULL;


  int k=0;
  char* token;
  char tempstr[5000];
  char tempstr1[5000];
  long  no_of_genomes=0;
  int flag=0;

  long offset=0;
  long maxfragmentlength=0;
  long maxfragmentlength1=0;
  long maxfragmentlength2=0;
  int j;
  long db_length=0;
  long query_length=0;
  if(argc!=5){
    printf("Invalid options\n");
    printf("Usage: gt2chainer -i macth_filename -o outfilename\n");
    printf("Example:\n ");
    printf(" >gt2chainer -i match_file\n");
    return(FAIL);
  } 
  
  // checking options
  for(j=1;j<argc;j++){
    if (strcmp("-p",argv[j])==0){
	j++;
	RC_FLAG=0;
      continue;
    }
    else if (strcmp("-m",argv[j])==0){
	j++;
	RC_FLAG=1;
      continue;
    }
    else if (strcmp("-i",argv[j])==0){
	j++;
	//strcpy(filename,argv[j]);
	filename=argv[j];
	fptr=FOPEN(filename,"r");
	if(fptr==NULL){
	    printf("Invalid option, or file not found \n");
	    return(FAIL);
	}
	if(fptr!=NULL)fclose(fptr);	
	continue;
    }

    else if (strcmp("-o",argv[j])==0){
	j++;
	//strcpy(filename2,argv[j]);
	filename2=argv[j];
	fptr2=FOPEN(filename2,"w");
	if(fptr2==NULL){
	    printf("Invalid option, or file not found \n");
	    return(FAIL);
	}
	if(fptr2!=NULL)fclose(fptr2);	
	continue;
    }
    
    else{
	printf("Error: unknown option: %s",argv[j]);
	exit(2);
    }
  }

  


  if(1){ 
      // filetype = 0 a multimat file
      // parsing the multimat file
   
      no_of_genomes=2;
    
      fptr=FOPEN(filename,"r");
      if(fptr==NULL){
	  printf("FILE not found \n");
	  return(FAIL);
      }

      fptr2=FOPEN(filename2,"w");
      if(fptr2==NULL){
	  printf("FILE not found \n");
	  return(FAIL);
      }
      fprintf(fptr2,">CHA %d\n",(int)no_of_genomes);  


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
	    
	      // print the data with respect to the length
	      if(count==0){
		  pos1=atoi(token);
		  fprintf(fptr2,"[%lu,%lu] ",pos1,pos1+length1-1);
	      }
	      if(count==2){
		  pos1=atoi(token);
		  fprintf(fptr2,"[%lu,%lu] ",pos1,pos1+length1-1);
	      }
	      
	      count++;	
	  }while(token);      
	  //flag=1;
	  fprintf(fptr2,"\n"); 
	  if(count-1!=no_of_genomes){
	      printf("Error in input file\n");
	      break;
	  }
      }
      if(fptr2!=NULL)fclose(fptr2);
  }
  if(fptr2!=NULL)fclose(fptr);    
  return(SUCCESS);
}
