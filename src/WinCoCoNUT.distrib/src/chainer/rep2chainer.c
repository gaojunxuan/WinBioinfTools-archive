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
//  double  evalue;
//  double doubleval;
//  unsigned int intval;
  FILE* fptr=NULL;
  FILE* fptr2=NULL;
  FILE* fptr3=NULL;
  FILE* logfptr=NULL;
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
    printf("Usage: rep2chainer -{filetype} filename \n");
    printf("filetype=vm if file is in vmatch format generated with the option absolute\n");
    printf("filetype=mm if file is in multimat format\n");
    printf("Output:\n ");
    printf("  -One file: 'filename.pp' in case of multimat input file \n");
    printf("  -Two files: (1) 'filename.pp'  containing direct matches.\n");  
    printf("               (2)'filename.pm' containing  palindromic matches. \n");
    return(FAIL);
  } 
  
  // checking options
  for(j=1;j<argc;j++){
    if (strcmp("-vm",argv[j])==0){
      continue;
    }
    else if(strcmp("-mm",argv[j])==0){
      continue;
    }
    else{
      strcpy(filename,argv[j]);
      filetype=1;
      fptr=FOPEN(filename,"r");
      if(fptr==NULL){
	printf("Invalid option, or file not found \n");
	return(FAIL);
      }
      fclose(fptr);
    }
  }
  if (strcmp("-vm",argv[1])==0){
    strcpy(filename,argv[2]);
    filetype=1;
    fptr=FOPEN(filename,"r");
    if(fptr==NULL){
      printf("FILE not found \n");
      return(FAIL);
    }
  }
  else if (strcmp("-mm",argv[1])==0){
    strcpy(filename,argv[2]);
    filetype=0;
    fptr=FOPEN(filename,"r");
    if(fptr==NULL){
      printf("FILE not found");
      return(FAIL);
    }
  }
  else{
    strcpy(filename,argv[1]);
     fptr=FOPEN(filename,"r");
     if(fptr==NULL){
       printf("Invalid Options or FILE not found");
       return(FAIL);
     }
     else{
       if (strcmp("-vm",argv[2])==0){	 
	 filetype=1;
       }
       if (strcmp("-mm",argv[1])==0){
	 filetype=0;
       }
       else{
	  printf("Invalid Options: filetype is not specified");
	  fclose(fptr);
	  return(FAIL);
       }
     }
  }






  if(filetype==1){
      // compute statistics
      for(k=0;k<128;k++){
	  genomesizearray[k]=0;
      }
      while(feof(fptr)==0){
	  strcpy(line,"");
	  line[0]=0;
	  fgets(line,MaxFragmentLine,fptr);
	  strcpy(tmp_line,line);
	  token = strtok (line, " " );  
	  if(token==NULL)continue;
	  if (strcmp(token,"#")==0) {
	      token = strtok ( NULL, " =");
	      
	      if((strcmp(token,"databaselength")!=0)&&(strcmp(token,"querylength")!=0)){
		  continue;
	      }
	      else{
		  //printf("%s\n",token);
		  if(strcmp(token,"databaselength")==0){
		      //printf("%s\n",token);
		      token = strtok ( NULL, " =");
		      //printf("Num %s\n",token);
		      db_length=atol(token);
		      query_length=db_length;
		      break;
		      token = strtok (NULL, " ");
		      if(token==NULL) continue;
		      if(strcmp(token,"(including")==0){
			  token = strtok ( NULL, " ");
			  //  printf("Sep %s\n",token);  
			  //db_length=db_length-atol(token);
		      }
		      break;
		  }
	      }
	  }
	  
	  
	  	      	 	  
      }
      //printf("db_len %lu query_len %lu \n",db_length,query_length);

      if((db_length==0)){
	  printf("Error, unable to locate genome lengths, use verbose mode\n");
	  fclose(fptr);
	  return 0;
      }

      while(feof(fptr)==0){
	  
	  strcpy(line,"");
	  line[0]=0;
	  fgets(line,MaxFragmentLine,fptr);
	  count=0;
	  token = strtok (line, " " );  
	  if(token==NULL)continue;
	  if (strcmp(token,"#")==0) {
	      //printf ( "The tokens first are " );
	      //printf ( "%s+", token );
	      continue;
	  }
	  do{
	      if(count==0){      
		  //fprintf(fptr2, "%d ", atoi(token));
		  length1=atoi(token);
		  maxfragmentlength1=max(maxfragmentlength1,length1);
	      }
	      token = strtok ( NULL, " ");
	      count++;
	      
	      if(count==1){
		  pos1=atoi(token);
	      }
	      
	      if(count==2){
		  if(strcmp(token,"D")==0)
		      type=1;
		  else type=0;
	      }
	      // if(count==6){  // used in case of pal
	      if(count==3){	  
		  length2=atoi(token);  
		  maxfragmentlength2=max(maxfragmentlength2,length2);
	      }
	      if(count==4){
		  
		  pos2=atoi(token);    
	      }
	      
	  }while(token);
	  if(type==1){ // Direct	      
	      genomesizearray[0]=max(genomesizearray[0],pos1+length1);
	      genomesizearray[1]=max(genomesizearray[1],pos2+length2);
	      if(open_forward_file_flag==1){
		  open_forward_file_flag=0;
		  sprintf(tempstr,"%s.pp",filename);
		  fptr2=FOPEN(tempstr,"w");
		  if(fptr2==NULL){
		      printf("unable to open file \n"); return(FAIL);
		  } 
	      
	      }
	      
	  }
	  else{ //palindrome
	      if(open_reverse_file_flag==1){
		  open_reverse_file_flag=0;
		  sprintf(tempstr,"%s.pm",filename);
		  fptr3=FOPEN(tempstr,"w"); 
		  if(fptr3==NULL){
		      printf("unable to open file \n"); return(FAIL);
		  }
	      }
	      genomesizearray[0]=max(genomesizearray[0],pos1+length1);
	      genomesizearray[1]=max(genomesizearray[1],pos2+length2);
	      //genomesizearray[0]=max(genomesizearray[0],pos1+length1);
	      //genomesizearray[1]=max(genomesizearray[1],pos2+length2);
	  }
      } // endwhile
      fclose(fptr);
      fptr=FOPEN(filename,"r");
      //
      
      sprintf(tempstr,"%s.info",filename);

      logfptr=FOPEN(tempstr,"w");       
      genomesizearray[0]=db_length;
      genomesizearray[1]=query_length;

      fprintf(logfptr,"%d %lu %lu\n",2,db_length,query_length);  
      sprintf(tempstr1,"%s.pp",filename);
      fprintf(logfptr,"%s\n",tempstr1); 
      sprintf(tempstr1,"%s.pm",filename);
      fprintf(logfptr,"%s\n",tempstr1); 
      fclose(logfptr);
      
      maxfragmentlength=max(maxfragmentlength1,maxfragmentlength2);
      offset=genomesizearray[1];
      
      
      if(fptr2!=NULL)fprintf(fptr2,">CHA %d\n",(int)2);  
      if(fptr3!=NULL)fprintf(fptr3,">CHA %d\n",(int)2);  
      while(feof(fptr)==0){
	  strcpy(line,"");
	  line[0]=0;
	  fgets(line,MaxFragmentLine,fptr);
	  count=0;
	  token = strtok (line, " " );  
	  if(token==NULL)continue;
	  if (strcmp(token,"#")==0) {
	      //printf ( "The tokens first are " );
	      //printf ( "%s+", token );
	      continue;
	  }
	  do{
	      if(count==0){      
		  //fprintf(fptr2, "%d ", atoi(token));
		  length1=atoi(token);
	      }
	      token = strtok ( NULL, " ");
	      count++;
	      //      if(count==2){ // used in case of palindroms
	      if(count==1){
		  pos1=atoi(token);
		  //fprintf(fptr2, "%d ", atoi(token));
		  //printf("%d %d %d\n", length,(pos1),(pos2)); 
	      }
	      //if(count==3){  // used in case of pal
	      if(count==2){
		  if(strcmp(token,"D")==0)
		      type=1;
		  //fprintf(fptr3, "%d\n", 0);
		  else type=0;
		  //fprintf(fptr3, "%d\n", 1);
		  //printf("%d %d %d\n", length,(pos1),(pos2)); 
	      }
	      // if(count==6){  // used in case of pal
	      if(count==3){	  
		  length2=atoi(token);  
	      }
	      if(count==4){
		  
		  pos2=atoi(token);
		  //fprintf(fptr2, "%d\n", atoi(token));
		  //printf("%d %d %d\n", length,(pos1),(pos2));    
	      }
	      
	  }while(token);
	  if(type==1){
	      //fprintf(fptr2,"%u %u %u\n",length1,pos1,pos2);
	      // the CHAINER format
	      if(fptr2!=NULL)fprintf(fptr2,"#%.2f\n",(float)length1);
	      if(fptr2!=NULL)fprintf(fptr2,"[%lu,%lu] [%lu,%lu]\n",pos1,pos1+length1-1,pos2,pos2+length2-1);
	      
	  }
	  else{
	      //fprintf(fptr3,"%u %u %u\n",length1,pos1,pos2);
	      if(fptr3!=NULL)fprintf(fptr3,"#%.2f\n",(float)length1);
	      if(fptr3!=NULL)
		  fprintf(fptr3,"[%lu,%lu] [%lu,%lu]\n",pos1,pos1+length1-1,offset-(pos2+length2-1)-1,offset-pos2-1);
	      //fprintf(fptr3,"[%d,%d] [%d,%d]\n",pos1,pos1+length1-1,pos2,pos2+length2-1);
	  }
      }
      if(fptr2!=NULL)fclose(fptr2);
      if(fptr3!=NULL)fclose(fptr3);
  }
  else{ // filetype = 0 a multimat file
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
