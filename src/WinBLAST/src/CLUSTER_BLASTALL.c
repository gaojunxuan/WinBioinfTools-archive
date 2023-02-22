//############################################################################
//BLASTALL.c
//It needs MPICH2 to be installed before. 
//This is Done By: Hisham Adel Hasan
//Research Assistant, Nile University,Cairo,Egypt.
//Bioinformatics group
//E-mail:hisham.mohamed@nileu.edu.eg
//      :hosham2004@yahoo.com
//Date:15 Dec 2008
//############################################################################
#define _CRT_SECURE_NO_DEPRECATE
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
void masterwork();
void slavework();
void getresults();
void sendresults();
//*****************************************************************************
 int rank,size;        //rank of each node  size-Number of nodes
        MPI_Status status;
        MPI_Request request;
        FILE *infile;
        FILE *outfile;  
        int nvf=0;
	char temp;
	char *b;
	char *sub;
	int i;
	FILE *BlastOut;
	int count=0;
	char *Blarr;
	int Rec;
	char filename[100];
	int error;
	int Ready=0;
        double startwtimeBo= 0.0, endwtimeBo;
         double optime=0.0;
        double optime1=0.0;
        double totaltime3=0.0;
        double startCtime= 0.0, endCtime;
        struct stat file_status;
	int Filesize;
	char BlastOrder[200];
	char filePath[200];
	char querryPath[200];
        int selectt=0;
		int hh;
		int inca;
		double max;
		int kkkk;
		double alltend;
		int x;
		char name[100];
//  clock_t allt;
double allt=0.0;
double OpCo=0.0;

//*****************************************************************************
main (int argc,char* argv[])
{
double *timeee;

        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm_size(MPI_COMM_WORLD,&size);
timeee=(double *)malloc(size* sizeof(double));
	if(rank==0)
	{
            masterwork();
printf("\nSelect:\n\n");                //to get the size of the data base
printf("\n\n-Press 1 if DataBase is Protein and Query is Protien:");
printf("\n\n-Press 2 if DataBase is Nucleotide and Query is Nucleotide:");
printf("\n\n-Press 3 if DataBase is Protein and Query is Nucleotide:");
printf("\n\n-Press 4 if DataBase is Nucleotide and Query is Protein:");
printf("\n\nEnter number:"); 
fflush (stdout);
scanf("%d",&selectt);

	}//end if rank ==0

	MPI_Bcast(&nvf,1,MPI_INT,0,MPI_COMM_WORLD);  //BroadCast the querry size to all nodes
        MPI_Bcast(&selectt,1,MPI_INT,0,MPI_COMM_WORLD);

        if(rank!=0)
       {
         b=(char *)malloc(nvf * sizeof(char));  //memory allocation for the querry array
       }

        MPI_Bcast(b,nvf,MPI_CHAR,0,MPI_COMM_WORLD);  //BroadCast the querry to all nodes
        MPI_Bcast(&Filesize,1,MPI_INT,0,MPI_COMM_WORLD);
       if(rank!=0)
        {
         slavework();

        }//rank!=0



MPI_Barrier(MPI_COMM_WORLD);
//Apply BLASTALL
if (rank==0)
{
if(selectt==1)
{
	sprintf(BlastOrder,"D:\\Blastout\\bin\\blastall -p blastp -d C:\\BLASTDB\\data.nt -i %s -z %d -o C:\\Blast-result\\OutMaster.txt",querryPath,Filesize);

}
else if (selectt==2)
{

sprintf(BlastOrder,"D:\\Blastout\\bin\\blastall -p blastn -d C:\\BLASTDB\\data.nt -i %s -z %d -o C:\\Blast-result\\OutMaster.txt",querryPath,Filesize);

}
else if(selectt==3)
{

sprintf(BlastOrder,"D:\\Blastout\\bin\\blastall -p blastx -d C:\\BLASTDB\\data.nt -i %s -z %d -o C:\\Blast-result\\OutMaster.txt",querryPath,Filesize);
}
else if(selectt==4)
{
sprintf(BlastOrder,"D:\\Blastout\\bin\\blastall -p tblastn -d C:\\BLASTDB\\data.nt -i %s -z %d -o C:\\Blast-result\\OutMaster.txt",querryPath,Filesize);

}
       //time start before start search
    //sprintf(BlastOrder,"blastall -p blastn -d ./BLASTDB/data.nt -i %s -z %d -o ./Blast-result/OutMaster.txt",querryPath,Filesize);
//startwtimeBo=endwtimeBo=0;
startwtimeBo=endwtimeBo=0.0;
startwtimeBo= MPI_Wtime(); 
//clock_t startop = clock();
    error=system(BlastOrder);
//clock_t endop = clock();

//optime=(double)(endop-startop);
   
 endwtimeBo= MPI_Wtime();
 timeee[0]=endwtimeBo-startwtimeBo;
//timeee[0]=optime/CLOCKS_PER_SEC;
//printf("\n%f\n",timeee[0]);
    fprintf(stderr, "%s\n", strerror(errno)); 
//totaltime1=endwtimeBo-startwtimeBo;
//error=system("blastall -p blastn -d ./BLASTDB/data.nt -i querry.txt -z 3164247 -o /home/hisham/Blast-result/OutMaster.txt");

if(error!=0)
{

printf("\nCannot apply Search on the Node 0....\n");
fflush (stdout);
}

inca=0;
//double timeee[size];
for(hh=1;hh<size;hh++)                 //get operation time
{
MPI_Recv(&optime1,1,MPI_DOUBLE,hh,321,MPI_COMM_WORLD,&status);
//optime1=optime1/CLOCKS_PER_SEC;
timeee[hh]=optime1;

//totaltime1=totaltime1+totaltime2;
//totaltime2=0;
}


max=timeee[0];
kkkk=0;
while(kkkk<size)
{

if(timeee[kkkk]>max)
{
max=timeee[kkkk];

}
kkkk++;
}

printf("\nMax value %f\n",max);
fflush (stdout);

//Start receving Querry Result

//clock_t alltend = clock();

alltend=MPI_Wtime();

//OpCo=alltend-allt;
getresults();
//totaltime2=0;
/*
for(hh=1;hh<size;hh++)                 //get comm time
{
MPI_Recv(&totaltime3,1,MPI_DOUBLE,hh,432,MPI_COMM_WORLD,&status);
totaltime2=totaltime2+totaltime3;
//totaltime3=0;
}
*/
//printf("Time elapsed for querry: %f\n", (totaltime1 / CLOCKS_PER_SEC));
printf("\nTime for query = %f\n\n", max);
fflush (stdout);
//printf("\nCommunication time =%f\n\n",OpCo-max);
//printf("\nTotal time = %f \n\n",OpCo);

 } //end if rank =0



//*****************************************************
//Apply Blastall then send result to Master Node
if (rank !=0)
{  

if(selectt==1)
{
	sprintf(BlastOrder,"D:\\Blastout\\bin\\blastall -p blastp -d C:\\BLASTDB\\data.nt -i C:\\querry.txt -z %d -o C:\\out.txt",Filesize);
//fflush (stdout);
}
else if (selectt==2)
{
	sprintf(BlastOrder,"D:\\Blastout\\bin\\blastall -p blastn -d C:\\BLASTDB\\data.nt -i C:\\querry.txt -z %d -o C:\\out.txt",Filesize);
//fflush (stdout);
}
else if(selectt==3)
{

	sprintf(BlastOrder,"D:\\Blastout\\bin\\blastall -p blastx -d C:\\BLASTDB\\data.nt -i C:\\querry.txt -z %d -o C:\\out.txt",Filesize);
//fflush (stdout);
}
else if(selectt==4)
{
	sprintf(BlastOrder,"D:\\Blastout\\bin\\blastall -p tblastn -d C:\\BLASTDB\\data.nt -i C:\\querry.txt -z %d -o C:\\out.txt",Filesize);
//fflush (stdout);
}

    //sprintf(BlastOrder,"./blast-2.2.18/bin/blastall -p blastn -d ./db/data.nt -i querry.txt -z %d -o out.txt",Filesize);
//startwtimeBo=endwtimeBo=0;
 //startwtimeBo=endwtimeBo=0.0;
 startwtimeBo= MPI_Wtime();
//clock_t startop = clock();
    error=system(BlastOrder);
	//printf("\nError: %s ", strerror(errno));
//printf("order= %s\n, My rank=%d",BlastOrder,rank);

//fflush (stdout);
//clock_t endop = clock();
endwtimeBo= MPI_Wtime();
//optime=endop-startop;
  //
optime=endwtimeBo-startwtimeBo;
   fprintf(stderr, "%s\n", strerror(errno)); 
//totaltime1=endwtimeBo-startwtimeBo;

if(error!=0)
{
	MPI_Get_processor_name(name,&x);
printf("\nCannot apply Search on the Node %d ....My name=%s\n",rank,name);

		
fflush (stdout);
exit(8);

}
if (error==0)
{
	MPI_Get_processor_name(name,&x);
printf("\nSearch done on the Node %d ....My name=%s\n",rank,name);

}		
fflush (stdout);
MPI_Send(&optime,1,MPI_DOUBLE,0,321,MPI_COMM_WORLD);
//error=system("./blast-2.2.18/bin/blastall -p blastn -d ./db/data.nt -i querry.txt -z 3164247 -o out.txt");


sendresults();


}//if rank!=0

//***************************************************************


MPI_Finalize();
}//end main

//*****************MasterWork()*******************************************
void masterwork()
{

	system("mkdir C:\\Blast-result");
	    printf("\n\n");
		fflush (stdout);
	printf("\nPlease Enter the Path of the Database:");                //to get the size of the data base
	fflush (stdout);
	scanf("%s",filePath);
	printf("\n\n");
	fflush (stdout);
	printf("\nPlease Enter the Path of the querry File:");             //path of querry file
	fflush (stdout);
	scanf("%s",querryPath);
printf("\n%s",querryPath);  
fflush (stdout);
          if(stat(filePath, &file_status) != 0)                            //get size of database file
          {
            perror("\nCould not find the file\n");
			fflush (stdout);
            exit(8);
          }
              Filesize=file_status.st_size;

               infile=fopen(querryPath,"r");   //Start to read the querry
         if(infile==NULL)
	{  
            printf("\nError:Cannot Open The File.....\n");
			fflush (stdout);
            exit(8);

	}//end if infile=null

while(!feof(infile))                              //count number of elements in the file
{
fscanf(infile,"%c",&temp);    //read element by

nvf=nvf+1;  

}//end while
fclose(infile);

infile=fopen(querryPath,"r");

b=(char *)malloc(nvf * sizeof(char));  //make an array b[] to save querry in it 

for(i=0;i<nvf;++i)
{
  fscanf(infile,"%c",&b[i]);  //start reading elements from file and put them in b array

}//end for
fclose(infile);

}//end master Work
//********************************************************************
//*****************Slave work()****************************************
void slavework()
{
	outfile=fopen("C:\\querry.txt","w");      //save the querry to a file
  
 for(i=0;i<nvf-1;i++)  //To save data from 0 to nvf
    {
	fprintf(outfile,"%c",b[i]);
	fflush (stdout);

    }//end for loop

fclose(outfile);

}
//*****************end SlaveWork()*************************************
//******************getresults()****************************************
void getresults()
{

int numofnodes=1;
int j=-1;
int Rerr;
//*****************************************************************************************************
//The Following code REcv data from Nodes
//When the node finish it send to the master node it's rank then send to the master node that it is ready
//then start sending data to master node
//
//******************************************************************************************************
while(numofnodes!=size)        //still in while loop until numofnodes=size
{
    Rerr= MPI_Recv(&j,1,MPI_INT,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);  //get the rank of the node which finished 

 if(Rerr==MPI_SUCCESS)
 {
     //  MPI_Recv(&Ready,1,MPI_INT,j,25,MPI_COMM_WORLD,&status);                 //tell the master node that it is ready and it will send data now
     
     // if(Ready==123 &j!=-1)
if(j!=-1)
       {
	   MPI_Recv(&Rec,1,MPI_INT,j,90,MPI_COMM_WORLD,&status);            //Recv size of data

	   Blarr=(char *)malloc(Rec * sizeof(char));                       //memory allocation for the data

	   MPI_Recv(Blarr,Rec,MPI_CHAR,j,98,MPI_COMM_WORLD,&status);       //Recv the data in Blarr array

	   sprintf(filename,"C:\\Blast-result\\OutSlave%d.txt",j);

          BlastOut=fopen(filename,"w");

           for(i=0;i<Rec-1;i++)  //To save data in file
            {
	       fprintf(BlastOut,"%c",Blarr[i]);

            }//end for loop

 fclose(BlastOut);  //close the file

 // Ready=0;

       }
       numofnodes++;    //increase num of nodes

 }
  

}//end while

//***************sequentional code******************************************
//************************************************************************
// the following code Recv from node 1 then 2 then 3 then .........
//which is not efficient
//
//************************************************************************
/*
int j;
for(j=1;j<size;j++)
{
    MPI_Recv(&Rec,1,MPI_INT,j,90,MPI_COMM_WORLD,&status);  
Blarr=(char *)malloc(Rec * sizeof(char));

//printf("\n'm hereee master %d\n",Rec);
MPI_Recv(Blarr,Rec,MPI_CHAR,j,98,MPI_COMM_WORLD,&status);

sprintf(filename,"/home/hisham/Blast-result/OutSlave%d.txt",j);

BlastOut=fopen(filename,"w");

 for(i=0;i<Rec-1;i++)  //To save data from 0 to nvalues+inc
    {
	fprintf(BlastOut,"%c",Blarr[i]);
//	printf("%c",revbuff[i]);

    }//end for loop
 fclose(BlastOut);



}//end for loop
*/
//*********************************************************************

}
//***********************end getresults()******************************************
//***********************Start sendresults()****************************************
void sendresults()
{
	BlastOut=fopen("C:\\out.txt","r");                  //open result file

while(!feof(BlastOut))                           //determine number of elements in the file
{
// printf("\n'm node %d\n",rank);   
fscanf(BlastOut,"%c",&temp);   

count=count+1;  

}//end while

fclose(BlastOut);
BlastOut=fopen("C:\\out.txt","r");

Blarr=(char *)malloc(count * sizeof(char));
for(i=0;i<count;++i)
{

  fscanf(BlastOut,"%c",&Blarr[i]);  //start reading elements from file and put them in b array

}//end for
fclose(BlastOut);

//Ready=123;


error=MPI_Send(&rank,1,MPI_INT,0,20,MPI_COMM_WORLD);                   //send rank of the node to tell master node that 'm ready
if(error==MPI_SUCCESS)
{
  
    //  MPI_Send(&Ready,1,MPI_INT,0,25,MPI_COMM_WORLD);        //tell master node that it is ready
//startCtime=MPI_Wtime();
//clock_t start = clock();
    MPI_Send(&count,1,MPI_INT,0,90,MPI_COMM_WORLD);       //send size of the array
    MPI_Send(Blarr,count,MPI_CHAR,0,98,MPI_COMM_WORLD);   //send the array of data
//totaltime2=(double)clock() - start;
//endCtime=MPI_Wtime();
}//end if error



//totaltime2=endCtime-startCtime;
//MPI_Send(&totaltime2,1,MPI_DOUBLE,0,432,MPI_COMM_WORLD); 

}
//**********************end sendresults()**********************************************
