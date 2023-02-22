//############################################################################
//Parallel sequence Alignment
//It needs MPICH2 to be installed before. 
//This is Done By: Hisham Adel Hasan
//Research Assistant, Nile University,Cairo,Egypt.
//Bioinformatics group
//E-mail:hisham.mohamed@nileu.edu.eg
//      :hosham2004@yahoo.com
//Date:15 Dec 2008
//############################################################################
#include <iostream>
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include "mpi.h"
#include<string>
#include "PSA.h"
using namespace std;
using std::string;

int main(int argc,char* argv[]) 
{
    int rank;
    int size;
    int provided;
    PSA SA;
    // MPI_Init(&argc,&argv);
    MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0)
    {
	SA.openfiles(size);
	//cout<<"Rank =0 Number of diagonals = "<<SA.NumOfDiag<<endl;
	//cout<<"Rank =0 Each Node will take = "<<SA.DiaPN<<endl;
        //cout<<"Rank =0 Last Node Will take = "<<SA.DiaLN<<endl;
    }
    if(rank!=0)
    {

	SA.Recv(0,rank,size);
	//cout<<"Rank =0 Number of diagonals = "<<SA.NumOfDiag<<endl;
	//cout<<"Diag for Node "<<rank<<" = "<<SA.DiaPN<<endl;
	//cout<<"Diag for Node "<<rank<<" = "<<SA.DiaLN<<endl;
      
    }
    SA.Broadcast( rank);

//    cout<<"File 1 size "<<SA.file1s<<"I am Rank "<<rank<<endl;
    //  cout<<"File 2 size "<<SA.file2s<<"I am Rank "<<rank<<endl;

    SA.Allocate(rank,size);
MPI_Barrier(MPI_COMM_WORLD);
SA.Calc(rank,size);


    SA.gatherdata(rank,size);
/*
MPI_Barrier(MPI_COMM_WORLD);
    if(rank==2)
    {
	int xx;
	xx=SA.getposition(5,0,0);
	    cout<<"XX = "<<xx<<endl;
    }
*/

    MPI_Finalize();

}//end main
