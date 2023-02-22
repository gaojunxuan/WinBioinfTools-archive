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
#include "mpi.h"
class PSA
{
private:
void countt(int key,int rank,int size);
bool found;
bool found2;
int position;
int position2;
int position3;
int position4;
int *arr;
int Max;
int* askp;
void getmax(int arr[],int& max);
void MaxSave(int array[],int pos);
void Value(int r,int c,int match,int missmatch,int gab,int posi,int temp1,int temp2,int temp3);
int STARTD;
int ENDD;
int STARTR;
int STARTC;
int relvD;
int Stap;

long double Glocation2;



public:
int posi;
 int* resl;
int savearr;
    char* file1;
    char* file2;
    int file1s;
    int file2s;
    int NumOfDiag;
    int DiaPN;
    int DiaLN;
    int r;
    int c;
    int DCalc;
    MPI_Status status;
    PSA();
    void openfiles(int size);
    void Broadcast(int rank);
    void Recv(int Case,int rank,int size);
    void Allocate(int rank,int size);
    void Calc(int rank,int size);
    void gatherdata(int rank,int size);
//void *exchange (void* rptr);


int getposition(int r,int c,int position);
int reNP(int r,int c,int key);


}; 
 
