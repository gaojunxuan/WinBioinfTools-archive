#ifdef MACHINE_64
#define FOPEN fopen64
#define VERSION 64
#else
#define FOPEN fopen
#define VERSION 32
#endif

#define MaxFragmentLine 4999      // Max. length of a line of any file parsed. 
#define MyMIN_VALUE -2000000000   // Min. value of score
#define MyMAX_VALUE 2000000000    // Max value of score
#define MaxFileNameSize 1499      // Max. length of a file name  
#define MaxGenomes 128            // Max number of genomes read
