#include <iostream>
#include <fstream>
#include <vector>
#include <fcntl.h>
#include <string.h>
#include <dirent.h>


#include <sys/time.h>

using namespace std;
char complement(char X);


int main (int argc, char* argv[]) {
    
    
    /* default Werte setzen */
 
    
    /* keine Argumente -> Hilfe */
    if (argc == 1)
    {
	cerr << "usage: " << argv[0] << "input_fasta_file\n\n";
	cerr << " The program computes reverse complement of a sequence\n";
	cerr << " given in a multifasta file\n";
	exit(-1);
     }
   
    //base_dir_name.append((char*)dirname(argv[0]));
    
    /* Kontrolle der uebergebenen Argumente */
    char* filename=argv[1];
    string gdna="";
    string zeile;
    ifstream in (filename, ios::in);
    
    if (!in) {
	cerr << "Fehler beim Offnen des Chain File!\n";
	exit (-1);
    }
    getline(in, zeile);

    string fasta_header=zeile;

    while (getline (in, zeile)) {
	for(unsigned long i=0;i<zeile.size();i++){
	    if((zeile[i]!='\n')&&(zeile[i]!=' ')){
		char tmparr[3];
		sprintf(tmparr,"%c",zeile[i]);
		string strtmp=tmparr;
		gdna.append(strtmp);
	    }
	}
    }
    in.close();
    
    cout <<fasta_header<<"\n";
 
    unsigned long dlength=   gdna.size();
    int j=0;
    for(unsigned long i= dlength-1;i>=0;i--){
	if(j==60){
	    cout << "\n";
	    j=0;
	}
	cout << complement(gdna[i]);
	j++;
	if(i==0) break;
    }
    cout <<"\n";
    gdna.empty();
    exit(0);

}

char complement(char X){
    if((X=='A')||(X=='a')){
	return('T');
    }
    if((X=='C')||(X=='c')){
	return('G');
    }
    if((X=='G')||(X=='g')){
	return('C');
    }
    if((X=='T')||(X=='t')){
	return('A');
    }
    else{
	return(X);
    }

}
