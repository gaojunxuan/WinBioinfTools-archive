#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/mman.h>
#include <ctype.h>

#include <fcntl.h>

#include "align_chain.hpp"
#include "fragment.hpp"
#include "math.h"

#include "memmap.hpp"

#include "alignment.hpp"
#include "myglobaldef.h"

using namespace std;

#define TMPFILETEMPLATE "_XXXXXX"
#define MAXFILENAMELENGTH (strlen (TMPFILETEMPLATE) + 5)
#define MAXFILELINELENGTH 1024
#define FILESUFFIXFASTA ".fna"
#define FILESUFFIXGDE   ".gde"
#define FILESUFFIXDND   ".dnd"
#define SEPARATOR       UCHAR_MAX         // separator symbol in multiple seq

#define MY_MAX 2000000000
#define MY_MIN -1

#define MAXCOMMANDLINELENGTH 512

#ifndef CLUSTALWNAME
#define CLUSTALWNAME "clustalw"
#endif

#define CHECKZERO(INT)\
  if (INT != 0)\
  {\
    return INT;\
  }


/*EE

  These parameters should not be changed as the function
  \texttt{computemultialignment} depends on them.

*/
#define CLUSTALWBASICPARAMETERS "-batch -align -type=DNA"
#define CLUSTALWOUTPUTPARAMETERS "-output=GDE -outorder=INPUT"
//#define CLUSTALWOUTPUTPARAMETERS "-output=CLUSTAL  -outorder=INPUT"

/*EE

  These parameters can be changed to finetune \textit{CLUSTAL\,W}.  A
  simple example with one parameter is \texttt{"-dnamatrix=CLUSTALW"},
  more parameters are possible.

*/
#ifndef CLUSTALWTUNINGPARAMETERS
#define CLUSTALWTUNINGPARAMETERS "-dnamatrix=IUB"
#endif


int id_swap = 1;

//extern int luecke;
//---------------------------------------------
align_chain::align_chain (){

    
    total_fragment_lengths_in_genomes=NULL;       	
    total_aligned_gap_lengths_in_genomes=NULL;   
    total_unaligned_gap_lengths_in_genomes=NULL;     
    chain_boundaris=NULL;
  
}


/* Konstrukteur: wird gerufen, wenn Headerzeile gelesen wurde */
align_chain::align_chain (long in_number_of_genomes,int* in_genome_files_handles,long* in_start_offset_array) {
	/* Id steht in den Fragmentzeilen */
	id = 0;
	/* Name wird erst ermittelt, wenn die Sequenz eingelesen wird */
	name = "";

//	threshold_alignable_gap= 2000000000;
	threshold_alignable_gap= 1000;
	number_of_genomes=in_number_of_genomes;
	chain_boundaris=new interval[number_of_genomes];
	for(long i=0;i<number_of_genomes;i++){	    
	    chain_boundaris[i].start=INT_MAX;
	    chain_boundaris[i].end=-1;
	    
	}
	genome_files_handles=in_genome_files_handles;
	start_offset_array=in_start_offset_array; // to skip header line of fasta file
	boundary_threshold=2*getpagesize();
	//boundary_threshold=2000000000;
	
	// For statistics
	unaligned_gap_counter=0;
	total_fragment_lengths_in_genomes=new long [number_of_genomes];
	total_aligned_gap_lengths_in_genomes=new long [number_of_genomes];
	total_unaligned_gap_lengths_in_genomes=new long [number_of_genomes];
	for(long i=0;i<number_of_genomes;i++){
	    total_fragment_lengths_in_genomes[i]=0;
	    total_aligned_gap_lengths_in_genomes[i]=0;
	    total_unaligned_gap_lengths_in_genomes[i]=0;
	}
	seq_no_digit=(long)log10((double)number_of_genomes);
	total_identity=0;
	output_mode=0;
	fragment_output_mode=0;
	genome_size_array=NULL;
	//direction="";
	//fasta_line_size=70;
	contig_offset_array=NULL;
	palindrome_flag=0;
	relative_flag=0;
	
}


/* Destrukteur */
align_chain::~align_chain () {
    
    for (unsigned long i = 0; i < fragmente.size (); i++){	
	delete[] fragmente[i]->positions;
	delete fragmente[i];
    }	
    fragmente.clear ();
    if(total_fragment_lengths_in_genomes!=NULL){
	delete[] total_fragment_lengths_in_genomes;total_fragment_lengths_in_genomes=NULL;
    }
    if(total_aligned_gap_lengths_in_genomes!=NULL){
	delete[] total_aligned_gap_lengths_in_genomes;
	total_aligned_gap_lengths_in_genomes=NULL;
    }
    if(total_unaligned_gap_lengths_in_genomes!=NULL){
	delete[] total_unaligned_gap_lengths_in_genomes;
	total_unaligned_gap_lengths_in_genomes=NULL;
	
    }
    if(chain_boundaris!=NULL){
	delete[] chain_boundaris;
	chain_boundaris=NULL;
    }
    
}

void align_chain::set_contig_offset_array(int in_relative_flag,long* in_contig_offset_array,long* in_contig_id_array,contigs* in_all_contigs_array,int in_palindrome_flag){
    
    palindrome_flag=in_palindrome_flag;
    relative_flag=in_relative_flag;
    contig_offset_array=in_contig_offset_array;
    contig_id_array=in_contig_id_array;
    all_contigs_array=in_all_contigs_array;
}


/* Parst die eingelesene Zeile und splittet sie in Id, c_tsart, c_ende,
	g_start und g_ende auf
	line   [in]   string   aufzuteilende Zeichenkette
	sep [in]   char*    Menge von Trennzeichen
	words  [out]  vector<string> Vektor mit den getrennten Substrings */


void align_chain::fsplit (string line, const char* sep, vector<string> & words) {
		
	string::size_type a = 0, e;
	while ((a = line.find_first_not_of (sep, a)) != string::npos) {
		e = line.find_first_of (sep, a);
		if (e != string::npos) {
			words.push_back (line.substr (a, e-a));
			a = e + 1;
		}
		else {
			words.push_back (line.substr (a));
			break;
		}
	}
}


int align_chain::generatefilenames (char *fasta, char *gde, char *dnd){
/*
  char my_template[MAXFILENAMELENGTH];
  char *prefix;

  sprintf (my_template, "%s", TMPFILETEMPLATE);
  prefix = mktemp (my_template);
  if (prefix == NULL)
  {
      cerr << "Error: cannot generate temporary filename for template" << my_template << "\n";
    return -1;
  }

  fasta = strcpy (fasta, prefix);
  fasta = strcat (fasta, FILESUFFIXFASTA);

  gde = strcpy (gde, prefix);
  gde = strcat (gde, FILESUFFIXGDE);

  dnd = strcpy (dnd, prefix);
  dnd = strcat (dnd, FILESUFFIXDND);
  return 0;
*/
    return 0;
}


/* Fragment hinzufuegen: wird aufgerufen, wenn ein Fragmentzeile eingelesen
   wurde */
void align_chain::addFragment (string frag) {
	vector<string> words;
	std::string::size_type x = 0;
	
	/* zuerst Leerzeichen durch ',' ersetzen: ',' sind Trenner */
	while ((long)(x = frag.find (" ", 0)) != -1) {
		frag.replace (x, 1, ",");
	}

	x = 0;
	/* Dann Klammern ('[', ']') entfernen */
	while ((long)(x = frag.find_first_of ("][", 0)) != -1)
		frag.replace (x, 1, "");

	/* Dann kann gesplittet werden */
	
	fsplit (frag, ",", words);
	
	long w_size=words.size();
	//number_of_genomes=w_size/2; // will be modified so that the number of genomes is given in the prog. input
	//id = atoi (words[0].c_str ());
	// debug mohamed
	
	long frag_id = atoi(words[w_size-1].c_str());
	
	interval* positions= new interval [number_of_genomes];		

	long j=0;
	for(long i=0;i<number_of_genomes;i++){
	    positions[i].start=atoi(words[j].c_str());
	    positions[i].end=atoi(words[j+1].c_str());
	    j=j+2;
	}
	fragmente.push_back(new fragment(positions,frag_id));		
	
	// here set the chain boundaries
	for(long i=0;i<number_of_genomes;i++){	    	    
	    if(chain_boundaris[i].start>positions[i].start){
		chain_boundaris[i].start=positions[i].start;
	    }
	    if(chain_boundaris[i].end<positions[i].end){
		chain_boundaris[i].end=positions[i].end;
	    }
	    
	}	    
	words.clear();
}




/* ermittelt die gesamte Sequenz der align_chain oder die entsprechende gDNA-
	Subsequenz. Die Sequenz wird dabei aus dem entsprechenden cDNA File
	oder gDNA Chromosom File mit Memory Mapping gelesen */
//-------------------------------- Neu Get Sequ --------------------------------//
void align_chain::GetSeqFromFastaFile (string* in_gdna,long start, long ende, int file_handle, long in_start_offset, long fasta_line_size) {
	
        *in_gdna="";
        
        /* wird zurueckgegeben */
	//string seq;
	/* Seitengroesse im RAM */
	long pagesize = getpagesize ();
	
	/* gDNA Sequenz holen */
	
	/* so lang ist die Subsequenz */
	long seq_len = ende - start + 1;

/*--------- dieser Teil muss so sein, weil es sonst nicht funktioniert---*/

	/* keine Ahnung warum, aber das muss so sein */
	ende += in_start_offset;

	
	
	/* so viel Newlinezeichen sind in der Sequenz vor der
	   gesuchten Subsequenz */
	long nl = start / fasta_line_size;
	/* Startpunkt um Anzahl NewLines verschieben */
	start += nl;
	/* gleiches gilt fuer das Ende */
	nl = ende / fasta_line_size;
	ende += nl;
	
/*--- ich hab aber keine Ahnung warum das so ist -----*/

	/* soviel muss vorne weg */
	long offset_page = start % pagesize;
	/* solang ist die Sequenz im Arbeitsspeicher */
	long len = ende - start + 1 + offset_page;
	/* offset im File muss Vielfaches von pagesize sein */
	long offset_seq = pagesize * (start / pagesize);
	/* enthaelt den Inhalt, des in den speicher kopierten
	   Filesbereichs */
	long maplen = ((len / pagesize) + 1) * pagesize;
	
	
	/* Memory Mapping:
	   data: Adresse, an der der Mappingbereich im Speicher beginnt;
	   0: System bestimmt wo im Speicher der mapping Auschnitt steht;
	   maplen: soviel Speicher wird reserviert;
	   PROT_READ, MAP_SHARD: flags fuer Datei-/Speicherzugriff ;
	   fd_gdna: Filedescriptor (zu mappendes File);
	   offset_seq: nach sovielen Bytes im File wird mit dem Mappen
	   begonnen (muss Vielfaches von pagesize sein);
	*/
	char* data = (char*) mmap ((caddr_t)0, maplen, PROT_READ, MAP_SHARED, file_handle, offset_seq);
	
	/* der Startpunkt der Subsequenz muss ebenfalls um das Start
	   Offset verschoben werden */
	long i = offset_page + in_start_offset;
	
	long char_num=0;
	
	
	/* Sequenz in String speichern */
	//cout << "Sequence length in data " << maplen <<  ", i =" << i << " \n";
	//long my_seq_size=seq.size ();
	//long my_seq_size=in_gdna->size();
	
	

	while (in_gdna->size() < (unsigned long) seq_len) {
	    /* keine Newlines mitspeichern */
	    if (data[i] != '\n'){
		//seq.append (1, data[i]);
		in_gdna->append (1, data[i]);
		char_num++;
		//cout << data[i];
	    }
	    i++;
	}
	
	
	/* Speicherbereich freigeben. Achtung: data ist die Adresse,
	   die bei mmap zurueckgegeben wird! */
	munmap ((caddr_t)data, maplen);
	
//	string tmp_ss;
//	return tmp_ss;
	//return seq;
}

void align_chain::GetSeqFromFastaFileReverseComplement(string* in_gdna,long in_start, long in_ende, long genome_size,int file_handle, long in_start_offset,long fasta_line_size) {
	
        *in_gdna="";
	// adjust coordinates
	long start,ende;
	start=genome_size-in_ende-1;
	ende=genome_size-in_start-1;
	
	if(ende<start){
	    cerr << start << " " << ende << " Orignal: "<< in_start << " "<< in_ende  << "\n"; 
	}
        /* wird zurueckgegeben */
	//string seq;
	/* Seitengroesse im RAM */
	long pagesize = getpagesize ();
	
	/* gDNA Sequenz holen */
	
	/* so lang ist die Subsequenz */
	long seq_len = ende - start + 1;

/*--------- dieser Teil muss so sein, weil es sonst nicht funktioniert---*/

	/* keine Ahnung warum, aber das muss so sein */
	ende += in_start_offset;
	
	/* so viel Newlinezeichen sind in der Sequenz vor der
	   gesuchten Subsequenz */
	long nl = start / fasta_line_size;
	/* Startpunkt um Anzahl NewLines verschieben */
	start += nl;
	/* gleiches gilt fuer das Ende */
	nl = ende / fasta_line_size;
	ende += nl;
	
/*--- ich hab aber keine Ahnung warum das so ist -----*/

	/* soviel muss vorne weg */
	long offset_page = start % pagesize;
	/* solang ist die Sequenz im Arbeitsspeicher */
	long len = ende - start + 1 + offset_page;
	/* offset im File muss Vielfaches von pagesize sein */
	long offset_seq = pagesize * (start / pagesize);
	/* enthaelt den Inhalt, des in den speicher kopierten
	   Filesbereichs */
	long maplen = ((len / pagesize) + 1) * pagesize;
	
	
	/* Memory Mapping:
	   data: Adresse, an der der Mappingbereich im Speicher beginnt;
	   0: System bestimmt wo im Speicher der mapping Auschnitt steht;
	   maplen: soviel Speicher wird reserviert;
	   PROT_READ, MAP_SHARD: flags fuer Datei-/Speicherzugriff ;
	   fd_gdna: Filedescriptor (zu mappendes File);
	   offset_seq: nach sovielen Bytes im File wird mit dem Mappen
	   begonnen (muss Vielfaches von pagesize sein);
	*/
	char* data = (char*) mmap ((caddr_t)0, maplen, PROT_READ, MAP_SHARED, file_handle, offset_seq);
	
	/* der Startpunkt der Subsequenz muss ebenfalls um das Start
	   Offset verschoben werden */
	long i = offset_page + in_start_offset;
	
	long char_num=0;
	
	
	/* Sequenz in String speichern */
	//cout << "Sequence length in data " << maplen <<  ", i =" << i << " \n";
	//long my_seq_size=seq.size ();
	//long my_seq_size=in_gdna->size();
	
	

	while (in_gdna->size() < (unsigned long) seq_len) {
	    /* keine Newlines mitspeichern */
	    if (data[i] != '\n'){
		//seq.append (1, data[i]);
		in_gdna->append (1, data[i]);
		char_num++;
		//cout << data[i];
	    }
	    i++;
	}
	
	
	/* Speicherbereich freigeben. Achtung: data ist die Adresse,
	   die bei mmap zurueckgegeben wird! */
	munmap ((caddr_t)data, maplen);
	
	// compute reverse
	i=0;
	long in_seq_size=in_gdna->size ();
	long j=in_seq_size-1;
	char tmp;
	char* my_str=(char*)in_gdna->c_str();
	while(i<j){
	    tmp=my_str[j];
	    my_str[j]=my_str[i];
	    my_str[i]=tmp;
	    i++;
	    j--;
	}
	for(i=0;i<in_seq_size;i++){
	    
	    if(my_str[i]=='a'){
		my_str[i]='t';
		
	    }
	    else if (my_str[i]=='c'){
		my_str[i]='g';
	    }
	    else if (my_str[i]=='g'){
		my_str[i]='c';
	    }
	    else if (my_str[i]=='t'){
		my_str[i]='a';
	    }
	    else if (my_str[i]=='A'){
		my_str[i]='T';
	    }
	    else if (my_str[i]=='C'){
		my_str[i]='G';
	    }
	    else if (my_str[i]=='G'){
		my_str[i]='C';
	    }
	    else if (my_str[i]=='T'){
		my_str[i]='A';
	    }	    
	    
	}
	


//	string tmp_ss;
//	return tmp_ss;
	//return seq;
}






//-------------------------------- Neu Get Sequ --------------------------------//
string align_chain::GetSeqFromFastaFileReverseComplementOld (long in_start, long in_ende,long genome_size, int file_handle, long in_start_offset,long fasta_line_size) {

    long start,ende;
    
    // adjust coordinates
    start=genome_size-in_ende-1;
    ende=genome_size-in_start-1;

	/* wird zurueckgegeben */
	string seq;
	/* Seitengroesse im RAM */
	long pagesize = getpagesize ();

	/* gDNA Sequenz holen */
	
	/* so lang ist die Subsequenz */
	long seq_len = ende - start + 1;

/*--------- dieser Teil muss so sein, weil es sonst nicht funktioniert---*/

	/* keine Ahnung warum, aber das muss so sein */
	ende += in_start_offset;
	
	/* so viel Newlinezeichen sind in der Sequenz vor der
	   gesuchten Subsequenz */
	long nl = start / fasta_line_size;
	/* Startpunkt um Anzahl NewLines verschieben */
	start += nl;
	/* gleiches gilt fuer das Ende */
	nl = ende / fasta_line_size;
	ende += nl;
	
/*--- ich hab aber keine Ahnung warum das so ist -----*/

	/* soviel muss vorne weg */
	long offset_page = start % pagesize;
	/* solang ist die Sequenz im Arbeitsspeicher */
	long len = ende - start + 1 + offset_page;
	/* offset im File muss Vielfaches von pagesize sein */
	long offset_seq = pagesize * (start / pagesize);
	/* enthaelt den Inhalt, des in den speicher kopierten
	   Filesbereichs */
	long maplen = ((len / pagesize) + 1) * pagesize;
	

	/* Memory Mapping:
	   data: Adresse, an der der Mappingbereich im Speicher beginnt;
	   0: System bestimmt wo im Speicher der mapping Auschnitt steht;
	   maplen: soviel Speicher wird reserviert;
	   PROT_READ, MAP_SHARD: flags fuer Datei-/Speicherzugriff ;
	   fd_gdna: Filedescriptor (zu mappendes File);
	   offset_seq: nach sovielen Bytes im File wird mit dem Mappen
	   begonnen (muss Vielfaches von pagesize sein);
	*/
	char* data = (char*) mmap ((caddr_t)0, maplen, PROT_READ, MAP_SHARED, file_handle, offset_seq);
	
	/* der Startpunkt der Subsequenz muss ebenfalls um das Start
	   Offset verschoben werden */
	long i = offset_page + in_start_offset;
	
	long char_num=0;
	
	
	/* Sequenz in String speichern */
	//cout << "Sequence length in data " << maplen <<  ", i =" << i << " \n";
	while (seq.size () <(unsigned long) seq_len) {
	    /* keine Newlines mitspeichern */
	    if (data[i] != '\n'){
		seq.append (1, data[i]);
		char_num++;
		//cout << data[i];
	    }
	    i++;
	}
	
	
	/* Speicherbereich freigeben. Achtung: data ist die Adresse,
	   die bei mmap zurueckgegeben wird! */
	munmap ((caddr_t)data, maplen);
	
	// compute reverse
	i=0;
	long in_seq_size=seq.size ();
	long j=in_seq_size-1;
	char tmp;
	char* my_str=(char*)seq.c_str();
	while(i<j){
	    tmp=my_str[j];
	    my_str[j]=my_str[i];
	    my_str[i]=tmp;
	    i++;
	    j--;
	}
	for(i=0;i<in_seq_size;i++){
	    
	    if(my_str[i]=='a'){
		my_str[i]='t';
		
	    }
	    else if (my_str[i]=='c'){
		my_str[i]='g';
	    }
	    else if (my_str[i]=='g'){
		my_str[i]='c';
	    }
	    else if (my_str[i]=='t'){
		my_str[i]='a';
	    }
	    else if (my_str[i]=='A'){
		my_str[i]='T';
	    }
	    else if (my_str[i]=='C'){
		my_str[i]='G';
	    }
	    else if (my_str[i]=='G'){
		my_str[i]='C';
	    }
	    else if (my_str[i]=='T'){
		my_str[i]='A';
	    }	    
	    
	}
	
	return seq;
}






void align_chain::GetSeqFromFastaFileGeneralized(string* subseq,long in_start, long in_end, long genome_id) {
    	
    long fasta_line_size=fasta_line_size_array[genome_id];
    
    

    if(genome_size_array==NULL){
	    GetSeqFromFastaFile(subseq,in_start,in_end,genome_files_handles[genome_id],start_offset_array[genome_id],fasta_line_size);	    
    }
    else{
	    
	if(direction.c_str()[genome_id]=='p'){
	    GetSeqFromFastaFile(subseq,in_start,in_end,genome_files_handles[genome_id],start_offset_array[genome_id],fasta_line_size);		
	}
	else{
	    GetSeqFromFastaFileReverseComplement(subseq,in_start,in_end,genome_size_array[genome_id],genome_files_handles[genome_id],start_offset_array[genome_id],fasta_line_size);		
	}
    }

}

///////////////////////////////////////////////////////////////////////////////////
/* Wenn alle Fragmente der cDNA eingelesen wurden, wird die gesamte cDNA-
	Sequenz aus dem File gelesen und in cdna_seq gespeichert und der
	entsprechende gDNA Bereich wird gelesen und in gdna_subseq gespeichert
*/

void align_chain::setSeqs () {

    string subseq;
    for(int i=0;i<number_of_genomes;i++){

	GetSeqFromFastaFileGeneralized(&subseq,chain_boundaris[i].start,chain_boundaris[i].end,i);
	genomic_sequences[i].append(subseq);
/*	
	if(genome_size_array==NULL){
	    GetSeqFromFastaFile(&subseq,chain_boundaris[i].start,chain_boundaris[i].end,genome_files_handles[i],start_offset_array[i]);
	    genomic_sequences[i].append(subseq);
	}
	else{
	    //genomic_sequences[i].append(GetSeqFromFastaFile(&subseq,chain_boundaris[i].start,chain_boundaris[i].end,genome_files_handles[i],start_offset_array[i]));
	    if(direction.c_str()[i]=='p'){
		GetSeqFromFastaFile(&subseq,chain_boundaris[i].start,chain_boundaris[i].end,genome_files_handles[i],start_offset_array[i]);
		genomic_sequences[i].append(subseq);
	    }
	    else{
		GetSeqFromFastaFileReverseComplement(&subseq,chain_boundaris[i].start,chain_boundaris[i].end,genome_size_array[i],genome_files_handles[i],start_offset_array[i]);
		genomic_sequences[i].append(subseq);
	    }
	}
*/	
         
    }
    //gdna_subseq.append (GetSeq (g_start, g_ende,0));
    
    return; 


}


/* berechnet die Luecken zwischen den Fragmenten und bearbeitet diese
	entsprechender der zutreffenden Faelle */
void align_chain::process_chain (int in_output_mode,int in_fragment_output_mode,long in_aligned_gap_threshold) {
    
    
    int boundary_flag=0;
    int gap_flag=0;
    FILE* tmp_fptr=NULL;

    output_mode=in_output_mode;
    fragment_output_mode=in_fragment_output_mode;
    threshold_alignable_gap=in_aligned_gap_threshold;

    interval* gap_regions= new interval [number_of_genomes];
   /* alle Fragmente sind eingelesen -> Sequenzen holen und speichern */
    genomic_sequences=new string [number_of_genomes];
    
    cout << "\n \n# Chain no. " << id <<"\n";
    cout << "# Contigs ";
    

    for(int i=0;i<number_of_genomes;i++){
	
	if((palindrome_flag)&&(direction.c_str()[i]=='m')){
	    if(relative_flag==0){
		cout << 1 << " ";
	    }
	    else if(all_contigs_array[i].contig_size_list.size()>0){
		cout << all_contigs_array[i].contig_size_list.size()-contig_id_array[i]+1 << " ";
	    }
	    else{
		cout << 1 << " ";
	    }
	}
	else{
	    cout << contig_id_array[i] << " ";
	    
	}
    }
    cout << endl;

/*
    // displaying chain boundaries
    cout << "# Chain display options: ";
    if(palindrome_flag)
	cout << "palindrome";
    else
	cout << "w.r.t. reverse complement";

    if(relative_flag)
	cout << ", relative positions";
    else
	cout << ", absolute positions";
    cout << endl;
*/
    /*
			

	long j=0;
	for(long i=0;i<number_of_genomes;i++){
	    positions[i].start=atoi(words[j].c_str());
	    positions[i].end=atoi(words[j+1].c_str());
	    j=j+2;
	}
	fragmente.push_back(new fragment(positions,frag_id));
    */

    interval* tmp_positions= new interval [number_of_genomes];
    for(int i=0;i<number_of_genomes;i++){
	tmp_positions[i].start=chain_boundaris[i].start;
	tmp_positions[i].end=chain_boundaris[i].end;
	
    }
    cout << "# Boundaries: ";
    fragment* tmp_frag=new fragment(tmp_positions,0);
    report_fragment(tmp_frag,genome_files_handles[0],-1);    
    delete[] tmp_frag->positions;
    delete tmp_frag;
    tmp_frag=NULL;
    cout << endl;


    // check chain boundaries to decide if to retrieve the whole regions
    boundary_flag=0;
    for(int i=0;i<number_of_genomes;i++){
	
	if((chain_boundaris[i].end-chain_boundaris[i].start)>boundary_threshold){
	    boundary_flag=1;
	    break;
	}
    }
    
    if(boundary_flag==0){
	
	setSeqs();
    }
    
    int two_genomes_flag=0;
    if(number_of_genomes>2){
	two_genomes_flag=0;
    }
    else{
	two_genomes_flag=1;
    }
/*
    //For deubug
    FILE* tmpf;
    tmpf=fopen("tmpseq","a");
    
    for(int i=0;i<number_of_genomes;i++){
	fprintf(tmpf,">gi genome %d, range %d, %d \n", i, chain_boundaris[i].start,chain_boundaris[i].end);
	fprintf(tmpf,"%s\n",genomic_sequences[i].c_str());
    } 
    
    fclose(tmpf);

    printf("Chain %d retrieved\n",id);
    return;
*/     
    
    long chain_size=fragmente.size();
    vector <string> tmp_input_seq;

    /* nur bearbeiten, wenn auch Luecken vorhanden sind */
    if (chain_size > 1) {
	int f_index = 1;
	vector <fragment*>::iterator fragit = fragmente.begin();	
	
	/* solange abarbeiten, bis alle Luecken dran waren */
	
	//cerr << "Chain Size: " << chain_size <<"Last fragment start: " << fragmente[chain_size-1]->positions[1].end << endl;

	while (f_index < chain_size) {
	    fragit++;	    	    
	    // report fragment 
	    //cout << "\n";
	    report_fragment(fragmente[f_index-1],genome_files_handles[0],0); // mode, fragment, file handle, accumulated length
	    
	    // determine length of gaps
	    gap_flag=0;
	    for(int i=0;i<number_of_genomes;i++){
		if((fragmente[f_index]->positions[i].start - 
		    fragmente[f_index-1]->positions[i].end - 1)>
		   threshold_alignable_gap){
		    gap_flag=1;
		    //break;
		}
	    }
	    
	    // end of determine length of gaps
	    if(gap_flag==1){ // gap is not alignable
		unaligned_gap_counter++;
		report_unaligned_gap(fragmente[f_index-1],fragmente[f_index],genome_files_handles[0],0);
		f_index++;		
		continue;
	    }
	    
	    
	    if(two_genomes_flag==0){
		// retrieve ssequences, store them in tmpfile, run clustalw
		tmp_fptr=FOPEN("tmpffXg","w"); // to delete previous contents
		fclose(tmp_fptr);
		tmp_fptr=FOPEN("tmpffXg","a"); // ready to append
	    }
	    
	    for(int i=0;i<number_of_genomes;i++){
		// add header
		
		unsigned long st=fragmente[f_index-1]->positions[i].end + 1;
		unsigned long en=fragmente[f_index]->positions[i].start-1;
		if(two_genomes_flag==0){
		    fprintf(tmp_fptr,">gi%d range %lu, %lu\n",i,st,en);
		}
		else{
		}
		//cout << "\n"<< "start:"<< st << "end: "<< en << "\n";
		if((fragmente[f_index]->positions[i].start -fragmente[f_index-1]->positions[i].end-1)>0){
		    
		    if(boundary_flag==0){ // chain sequences are in main memory
			
			string gdna = "";
			
			gdna.append (genomic_sequences[i].substr (st-chain_boundaris[i].start,en-st+1));
			// format string in fasta format and stroe it in file
			if(two_genomes_flag==0){
			    echosymbolstring2file (tmp_fptr, (char*) gdna.c_str(), gdna.size());
			}
			else{
			    tmp_input_seq.push_back(gdna);
			}
			
			//cout << "\n"<< "start:"<< st << "end: "<< en << "\n";
			//cout << "String: \n"<< gdna.c_str() << "\n";
			
		    }
		    else{
			
			/*string gdna = GetSeqFromFastaFile(fragmente[f_index-1]->positions[i].end + 1,
							  fragmente[f_index]->positions[i].start-1,
							  genome_files_handles[i],start_offset_array[i]);
			*/
			string gdna=""; 
			/*
			GetSeqFromFastaFile(&gdna,fragmente[f_index-1]->positions[i].end + 1,
					    fragmente[f_index]->positions[i].start-1,
					    genome_files_handles[i],start_offset_array[i]);
			*/
			GetSeqFromFastaFileGeneralized(&gdna,fragmente[f_index-1]->positions[i].end + 1,
					    fragmente[f_index]->positions[i].start-1,i);
			
			//cout << "\n"<< "start:"<< st << "end: "<< en << "\n";
			//cout << "String: \n"<< gdna.c_str() <<"\n";
			// format string in fasta format and store it in file
			if(two_genomes_flag==0){
			    echosymbolstring2file (tmp_fptr, (char*) gdna.c_str(), gdna.size());		   
			}
			else{
			    tmp_input_seq.push_back(gdna);
			}
			
		    }
		    
		}
		else{
		    // write no thing
		    if(two_genomes_flag==0){
			fprintf(tmp_fptr,"_\n");
		    }
		    else{
			//gdna="_";
			tmp_input_seq.push_back("");
		    }
		}

	    } // end for :gap sequences retrieved from all genomes
	    if(two_genomes_flag==0){
		fclose(tmp_fptr);
	    }
	    
	    if(two_genomes_flag==0){
		callclustalw ("tmpffXg");
		Multialg tmp_alignment;	    
		readclustalwalignment (fragmente[f_index-1],fragmente[f_index],&tmp_alignment,"tmpffXg.gde",number_of_genomes,0);
	    }
	    else{
		    
	      vector<string> tmp;
	      alignment* align = new alignment (tmp_input_seq[0],tmp_input_seq[1]);
	      align -> constructAlignedSeq ();
	      align -> getAlignedSeq (tmp);
	      /// mohamed
	      unsigned char *allseq = NULL;
	      //allseq=new unsigned char [tmp[0].size()+tmp[1].size+1];
	      long tmp0Size=tmp[0].size();
	      char sep_str[10];
	      sprintf(sep_str,"%c",SEPARATOR);
	      string seperator_string;
	      seperator_string.append(sep_str);
	      //string tmp2=tmp[0]+seperator_string+tmp[1];
	      tmp[0].append(seperator_string+tmp[1]);
	      //allseq=(unsigned char*)tmp2.c_str();
	      allseq=(unsigned char*)tmp[0].c_str();
	      if(output_mode==0){          
		  //output_formatted_alignment_fasta(fragmente[f_index-1],fragmente[f_index],allseq,tmp0Size, tmp2.size(),0);
		  output_formatted_alignment_fasta(fragmente[f_index-1],fragmente[f_index],allseq,tmp0Size, tmp[0].size(),0);
		  //output_formatted_alignment_fasta(in_frag1,in_frag2,allseq, lenalg, numofbytes,accumulated_align_len);
	      }
	      else{       
		  output_formatted_alignment_single_line(fragmente[f_index-1],fragmente[f_index],allseq,tmp0Size, tmp[0].size(),0);
		  //output_formatted_alignment_single_line(fragmente[f_index-1],fragmente[f_index],allseq,tmp0Size, tmp2.size(),0);
	      }
	      delete align;
	      tmp_input_seq.clear();
	      tmp.clear();
	     
	    }
	    
		f_index++; 	    
		
	
	       	    
	} // end while over all fragments in the chain
	// report last fragment 
	//cout << "\n";
	report_fragment(fragmente[f_index-1],genome_files_handles[0],0); // mode, fragment, file handle, accumulated length
    } // end if condition so that the chain contain more than one fragment
    else if(chain_size==1){
	// report the single fragment	
	//cout << "\n";
	report_fragment(fragmente[0],genome_files_handles[0],0);
    }
    report_statistics();

    if(gap_regions!=NULL){
	delete[] gap_regions;//= new interval [number_of_genomes];
	gap_regions=NULL;
    }
    
    if(genomic_sequences!=NULL){
	 for(int i=0;i<number_of_genomes;i++){
	     genomic_sequences[i].empty();
	 }
	delete[] genomic_sequences;//=new string [number_of_genomes];
	genomic_sequences=NULL;
    }
    
}






int align_chain::callclustalw (const char *fasta){

    char command[2*MAXCOMMANDLINELENGTH];

    
	 sprintf (command, "/bin/sh -c \"%s %s %s %s -infile=%s\" > /dev/null",
	    clustalw_name, CLUSTALWBASICPARAMETERS, CLUSTALWOUTPUTPARAMETERS,
	    CLUSTALWTUNINGPARAMETERS, fasta);

    // replaced with clustalw name
	// sprintf (command, "/bin/sh -c \"%s %s %s %s -infile=%s\" > /dev/null",
	//    CLUSTALWNAME, CLUSTALWBASICPARAMETERS, CLUSTALWOUTPUTPARAMETERS,
	//    CLUSTALWTUNINGPARAMETERS, fasta);
    /*
     
    sprintf (command, "/bin/sh -c \"%s %s %s %s -infile=%s\"",
	     CLUSTALWNAME, CLUSTALWBASICPARAMETERS, CLUSTALWOUTPUTPARAMETERS,
	     CLUSTALWTUNINGPARAMETERS, fasta);
    */
    //cout << command << "\n";
    int retvalue=system (command);
    //cout << "Return Value = " << retvalue << "\n";
    return retvalue;
}

void align_chain::echosymbolstring2file (FILE *outf,  char* w, long wlen)
{
  long i;
  for (i = 0; i < wlen; i++)
  {
      if(((i%fasta_line_size_array[0])==0)&&(i>0)){ // for fasta format
	fputc('\n',outf);  
      }  
    fputc (w[i], outf);
  }
  fputc('\n',outf);
}



int align_chain::readclustalwalignment (fragment* in_frag1, fragment* in_frag2,Multialg *alignment, char* gde, long numofsequences, long accumulated_align_len)
{
  unsigned char *file, *allseq = NULL;
  long numofbytes, i, lenalg;
  //int retval;

/*
// for debug
  cout << "HALOOOOOOOOOOOOOOOO "<<endl;
  
  int fd;
  struct stat buf;

  if((fd = open(gde,O_RDONLY)) == -1)
  {
     printf("cannot open \"%s\"",gde);
     return -1;
  }
  if(fstat(fd,&buf) == -1)
  {
     printf("cannot access file status for \"%s\"",gde);
     return -2;
  }
  cout << "FILE SIZE "<< buf.st_size <<endl;
  close (fd);
  return 0;
*/
  // for debug
  
  CONSTRUCTMEMMAP();

//  memmap* memmap_obj;

  file = (unsigned char *) CREATEMEMORYMAP (gde, 0, &numofbytes);
  if (file == NULL)
  {
    printf("Error: cannot create memory map for file %s", gde);
    return -1;
  }


/*
// temporary store allseq
  FILE* allseqf=fopen("AllSeq.align","a");
  fprintf(allseqf,"%s\n",file);
  fclose(allseqf);
  // end testing
  
*/


  allseq=new unsigned char [numofbytes];
  readmultiplegdefileagain (allseq, file, numofbytes);
// first sequence
  for (i = 0; allseq[i] != SEPARATOR; i++)
      ; // just count
  alignment->alignedlength = lenalg = i;


// temporary store allseq
//  FILE* allseqf=fopen("AllSeq.align","a");
//  fprintf(allseqf,"%s\n",allseq);

// format ouptut and compute statistics
  if(output_mode==0){
          
      output_formatted_alignment_fasta(in_frag1,in_frag2,allseq, lenalg, numofbytes,accumulated_align_len);
  }
  else{       
      output_formatted_alignment_single_line(in_frag1,in_frag2,allseq,lenalg,numofbytes,accumulated_align_len);
  }
// end formating output


//  fclose(allseqf);
// end testing
  

  
  delete[] allseq;
  allseq=NULL;

  DELETEMEMORYMAP (file);
  return 0;
}



/*EE
  The following function reads in just the sequence part of a multiple
  gde file. It does not store the positions of the marking symbols and 
  discards the sequence descriptions. Moreover, the input sequence
  is not transformed. It is assumed that the file is in gde format.
  No check for this is performed. 
*/

void align_chain::readmultiplegdefileagain(unsigned char *newseq, unsigned char *input,long inputlen)
{

  unsigned char *inputptr, *newptr = newseq;
  int indesc = 0;

  for(inputptr = input; inputptr < input + inputlen; inputptr++)
  {
    if(indesc)
    {
      if(*inputptr == '\n')
      {
        indesc = 0;
      }
    } else
    {
      if(*inputptr == '#')
      {
        if(inputptr > input)
        {
	    *newptr++ = SEPARATOR; // eliminate
	    //*newptr++ = 'W'; // eliminate
        }
        indesc = 1;
      } else
      {
        if(!isspace((int) *inputptr))
        {
          *newptr++ = *inputptr;
        }
      }
    }
  }
  *newptr++ =0;
}


unsigned char * align_chain::getsequencebynumber (long number, long lenalg, unsigned char *all, long numofbytes)
{
  if ((number+1) * lenalg > numofbytes)
  {
    return NULL;
  }
  // (lenalg+1) skips SEPARATOR
  return all + (number * (lenalg+1));
}


int align_chain::output_formatted_alignment_single_line(fragment* in_frag1, fragment* in_frag2,unsigned char* allseq, long lenalg, long numofbytes, long accumulated_align_len){
    
    long j=0;
    long i=0;
    unsigned char* start_ptr_align_j;
    char* outstr=new char[seq_no_digit+10];    

    
    //long fasta_line_size=fasta_line_size_array[0];
    long fasta_line_size=60;
    cout << "\n";
//    fprintf(allseqf,"Gap:\n");
    // output gap aligned regions
    for(i=0;i<number_of_genomes;i++){
	long len=in_frag2->positions[i].start - in_frag1->positions[i].end -1;
	if((palindrome_flag==0)||(direction.c_str()[i]=='p')){
	    
	    cout <<  len << ":" << in_frag1->positions[i].end+1-contig_offset_array[i] << "-" << in_frag2->positions[i].start-1-contig_offset_array[i] << " ";
	}
	else{
	    	   
	    if(relative_flag==0){
		cout <<  len << ":" << genome_size_array[i]-(in_frag2->positions[i].start-1)-1 << "-" <<  genome_size_array[i]-(in_frag1->positions[i].end+1)-1 << " ";
	    }
	    else{ // palindrome - relative position and the strand is '-'
		long rel_pos_start=in_frag1->positions[i].end+1-contig_offset_array[i];
		long rel_pos_end=in_frag2->positions[i].start-1-contig_offset_array[i];
		// not that this w.r.t. the reversed contigs, we should reverse contig id.
		long ctg_id=all_contigs_array[i].contig_size_list.size()-contig_id_array[i]+1;
		// get contig size
		long ctg_size=all_contigs_array[i].contig_size_list[ctg_id-1];
		cout << len <<  ":" <<ctg_size - rel_pos_end-1 << "-" << ctg_size - rel_pos_start-1  << " ";
	    }
	}

	    
	total_aligned_gap_lengths_in_genomes[i]=total_aligned_gap_lengths_in_genomes[i]+len;
    }
    cout << "\n";
    for(j=0;j<number_of_genomes;j++){
	
	sprintf(outstr,"Seq %*lu:",seq_no_digit,(j+1));
	cout << outstr;

	start_ptr_align_j=getsequencebynumber (j, lenalg, allseq, numofbytes);		 
	if(j==0){                // write the alignment as it is, as if reference seq.
	    for (i = 0; i<lenalg; i++){  
		if(((i+accumulated_align_len)%fasta_line_size)==0){
		    //if(start_ptr_align_j[i]!='\n'){
		    //fprintf(allseqf,"\n");
		    //cout << "\n";
		}
		else if(start_ptr_align_j[i]!='-'){
		    //fprintf(allseqf,"%c",start_ptr_align_j[i]); 
		    cout << start_ptr_align_j[i];
		}
		else{
		    //fprintf(allseqf,"%c",'_');
		    cout << '_';
		}
		
	    }
	    
	}
	else{    
	    // compare character, if not a gap, to the first sequence 
	    for (i = 0; i<lenalg; i++){  
		if((((i+accumulated_align_len)%fasta_line_size)==0)){
		//if(start_ptr_align_j[i]!='\n'){
		    //fprintf(allseqf,"\n");
		    //cout << "\n";
		}
		else if(start_ptr_align_j[i]=='-'){
		    //fprintf(allseqf,"%c",start_ptr_align_j[i]);
		    //fprintf(allseqf,"%c",'_');
		    cout << '_';
		}
		else if(start_ptr_align_j[i]!=allseq[i]){
		    //fprintf(allseqf,"%c",start_ptr_align_j[i]);
		    cout << start_ptr_align_j[i];
		}
		else{
		    //fprintf(allseqf,"%c",'.');
		    cout << '.';
		}
	      }
	}
	cout << "\n";
    }
    
//    fprintf(allseqf,"\n");
    cout << "\n";
    delete[] outstr;
    return 0;
// we could return some statistics later
    
}


int align_chain::output_formatted_alignment_fasta(fragment* in_frag1, fragment* in_frag2,unsigned char* allseq, long lenalg, long numofbytes, long accumulated_align_len){
    
    long j=0;
    long i=0; // i works as offset in allseq
    unsigned char* start_ptr_align_j;
    long k = 0;        
    int term_flag=1;
    char* outstr=new char[seq_no_digit+10];

    //long fasta_line_size=fasta_line_size_array[id];
    long fasta_line_size=60;

      cout << "\n";
//    fprintf(allseqf,"Gap:\n");
    // output gap aligned regions
    for(i=0;i<number_of_genomes;i++){
	long len=in_frag2->positions[i].start - in_frag1->positions[i].end -1;
	if((palindrome_flag==0)||(direction.c_str()[i]=='p')){
	    
	    cout <<  len << ":" << in_frag1->positions[i].end+1-contig_offset_array[i] << "-" << in_frag2->positions[i].start-1-contig_offset_array[i] << " ";
	}
	else{
	    if(relative_flag==0){
		cout <<  len << ":" << genome_size_array[i]-(in_frag2->positions[i].start-1)-1 << "-" <<  genome_size_array[i]-(in_frag1->positions[i].end+1)-1 << " ";
	    }
	    else{
		long rel_pos_start=in_frag1->positions[i].end+1-contig_offset_array[i];
		long rel_pos_end=in_frag2->positions[i].start-1-contig_offset_array[i];
		// not that this w.r.t. the reversed contigs, we should reverse contig id.
		long ctg_id=all_contigs_array[i].contig_size_list.size()-contig_id_array[i]+1;
		// get contig size
		long ctg_size=all_contigs_array[i].contig_size_list[ctg_id-1];
		cout << len <<  ":" <<ctg_size - rel_pos_end-1 << "-" << ctg_size - rel_pos_start-1  << " ";
	    }

	}
	total_aligned_gap_lengths_in_genomes[i]=total_aligned_gap_lengths_in_genomes[i]+len;
    }
    i=0;
    cout << "\n";
    while(term_flag==1){	
	for(j=0;j<number_of_genomes;j++){
	    start_ptr_align_j=getsequencebynumber (j, lenalg, allseq, numofbytes);	    
	    if(j==0){                // write the alignment as it is, as if reference seq.		
		sprintf(outstr,"Seq %*lu:",seq_no_digit,(j+1));
		cout << outstr;
                //cout << "Seq " << (j+1) <<": ";
		
		for (k = 0; ((k+i<lenalg)&&(k<fasta_line_size)); k++){  		    
		    if((k>0)&& (((k+i+accumulated_align_len)%fasta_line_size)==0)) {
			//fprintf(allseqf,"\n");	
			cout << "\t"<< (k+i);
			break;
		    }
		    else if(start_ptr_align_j[k+i]!='-'){
			//fprintf(allseqf,"%c",start_ptr_align_j[k+i]); 
			cout << start_ptr_align_j[k+i];
		    }
		    else{
			//fprintf(allseqf,"%c",'_');
			cout << '_';
		    }		    
		    // check identity
		    check_identity((k+i),lenalg, allseq, numofbytes);
		    //
		}
		if(k==fasta_line_size){
		    cout << "\t"<< (k+i);
		}
		
	    }
	    else{  
		sprintf(outstr,"Seq %*lu:",seq_no_digit,(j+1));
		cout << outstr;
		//cout << "Seq " << (j+1) <<": ";
		// compare character, if not a gap, to the first sequence 
		for (k= 0; ((k+i<lenalg)&&(k<fasta_line_size)); k++){  
		    if((k>0)&& (((k+i+accumulated_align_len)%fasta_line_size)==0)){
			//fprintf(allseqf,"\n");
			break;
		    }
		    else if(start_ptr_align_j[k+i]=='-'){
			//fprintf(allseqf,"%c",start_ptr_align_j[i]);
			//fprintf(allseqf,"%c",'_');
			cout << '_';
		    }
		    else if(start_ptr_align_j[k+i]!=allseq[k+i]){
			//fprintf(allseqf,"%c",start_ptr_align_j[k+i]);
			cout << start_ptr_align_j[k+i];
			
		    }
		    else{
			//fprintf(allseqf,"%c",'.');
			cout << '.';
			
		    }
		}
	    }
	    //fprintf(allseqf,"\n");
	    cout << "\n";
	} // end for j
	i=i+k; // update offset in allseq	
	//fprintf(allseqf,"\n");	
	cout << "\n";
	if(i>=lenalg)break;
    }
    
    delete[] outstr;
    return 0;
// we could return some statistics later
    
}

// function to report fragment
void align_chain::report_fragment(fragment* in_frag,int file_handle, long accumulated_len){ 
// mode, fragment, file handle, accumulated length
    long i;
    
    //cout << "\n";
    for(long i=0;i<number_of_genomes;i++){
	long len=in_frag->positions[i].end - in_frag->positions[i].start +1;
	if(accumulated_len!=-1){
	    total_fragment_lengths_in_genomes[i]=total_fragment_lengths_in_genomes[i]+len;
	}
	if((palindrome_flag==0)||(direction.c_str()[i]=='p')){
	    //cout << "CTOF: "<< contig_offset_array[i]<< ", len ";
	    cout << len << ":" << in_frag->positions[i].start-contig_offset_array[i] << "-" << in_frag->positions[i].end-contig_offset_array[i] << " ";
	}
	else{
	    //cout << len << ":" << in_frag->positions[i].start-contig_offset_array[i] << "-" << in_frag->positions[i].end-contig_offset_array[i] << " ";
	    if(relative_flag==0){
		cout <<  len << ":" << genome_size_array[i]-(in_frag->positions[i].end)-1 << "-" <<   genome_size_array[i]-(in_frag->positions[i].start)-1 << " ";
	    }
	    else{
		long rel_pos_start=in_frag->positions[i].start-contig_offset_array[i];
		long rel_pos_end=in_frag->positions[i].end-contig_offset_array[i];
		// not that this w.r.t. the reversed contigs, we should reverse contig id.
		long ctg_id=all_contigs_array[i].contig_size_list.size()-contig_id_array[i]+1;
		// get contig size
		long ctg_size=all_contigs_array[i].contig_size_list[ctg_id-1];
		cout << len <<  ":" <<ctg_size - rel_pos_end-1 << "-" << ctg_size - rel_pos_start-1  << " ";
	    }
	}
    }
        
    if(accumulated_len==-1){
	cout << "\n";
	return;
    }

    //long fasta_line_size=fasta_line_size_array[0];
    long fasta_line_size=60;
    
    if(fragment_output_mode == 1){ // report seqeunce
	
	cout << "\nExact: \n";
	
/*	string gdna = GetSeqFromFastaFile(in_frag->positions[0].start,
					  in_frag->positions[0].end,file_handle,start_offset_array[0]);
*/	
	string gdna;
//	GetSeqFromFastaFile(&gdna,in_frag->positions[0].start,
//			    in_frag->positions[0].end,file_handle,start_offset_array[0]);
	GetSeqFromFastaFileGeneralized(&gdna,in_frag->positions[0].start,in_frag->positions[0].end,0);
	
	
	long wlen=gdna.size();
	char* w=(char*) gdna.c_str();
	if(output_mode==0){ // mga fasta
	    
	    for (i = 0; i < wlen; i++)
	    {
		if(((i%fasta_line_size)==0)&&(i>0)){ // for fasta format
		    cout << '\n';  
		}  
		cout << w[i];
	    }
	    cout << '\n';
	    
	}
	else{ // lineal
	    
	    for (i = 0; i < wlen; i++)
	    {		
		cout << w[i];
	    }
	    cout << '\n';
	    
	}
    }    

    cout << "\n";
}

void align_chain::report_unaligned_gap(fragment* in_frag1,fragment* in_frag2,int file_handle, long accumulated_len){ 
// mode, fragment, file handle, accumulated length
    cout << "\n";
    long min_len=MY_MAX;
    long max_len=MY_MIN;
    long len;



//    if(output_mode == 0){ // default DNA Mode
	for(long i=0;i<number_of_genomes;i++){
	    len=in_frag2->positions[i].start - in_frag1->positions[i].end -1;
	    if((palindrome_flag==0)||(direction.c_str()[i]=='p')){
		cout <<  len << ":" << in_frag1->positions[i].end+1-contig_offset_array[i] << "-" << in_frag2->positions[i].start-1-contig_offset_array[i] << " ";
	    }
	    else{
		//cout <<  len << ":" << in_frag1->positions[i].end+1-contig_offset_array[i] << "-" << in_frag2->positions[i].start-1-contig_offset_array[i] << " ";
		if(relative_flag==0){
		    cout <<  len << ":" << genome_size_array[i]-(in_frag2->positions[i].start-1)-1 << "-" <<  genome_size_array[i]-(in_frag1->positions[i].end+1)-1 << " ";
		}
		else{
		    long rel_pos_start=in_frag1->positions[i].end+1-contig_offset_array[i];
		    long rel_pos_end=in_frag2->positions[i].start-1-contig_offset_array[i];
		    // not that this w.r.t. the reversed contigs, we should reverse contig id.
		    long ctg_id=all_contigs_array[i].contig_size_list.size()-contig_id_array[i]+1;
		    // get contig size
		    long ctg_size=all_contigs_array[i].contig_size_list[ctg_id-1];
		    cout << len <<  ":" <<ctg_size - rel_pos_end-1 << "-" << ctg_size - rel_pos_start-1  << " ";
		}
	    }
	    if(min_len>len)min_len=len;
	    if(max_len<len)max_len=len;
	    total_unaligned_gap_lengths_in_genomes[i]=total_unaligned_gap_lengths_in_genomes[i]+len;
	}
	cout << "\nGap "<< unaligned_gap_counter << ":\t" << min_len << " (min.length) + ";
	cout << (max_len - min_len) << " (diff.) = " << max_len << " (max.length)";
//    }    
    
    cout << "\n";
}

void align_chain::check_identity(long position, long in_lenalg,unsigned char* in_allseq,long in_numofbytes){
    
    //unsigned char ref_char=in_allseq[position];  
    unsigned char ref_char=getsequencebynumber (0, in_lenalg, in_allseq, in_numofbytes)[position];
    unsigned char* start_ptr_align_i;    
    if ((ref_char=='n')||(ref_char=='N'))return;
    for(long i=1;i<number_of_genomes;i++){
	 start_ptr_align_i=getsequencebynumber (i, in_lenalg, in_allseq, in_numofbytes);
	 if(ref_char!=start_ptr_align_i[position]){
	     return;
	 }
    }
    
    total_identity++;
}

void align_chain::report_statistics(){

    cout << "\n# Statistics:\n";
    cout << "# Number of unaligned gaps: " << unaligned_gap_counter << "\n";
    cout << "# Total Match length:\n";
    for(long i=0;i<number_of_genomes;i++){
	cout << "# Seq: "<<(i+1)<<" :" <<total_fragment_lengths_in_genomes[i]<< "\n";
    }
    
    cout << "# Total aligned gap length:\n";
    for(long i=0;i<number_of_genomes;i++){
	cout << "# Seq: "<<(i+1)<<" :" <<total_aligned_gap_lengths_in_genomes[i]<< "\n";
    }

    cout << "# Total unaligned gap length:\n";
    for(long i=0;i<number_of_genomes;i++){
	cout << "# Seq: "<<(i+1)<<" :" <<total_unaligned_gap_lengths_in_genomes[i]<< "\n";
    }
    cout << "# Total identity in aligned gaps: "<< total_identity << "\n";

    cout << "# Identity ratio of chain: ";
    
    long len;
    float per_cent;
    char tmp_char[10];
    for(long i=0;i<number_of_genomes;i++){
	len=chain_boundaris[i].end-chain_boundaris[i].start+1;
	per_cent= (float)(total_fragment_lengths_in_genomes[i]+total_identity)/(float)len;
	sprintf(tmp_char, "%f",per_cent);
	cout << tmp_char << " ";
    }
    cout << endl;
    
    
}
/////////////////////////////////////////////////////////////////////////////
// to check if the chains are in increasing ordered 
// that is if the program invert_chain is applied or not
int align_chain::check_chain_file_consistency(){
    
    if(fragmente.size()>1){
	if (fragmente[0]->positions[0].start > fragmente[1]->positions[0].start)
	    return (1);
	//else return (0);
    }
    
    // checking boundaries against file size

    long in_start;
    long in_end;
    
    int res;
    //cerr << "Chain boundary: " << chain_boundaris[0].start << " " << chain_boundaris[0].end <<"\n";
    if(genome_size_array==NULL){
	
	for(int i=0;i<number_of_genomes;i++){
	    long fasta_line_size=fasta_line_size_array[i];
	    in_start=chain_boundaris[i].start;
	    in_end=chain_boundaris[i].end;
	    res=check_chain_boundary_against_file_size(in_start,in_end,genome_files_handles[i],start_offset_array[i],fasta_line_size);
	    if(res!=0){
		cerr << "Original s: " << in_start << " Original e: " << in_end << " dirc: "<< direction.c_str()[i] << "\n";
		return(2);
	    }
	}
    }
    else{

	for(int i=0;i<number_of_genomes;i++){
	    long fasta_line_size=fasta_line_size_array[i];
	    in_start=chain_boundaris[i].start;
	    in_end=chain_boundaris[i].end;
	    
	    if(direction.c_str()[i]=='p'){
		res=check_chain_boundary_against_file_size(in_start,in_end,genome_files_handles[i],start_offset_array[i],fasta_line_size);
		if(res!=0){
		    cerr << "Original s: " << in_start << "Original e: " << in_end << " dirc: "<< direction.c_str()[i] << "\n";
		    return(2);
		}
	    }
	    else{
		long start=genome_size_array[i]-in_end-1;
		long end=genome_size_array[i]-in_start-1;
		res=check_chain_boundary_against_file_size(start,end,genome_files_handles[i],start_offset_array[i],fasta_line_size);
		if(res!=0) {
		    cerr << "Original s: " << in_start << "Original e: " << in_end << " dirc: "<< direction.c_str()[i] << "\n";
		    return(2);
		}
	    }
	}
    }
    
    return(0);
}

/////////////////////////////////////////////////////////////////////////////////
int align_chain::check_chain_boundary_against_file_size (long start, long ende, int file_handle, long in_start_offset, long in_fasta_line_size) {
	
        struct stat buf;
	long numofbytes;
	if(fstat(file_handle,&buf) == -1)
	{
	    printf("cannot access file status\n");
	    return -2;
	}
	numofbytes = (long) buf.st_size;

	
	long pagesize;
	pagesize = getpagesize ();
	
	/* gDNA Sequenz holen */
	
	/* so lang ist die Subsequenz */
	long seq_len;
	seq_len= ende - start + 1;

/*--------- dieser Teil muss so sein, weil es sonst nicht funktioniert---*/

	/* keine Ahnung warum, aber das muss so sein */
	ende += in_start_offset;

//	cerr << ende << " && " <<numofbytes << endl;
	
	/* so viel Newlinezeichen sind in der Sequenz vor der
	   gesuchten Subsequenz */
	long nl = start / in_fasta_line_size;
	/* Startpunkt um Anzahl NewLines verschieben */
	start += nl;
	/* gleiches gilt fuer das Ende */
	nl = ende / in_fasta_line_size;
	ende += nl;
	
	if(ende>numofbytes){
	    cerr << "Start: " << start << " End: " << ende << " file size " << numofbytes << "\n";
	    return(-1);
	}
	return(0);

}



