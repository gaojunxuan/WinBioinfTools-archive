#include <stdio.h>
//#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "math.h"

using namespace std;
//using std::cout;
//#include "assert.h"

#define max(a, b)  ((a)>=(b) ? (a) : (b))

#define FAIL 1
#define SUCCESS 0


vector <string> extension_array;
void fsplit (string line, const char* sep, vector<string> & words);

void create_files_for_all_suffixes(string file_name_prefix,long dimension);

int main(int argc,char *argv[]){

  long match_no=0;

  if(argc!=3)
    {
	cout << "Invalid options\n";
	cout << "Usage: multimat2chainer filename orientation\n";
	cout << "orientation: strand orientations except for first seq in +|- \n";    
	cout << "Example for 3 sequences with the -ve strand of the 3d:\n";
	cout << "> multimat2chainer Seqfor3genomes ++ \n";
	cout << "Output:\n ";
	cout << "-One file: 'filename.pm' in case of multimat input file \n";    
    return(FAIL);
  } 
  
  
  ifstream inputf(argv[1], ios::in);
  if (!inputf) {
      cerr << "Error in opening input file!\n";
      exit (1);
  }
  string input_file_str=argv[1];

  // create output file name
  // string in_suffix="";
  // in_suffix.append((char*)argv[2]);
  long number_of_genomes=atol(argv[2]);
  cout << "Input file: " << argv[1] << "\n";
  cout << "Number of genomes: " << number_of_genomes << endl;
    
  vector <long> query_size_array;
/*
  string suffix="p";  
  for(long i=0;i<in_suffix.size();i++){
      if(in_suffix.c_str()[i]=='+'){
	  suffix.append("p");
      }
      else{
	  suffix.append("m");
      }
  }

  string outfilename=argv[1];
  outfilename.erase(outfilename.size()-number_of_genomes+1);
  outfilename.append(".");
  outfilename.append(suffix);

  ofstream outputf(outfilename.c_str(), ios::out);
  if (!outputf) {
      cerr << "Error in opening file" << outfilename << "!\n";
      exit (1);
  }
  outputf << ">CHA " << number_of_genomes << "\n";
  outputf.close();
*/

  create_files_for_all_suffixes(input_file_str,number_of_genomes);
  
  //cout << "In Suffix: +" << in_suffix << "Stored Suf: "<< suffix << endl;
  string line;
  vector< string> words; 
  long length;
  vector<long> positions; 

  while(getline (inputf, line)){
      string direction="p";
      if(line.substr(0,1).compare("#")==0){	  
	  words.clear(); 
	  fsplit(line, " ", words);
	  if(words.size()>3){
	      //cout << line << endl;
	      if(words[1].substr(0,4).compare("file")==0){
		  //cout <<  "Query Seq. Size: " << words[words.size()-1] << endl;
		  query_size_array.push_back(atol(words[words.size()-2].c_str()));
		  cout <<  "Query Seq. Size: " << query_size_array[query_size_array.size()-1] << endl;
	      }
	  }
      }
      else{ // handling matches
	  words.clear(); 
	  fsplit(line, " ", words);words[0].c_str();
	  length=atol(words[0].c_str());
	  // get file extension and the exact ouput file
	  string tmp_suf=".";
	  for(long i=1;i<words.size();i++){
	      if(words[i].c_str()[0]!='-'){ // not reverse complement
		  tmp_suf=tmp_suf+"p";
	      }
	      else{
		  tmp_suf=tmp_suf+"m";
	      }
	  }
	  string outfilename=input_file_str;
	  outfilename.append(tmp_suf);
	  ofstream outputf(outfilename.c_str(), ios::app);
	  if (!outputf) {
	      cerr << "Error in opening file" << outfilename << "!\n";
	      exit (1);
	  }	  	  
	  outputf << "#" << length << endl;

	  for(long i=1;i<words.size();i++){
	      if(words[i].c_str()[0]!='-'){ // not reverse complement
		  //positions.push_back(atol(words[i].c_str()));
		  long tmpnum=atol(words[i].c_str());
		  outputf << "[" << tmpnum << "," << tmpnum+length-1 << "] ";
		  //direction.append("p");
	      }
	      else{		  
		  long tmp_pos1=atol(words[i].substr(1,words[i].size()-1).c_str());
		  long tmp_pos2=tmp_pos1+length-1;
		  tmp_pos1=query_size_array[i-1]-tmp_pos1-1;
		  tmp_pos2=query_size_array[i-1]-tmp_pos2-1;
		  outputf << "[" << tmp_pos2 << "," << tmp_pos1 << "] ";
		  //positions.push_back(tmp_pos);
		  //direction.append("m");
	      }
	  }
	  match_no++;
	  outputf << "\n";	 
	  outputf.close();
      }
  }   
  inputf.close();


  // creating info file
  string tmp_suf=".info";
  string outfilename=input_file_str;
  outfilename.append(tmp_suf);
  ofstream outputf(outfilename.c_str(), ios::out);
  if (!outputf) {
      cerr << "Error in opening file" << outfilename << "!\n";
      exit (1);
  }	  	  
  outputf << number_of_genomes <<" ";
  for(long i=0;i<number_of_genomes;i++){
      outputf<<query_size_array[i] <<" ";
  }
  outputf << "\n";
  for(long i=extension_array.size()-1;i>=0;i--){
      outputf<< input_file_str+extension_array[i] <<endl;
  }  
  outputf.close();


  return(0);
}


void fsplit (string line, const char* sep, vector<string> & words) {
		
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


void create_files_for_all_suffixes(string file_name_prefix,long dimensions){
    
    long no_files=(long)exp2((double)(dimensions-1));
    long buffer;
    int shift=1;
    //cerr << "No. of files: " << no_files << endl;
    extension_array.clear();
    for(long i=0;i<no_files;i++){
	buffer=i; 
	//cerr << "buffer: "<< buffer << endl;
	string tmp=".p";
	//for(long j=0;j<dimensions-1;j++){
	for(long j=dimensions-2;j>=0;j--){
	    //if(buffer & (1 << j)) {
	    if(buffer & (1 << j)) {
		//cerr << "1";
		tmp.append("p");
	    } else {
		//cerr << "0";
		tmp.append("m");
	    }	    
	}
	//cerr << tmp <<endl;
	extension_array.push_back(tmp);
    }
    
    for(long i=0;i<no_files;i++){
	string extension=extension_array[i];
	string synteny_file_name=file_name_prefix+extension;
	
	//FILE* fptr=fopen((char*)synteny_file_name.c_str(),"r");
	//if(fptr==NULL){
	FILE* fptr=fopen((char*)synteny_file_name.c_str(),"w");
	if(fptr==NULL){
	    cerr << "Error: Unable to open output file"<< synteny_file_name << endl;
	    exit(-1);
	}
	fprintf(fptr,">CHA %lu\n",dimensions);
	//cout << "haloooooooooooooooo\n"<< dimsnions << endl;
	fclose(fptr);


    }
}
