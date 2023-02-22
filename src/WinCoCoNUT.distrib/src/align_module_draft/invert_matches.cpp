#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


//char orig_matches[100];
//char* const neu_matches ("inv_org_matches.ctg");

void help (char* prog) {
	cout << "usage: " << prog << " <Chainer-Ausgabe-Datei>\n";
	
	exit (9);
}


int main (int argc, char* argv[]) {

    
        string output_file_name;
	if (argc != 2){
	    //strcpy (orig_matches, argv[1]);
	    //else
	    help (argv[0]);
	}
	
	vector <string> zeilen;
	vector <string>::iterator x;
	string zeile;
	int gelesen = 0, geschrieben = 0;

//	ifstream in ("../Sourcefiles/matchesChr19Fantom20.vm.chn.ctg", ios::in);
//	ifstream in ("../arabdopsis/matchchr1L20.vm.mat.chn.ctg2", ios::in);
	ifstream in (argv[1], ios::in);
	if (!in) {
		cerr << "Datei konnte nicht geoffnet werden\n";
		exit (1);
	}

//	ofstream out ("../Sourcefiles/inv_org_matches.ctg", ios::out);
//	ofstream out ("../arabdopsis/inv_org_matches2.ctg", ios::out);
	output_file_name.append(argv[1]);
	output_file_name.append(".ordered");
	ofstream out (output_file_name.c_str(), ios::out);
	if (!out) {
		cerr << "Datei konnte nicht angelegt werden\n";
		exit (2);
	}
	// read header
	getline (in, zeile);
	out << zeile << endl;
	while (getline (in, zeile)) {
	    zeilen.push_back (zeile);
		gelesen++;
	}

	while (!zeilen.empty ()) {
		int i = zeilen.size () - 1;
		x = zeilen.end ();
		while (zeilen[i].substr (0, 1).compare ("#") != 0) {
			i--;
			x--;
		}
		out << zeilen[i] << endl;
		geschrieben++;
		int j=zeilen.size ()-1;
		//while (i < zeilen.size ()) {
		//	out << zeilen[i] << endl;
		while (i < j) {
			out << zeilen[j] << endl;
			geschrieben++;
			j--;
		}

		while (x != zeilen.end ())
			zeilen.pop_back ();
		zeilen.pop_back ();
	}
	//cout << "gelesen     : " << gelesen << endl;
	//cout << "geschrieben : " << geschrieben << endl;

	in.close ();
	out.close ();

	if(zeilen.size()>0)
	    zeilen.clear();
	return 0;
}
