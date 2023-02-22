#ifndef fragment_h_INCLUDED
#define fragment_h_INCLUDED

#include <string>

using namespace std;

class interval{
public:
    long start;
    long end;
};


class fragment {
public:
    /* Intervallangaben */
    int c_start, g_start, c_ende, g_ende;
    class interval* positions;    
    
    /* vorgenommene Korrekturen am Fragment */
    string korrektur;
    /* Name der cDNA */
    string name;
    long id;
    
    /* Konstrukteur */
    fragment (int i, int j, int k, int l);
    fragment (class interval* in_position,long id);
    
    /* aendert Intervallgrenzen */
    void set (int i, int j, int k, int l);

	/* Ausgabe in Ausgabedatei */
	void print (int id);

	/* Dekonstrukteur */
	~fragment ();
};

#endif



