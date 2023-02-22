#include <iostream>
#include <fstream>
#include <vector>

#include "fragment.hpp"
//#include "alignment.hpp"
//#include "greedy_align.hpp"
//#include "bit_vektor.hpp"

using namespace std;

/* globale Variable fuer die ausgabedatei */
extern ofstream out;
/* Alignment-Verfahren */
extern int av;

//---------------zu Testzwecken---------------------
extern int gt_ag[];
extern int ct_ac[];
//--------------------------------------------------

/* Konstrukteur */
fragment::fragment (int i = -1, int j = -1, int k = -1, int l = -1) {
	c_start = i;
	c_ende = j;
	g_start = k;
	g_ende = l;
	korrektur = "";
	positions=NULL;
}   

fragment::fragment (interval* in_position,long in_id) {
    positions=in_position;
    id=in_id;
    korrektur = "";
}   


/* aendert Intervallgrenzen */
void fragment::set (int i = -1, int j = -1, int k = -1, int l = -1) {
	c_start = i;
	c_ende = j;
	g_start = k;
	g_ende = l;
}




/* Dekonstrukteur */
fragment::~fragment () {
    if(positions != NULL){
	//delete[] positions;
	//positions=NULL;
    }
}


