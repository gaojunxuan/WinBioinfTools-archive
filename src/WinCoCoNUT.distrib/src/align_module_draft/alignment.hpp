#ifndef alignment_h_INCLUDED
#define alignment_h_INCLUDED

#include <string>

using namespace std;


class alignment {
private:
	/* Koordinate in der Alignmentmatrix */
	struct koordinate {
		int x, y;
	};

	/* Koordinate rechts unten */
	struct koordinate l;

	int zeilen;
	int spalten;

	/* s: cDNA Sequenz, t: gDNA Sequnenz */
	string s;
	string t;

	/* alignierte Sequenzen */
	string s_aligned;
	string t_aligned;

	/* Alignment Matrix */
	int** matrix;

	/* anthält alle moeglichen alignierten Sequenzen */
	vector<string> alignments;

	/* score des Alignment aus s und t */
	int score;

	void Ausgabe_Alignment ();

	/* Alignment zurueckverfolgen und die Zeichen in die Liste speichern */
	void GetAlignment (struct koordinate k);

	/* berechnet den Score der beiden Zeichen */
	int Score (char a, char b);

	/* Berechnet das Maximum fuer eine Zelle der Matrix bei der
		Berechnung des Alignment-Scores */
	int Max (int i, int j);

public:
	/* Konstrukteur */
	alignment (string sx, string tx);

	/* gibt den score aus */
	int getScore ();

	/* berechnet alle moeglichen alignierten Sequenzen */
	void constructAlignedSeq ();

	/* gibt alle alignierten Sequenzen zurueck */
	void getAlignedSeq (vector<string> &alignedSeq);

	/* Dekonstrukteur */
	~alignment ();
};

#endif
