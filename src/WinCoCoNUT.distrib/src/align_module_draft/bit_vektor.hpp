#ifndef bit_vektor_h_INCLUDED
#define bit_vektor_h_INCLUDED

#include <string>
#include <vector>

using namespace std;


class bit_vektor {
private:
	/* Sequenzen, fuer die der Score berechnet werden soll */
	string s, t;
	/* Vorverarbeitete Teilsequenz fuer den aktuellen Block */
	int peq[4];
	/* Hilfe fuer Peq */
	vector<char> v;
	/* kodieren die vertikalen Delta fuer einen Block */
	int Pv, Mv;
	/* Score = Distanz; rScore = richtiger Score */
	int Score, rScore;
	/* betrachtete Stellen pro Block */
	int stellen;
	/* Maske um den Score anzupassen */
	int score_schablone;
	/* kodieren vor dem Durchlauf die horizontalen Deltas der ersten
		Zeile fuer den aktuellen Block und nach dem Durchlauf die
		horizontalen Deltas fuer den naechsten Block. Werden
		on-the-fly fuer den naechsten Block aktualisiert */
	vector<int> P_save, M_save;
	/* enthaelt alle Teilstrings fuer die Bloecke */
	vector<string> s_geteilt;
	/* Wortlaenge; hier 32 */
	int max_len;
	/* Anzahl Durchlaeufe bzw. Bloecke */
	int durchlaufe;

	/* Berechnet basis^exp */
	int hoch (int basis, int exp);

	/* Initialisiert alle Variablen, die nicht mehr geaendert werden
		oder initial eine besondere Form haben */
	void init ();

	/* Vorverarbeitung des Teilstrings von s pro Block 
		pro Block neuer Peq */
	void precompute (int k);

	/* initialisiert fuer neuen Block */
	void neu (int k);

	/* setzt den Score */
	void setScore (string x, string y, int dif);

	/* gibt Integer fuer Symbol c zurueck */
	int getPeq (char c);

public:
	/* Konstrukteur + Algorithmus */
	bit_vektor (string sx, string tx);

	/* gibt den Score zurueck */
	int get_Score ();

	/* Dekonstrukteur */
	~bit_vektor ();
};

#endif

