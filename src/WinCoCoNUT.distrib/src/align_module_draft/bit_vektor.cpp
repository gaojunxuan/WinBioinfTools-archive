#include "bit_vektor.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <ctype.h>


/*
 *	Berechnet basis^exp
 */
int bit_vektor::hoch (int basis, int exp) {
	int erg = 1;
	for (int i = 0; i < exp; i++)
		erg *= basis;

	return erg;
}


/*
 *	Initialisiert alle Variablen, die nicht mehr geaendert werden
 *	oder initial eine besondere Form haben
 */
void bit_vektor::init () {
	/* fuer Peq */
	v.push_back ('A');
	v.push_back ('C');
	v.push_back ('G');
	v.push_back ('T');

	/* Wortlaenge, hier 32 */
	max_len = 32;
	/* Anzahl Bloecke pro Spalte (s.length () DIV 32 + 1) */
	durchlaufe = 0;

	/* erster Score hat s.length Differenzen */
	Score = s.length ();

	/* Kopie von s */
	string s_tmp;
	s_tmp.append (s);
	/* Teilstrings fuer Bloecke ermitteln */
	while (s_tmp.length () > max_len) {
		/* Teilstring in Vector speichern */
		s_geteilt.push_back (s_tmp.substr (0, max_len));
		/* Teilstring aus Kopie loeschen */
		s_tmp.erase (0, max_len);
		/* pro Teilstring ein Durchlauf mehr */
		durchlaufe++;
	}
	/* falls was uebrig bleibt noch in Vector speichern */
	if (s_tmp.length () > 0)
		s_geteilt.push_back (s_tmp);
	/* ein Durchlauf findet immer statt */
	durchlaufe++;

	/* horizontale Deltas in der ersten Zeile eines Blocks
		enthaelt nach dem Durchlauf die horizontalen Deltas der 
		letzten Zeile eines Blocks und somit die erste fuer den
		naechsten Block */
	for (int i = 0; i < t.length (); i++) {
		/* initial immer horizontales Delta = +1 */
		M_save.push_back (0);
		P_save.push_back (1);
	}
}


/*
 *	Vorverarbeitung des Teilstrings von s pro Block
 *	pro Block ein neues Peq
 */
void bit_vektor::precompute (int k) {
	for (int i = 0; i < v.size (); i++) {
		/* initial 0 fuer jedes Symbol */
		peq[i] = 0;
		for (int j = 0; j < s_geteilt[k].length (); j++)
			/* der Position des Symbols entsprechend das Bit auf
				1 setzen */
			if (s_geteilt[k].at(j) == v[i])
				peq[i] = peq[i] | hoch(2, j);
	}
}


/*
 *	fuer neuen Block initialisieren
 */
void bit_vektor::neu (int k) {
	/* Spaltenvector der ersten Spalte immer alles auf +1 */
	Pv = hoch (2, s_geteilt[k].length ()) - 1;
	Mv = 0;

	/* soviele Stellen werden in diesem Block verwendet */
	stellen = Pv;
	/* fuer die Scoreanpassung wichtig */
	score_schablone = hoch (2, s_geteilt[k].length () - 1);
}


/*
 *	gibt die Integer fuer Symbol c zurueck
 */
int bit_vektor::getPeq (char c) {
	for (int i = 0; i < v.size (); i++)
		if (c == v[i])
			return peq[i];

	return -1;
}


/*
 *	setzt den Score (rScore: den richtigen Score)
 *	dif ist die Anzahl Unterschiede in x und y
 */
void bit_vektor::setScore (string x, string y, int dif) {
	/* Laenge vom laengeren String */
	int len;
	if (x.length () > y.length ())
		len = x.length ();
	else
		len = y.length ();

	/* neutrale Elemente werden als Match betrachtet, 
		darum Anzahl neutrale Elemente vom Score wieder
		abziehen */
	int a = 0, b = 0, anz = 0;
	while ((a = t.find_first_not_of ("ACGT", b)) != -1) {
		b = a + 1;
		anz++;
	}

	/* Score = Laenge von max(x.length, y.length) - 2*dif
				- anz.neutr.El. */
//	rScore = len - 2*dif - anz;
	rScore = 2*len - 3*dif - 2*anz;
}


/*
 *	Konstrukteur + sofortige Berechnung des Scores
 */
bit_vektor::bit_vektor (string sx, string tx) {
	/* beide Strings in Grossschreibung */
	s.append (sx);
	for (int i = 0; i < s.length (); i++)
		s.at(i) = toupper (s.at(i));

	t.append (tx);
	for (int i = 0; i < t.length (); i++)
	    t.at(i) = toupper (t.at(i));
	
	/* Werte initialisieren */
	init ();

	/* ein durchlauf fuer jeden Block / Teilstring von s */
	for (int k = 0; k < durchlaufe; k++) {
		/* fuer aktuellen Block neu initialisieren */
		neu (k);
		/* Peq fuer Block berechnen */
		precompute (k);

		/* Durchlauf fuer jeden Block durch t */
		for (int j = 0; j < t.length (); j++) {
/***************** vgl. Bit-Vektor-Alg. **********************/
			int Eq = getPeq (t.at (j));

			/* muss fuer geteilte s sein */
			if (M_save[j] == 1)
				Eq |= 1;

			int Xv = Eq | Mv;
			int Xh = ((((Eq & Pv) + Pv) ^ Pv) | Eq) & stellen;

			int Ph = (Mv | (~ (Xh | Pv))) & stellen;
			int Mh = Pv & Xh;

			/* speziell fuer die Scoreanpassung */
			int blubb1 = Ph & score_schablone;
			int blubb2 = Mh & score_schablone;

			/* speichert Eintrag fuer P_save und M_save */
			int save_M = 0;
			int save_P = 0;

			/* Score fuer Delta h = +1 anpassen */ 
			if (blubb1 == score_schablone) {
				/* Score nur im letzten Durchlauf / Block anpassen */
				if (k == durchlaufe - 1)
					Score += 1;
				save_P = 1;
			}
			/* Score fuer Delta h = -1 anpassen */
			else if (blubb2 == score_schablone) {
				/* Score nur im letzten Durchlauf / Block anpassen */
				if (k == durchlaufe - 1)
					Score -= 1;
				save_M = 1;
			}

			Ph <<= 1;
			Ph |= P_save[j];
			Ph &= stellen;

			Mh <<= 1;
			Mh |= M_save[j];
			Mh &= stellen;

			Pv = (Mh | (~ (Xv | Ph))) & stellen;
			Mv = Ph & Xv;
//*******************************************************************
			/* erstes Delta h fuer naechsten Block speichern */
			M_save[j] = save_M;
			P_save[j] = save_P;
		}
	}

	/* Score abspeichern */
	setScore (s, t, Score);

//	cout << "Score: " << Score << endl;
//	cout << "Score2: " << get_Score () << endl;
}


/*
 *	gibt den Score zurueck 
 */
int bit_vektor::get_Score () {
	return rScore;
}


/*
 *	Dekonstrukteur
 */
bit_vektor::~bit_vektor () {

}
