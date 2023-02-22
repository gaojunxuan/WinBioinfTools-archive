#include <iostream>
#include <vector>

#include "alignment.hpp"

using namespace std;


/* Alignment zurueckverfolgen und die Zeichen in die Liste speichern */
void alignment::GetAlignment(struct koordinate k) {
	struct koordinate akut;
	int w_akut, w_diag, w_hoch, w_links;
	int i;
//	int count = 0;

	/* Solange durchlaufen, bis eine Sequenz vollstaendig ist */
	while ((k.x != 0) && (k.y != 0)) {
		/* Werte der aktuellen Koordinate in der Matrix */
		akut.x = k.x;
		akut.y = k.y;
		w_akut = matrix[k.x][k.y];

		/* Werte der diagonalen Koordinate in der Matrix */
		w_diag = matrix[k.x-1][k.y-1];

		/* Werte der hoeheren Koordinate in der Matrix */
		w_hoch = matrix[k.x-1][k.y];

		/* Werte der linken Koordinate in der Matrix */
		w_links = matrix[k.x][k.y-1];

		/* kam der aktuelle Wert durch ein match (+2) oder ein mismatch (-1)
		   ueber die Diagonale zu stande */
		if (((w_diag == w_akut-2) && (s[akut.x-1] == t[akut.y-1])) ||
				((w_diag == w_akut+1) && (s[akut.x-1] != t[akut.y-1])) ||
				/* oder durch ein N auf t (0) */
				((w_diag == w_akut) && (t[akut.y-1] == 'N'))) {

			/* bei der diagonalen Koordinate gehts weiter */
			k.x--;
			k.y--;
			/* beide Zeichen in Listen speichern */
			s_aligned.insert (s_aligned.begin (), s[akut.x-1]);
			t_aligned.insert (t_aligned.begin (), t[akut.y-1]);
		}

		/* kam der aktuelle Wert durch einen gap der waagerechten Sequenz
		   zu stande (-1), oder durch einen gap und ein N auf t (0) */
		else if ((w_hoch == w_akut + 1) ||
				((w_hoch == w_akut) && (t[akut.y-1] == 'N'))) {

			/* bei der oberen Koordinate gehts weiter */
			k.x--;
			/* Zeichen der senkrechten Sequenz und das gapzeichen in Listen
			   speichern */
			s_aligned.insert (s_aligned.begin (), s[akut.x-1]);
			t_aligned.insert (t_aligned.begin (), '_');
		}

		/* kam der aktuelle Wert durch einen gap der senkrechten Sequenz
		   zu stande (-1), oder durch einen gap und ein N auf t (0) */
		else if ((w_links == w_akut + 1) ||
				((w_links == w_akut) && (t[akut.y-1] == 'N'))) {

			/* bei der linken Koordinate gehts weiter */
			k.y--;
			/* Zeichen der senkrechten Sequenz und das gapzeichen in Listen
			   speichern */
			s_aligned.insert (s_aligned.begin (), '_');
			t_aligned.insert (t_aligned.begin (), t[akut.y-1]);
		}

	} // If

	/* eine der beiden oder beide Sequenzen sind abgearbeitet
	   den Rest mit gaps alignieren */
	if (((k.x == 0) && (k.y > 0)) || ((k.x >0) && (k.y == 0))) {
		/* senkrechte Sequenz ist fertig */
		if (k.x == 0) {
			for (i=k.y; i>0; i--) {
				/* gaps und Zeichen speichern */
				s_aligned.insert (s_aligned.begin (), '_');
				t_aligned.insert (t_aligned.begin (), t[i-1]);
			}
		}
		/* waagerechte Sequenz ist fertig */
		else {
			for (i=k.x; i>0; i--) {
				/* Zeichen und gaps speichern */
				s_aligned.insert (s_aligned.begin (), s[i-1]);
				t_aligned.insert (t_aligned.begin (), '_');
			}
		}
	}

	/* alignierte Sequenzen in Vektor speichern */
	alignments.push_back (s_aligned);
	alignments.push_back (t_aligned);
}


/* berechnet den Score der beiden Zeichen */
int alignment::Score (char a, char b) {
	int test;

	/* gleiches Zeichen -> match */
	if (a == b)
		test = 2;
	/* neutrales Element gibt keine Bestrafung */
	else if (b == 'N')
		test = 0;
	/* verschiedene Zeichen -> mismatch oder gap */
	else
		test = -1;
	return test;
}

/* Berechnet das Maximum fuer eine Zelle der Matrix bei der
   Berechnung des Alignment-Scores */
int alignment::Max (int i, int j) {
	int maxi;

	/* match */
	int eins = Score(s[i-1], t[j-1]) + matrix[i-1][j-1];
	/* gap in waagerechter Sequenz */
	int zwei = Score(s[i-1], '_') + matrix[i-1][j];
	/* gap in senkrechter Sequenz */
	int drei = Score('_', t[j-1]) + matrix[i][j-1];

	/* Vergleich */
	if (eins > zwei)
		maxi = eins;
	else
		maxi = zwei;
	if (drei > maxi)
		maxi = drei;

	return maxi;
}


/* Konstrukteur */
alignment::alignment (string sx, string tx) {
	s_aligned.append ("");
	t_aligned.append ("");

	/* Strings speichern */
	s.append (sx);
	t.append (tx);
	/* Laengen der Sequenzen = Groesse der Matrix */
	zeilen = sx.size () + 1;
	spalten = tx.size () + 1;

	/* Koordinate rechts unten = Score */
	l.x = zeilen - 1;
	l.y = spalten - 1;

	/* dynamische Matrix erstellen, Array aus Pointern auf Arrays */
	matrix = new int*[zeilen];
	/* Arrays, auf die die Pointer des 1. Arrays zeigen */
	for (int i = 0; i < zeilen; i++ )
		matrix[i] = new int[spalten];

	/* Randwerte der Matrix: 1. Zeile */
	for (int j = 1; j < spalten; j++)
		matrix[0][j] = -j;

	/* 1. Spalte */
	for (int i = 0; i < zeilen; i++)
		matrix[i][0] = -i;

	/* Berechnung der inneren Werte der Matrix */
	for (int i = 1; i < zeilen; i++)
		for (int j = 1; j < spalten; j++)
			matrix[i][j] = Max(i, j);

	/* Score setzen */
	score = matrix[l.x][l.y];
}

/* gibt den score aus */
int alignment::getScore () {
	return score;
}

/* berechnet alle moeglichen alignierten Sequenzen */
void alignment::constructAlignedSeq () {
	GetAlignment (l);
}

/* gibt alle alignierten Sequenzen zurueck */
void alignment::getAlignedSeq (vector<string> &alignedSeq) {
	for (unsigned int i = 0; i < alignments.size (); i++)  {
		alignedSeq.push_back (alignments[i]);
	}
}


/* Dekonstrukteur */
alignment::~alignment () {
	/* erst die einlzelnen Zeilen loeschen */
	for (int i = 0; i < zeilen; i++) {
		if (matrix[i] != NULL)
			delete[] matrix[i];
	}
	/* dann den Vektor, der auf die Zeilen gezeigt hat */
	delete[] matrix;
	/* Vektor loeschen */
	alignments.clear ();
}
