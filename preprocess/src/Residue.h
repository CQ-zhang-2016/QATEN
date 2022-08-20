#ifndef RESIDUE_H_
#define RESIDUE_H_

#include "Atom.h"
#include "AtomRecord.h"

#include <vector>
#include <map>

class Chain;

using namespace std;

class Residue {
private:

	int Standardize();
	int SwapAtoms(const pair<string, string> &flippair, int doit);
	int SetAtomMap();

public:

	char type;
	char name[4];
	int seqNum;
	char insCode;
	char chainId;

	unsigned nAtoms;
	map<string, Atom*> atom;
	vector<Atom> atoms;

	Chain *chain;

	Residue();
	Residue(const Residue &source);
	Residue(const vector<AtomRecord::Atom> &aRec);
	~Residue();

	Residue& operator=(const Residue &source);

	vector<pair<string, string> > Flip(int doit = 0);

	//static double MinDist();

};

#endif /* RESIDUE_H_ */
