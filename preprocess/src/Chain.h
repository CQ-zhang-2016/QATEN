#ifndef CHAIN_H_
#define CHAIN_H_

#include <vector>
#include <string>
#include <utility>
#include <map>

#include "AtomRecord.h"
#include "Atom.h"
#include "Residue.h"
#include "kdtree.h"

using namespace std;

class Chain {

friend class MyLDDT;

private:

	kdtree *kd;

	map<pair<int, char>, Residue*> resmap; // {seqNum,insCode} --> Residue

public:

	unsigned nRes;
	unsigned nAtoms;
	vector<Residue> residues;
	vector<Atom*> atoms;

	Chain();
	Chain(const string &name);
	Chain(const Chain &source);
	Chain(const vector<Residue> &rvec);
	//Chain(const vector<string> &atoms_str);
	~Chain();

	Chain& operator=(const Chain &source);

	void Synchronize();
	void SetKD();

	void Save(const string &name) const;

	Residue* GetResidue(int n, char ins = ' ') const;
	Atom* GetAtom(int n, string name) const;
	Atom* GetAtom(const Atom &A) const;

	vector<pair<int, int> > GetBonds() const;
	vector<tuple<int, int, double> > GetContacts(double dmax) const;
	vector<tuple<int, int, double> > GetContacts(double dmax, unsigned top) const;

	int MaxBonds() const;
	int MaxAtoms() const;

	string GetSequence() const;

	// TODO: implement
	int SetAtomsLFR();

};

#endif /* CHAIN_H_ */
