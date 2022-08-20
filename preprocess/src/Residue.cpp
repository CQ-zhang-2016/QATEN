#include "Residue.h"
#include "Info.h"
#include "Topology.h"

#include <cstring>
#include <cassert>

using namespace std;

Residue::Residue() :
		type(23), seqNum(-1), insCode(' '), chainId(' '), nAtoms(0), chain(NULL) {

	strcpy(name, "UNK");

}

Residue::Residue(const Residue &source) :
		type(source.type), seqNum(source.seqNum), insCode(source.insCode), 
		chainId(source.chainId), nAtoms(source.nAtoms), atom(source.atom), 
		atoms(source.atoms), chain(source.chain) {

	strcpy(name, source.name);

	for (auto &A : atoms) {
		atom[A.name] = &A;
		A.residue = this;
	}

}

Residue::Residue(const vector<AtomRecord::Atom> &aRec) :
		seqNum(aRec.front().resNum), insCode(aRec.front().insCode),
		chainId(aRec.front().chainId), chain(NULL) {

	strcpy(name, aRec.front().resName);

	// check for alternative conformations
	int alt_count = 0;
	for (auto &A : aRec) {
		alt_count += (A.altLoc != ' ');
	}

	if (alt_count == 0) { // no alternative conformations
		for (auto &A : aRec) {
			atoms.push_back(A);
		}
	} else { // there are alternative conformations

		// save top-populated conformation for each atom
                map<string, const AtomRecord::Atom*> aRec_map;
                map<string, const AtomRecord::Atom*>::iterator it;
                const AtomRecord::Atom &A0 = aRec.front();
                aRec_map[A0.atomName] = &A0;

		for (const auto &A : aRec) {

			// save only atoms with the same residue name as the first atom
			if (strcmp(A0.resName, A.resName) != 0) {
				continue;
			}

			// save highest occupancy atoms
			it = aRec_map.find(A.atomName);
			if (it == aRec_map.end()) {
				aRec_map[A.atomName] = &A;
			} else {
				if (A.occup > it->second->occup) {
					it->second = &A;
				}
			}

		}

		// keep selected atoms
		for (auto &A : aRec_map) {
			atoms.push_back(*(A.second));
                }
	}

	Standardize();

	nAtoms = atoms.size();

}

Residue::~Residue() {

	/* */

}

Residue & Residue::operator =(const Residue & source) {

	assert(this != &source);

	type = source.type;
	strcpy(name, source.name);
	seqNum = source.seqNum;
	insCode = source.insCode;
	nAtoms = source.nAtoms;
	chainId = source.chainId;

	atom = source.atom;
	atoms = source.atoms;

	chain = source.chain;

	for (auto &A : atoms) {
		atom[A.name] = &A;
		A.residue = this;
	}

	return *this;

}

int Residue::Standardize() {

	type = -1;

	// if type is non-standard try conversion
	map<string, char>::const_iterator cit;
	cit = topology::CANONICAL20_MAP.find(name);
	if (cit == topology::CANONICAL20_MAP.end()) {
		return -1;
	}

	for (int i = 0; i < 20; i++) {
		if (cit->second == AAA1[i]) {
			strcpy(name, AAA3[i]);
			type = i;
		}
	}

	if (type < 0 || type >= 20) {
		return -1;
	}

	vector<Atom> atoms_new;
	const set<string> &ASET = topology::CANONICAL20_ATOMS.find(name)->second;

	// save existing atoms
	for (auto &A : atoms) {
		if (ASET.find(A.name) != ASET.end()) {
			atoms_new.push_back(A);
		}
	}
	atoms = atoms_new;

/*
	// save map for existing atoms
	atom.clear();
	for (auto &A : atoms) {
		atom[A.name] = &A;
		A.residue = this;
	}

	// fill in missing atoms
	int missing = 0;
	for (const auto &aname : ASET) {
		if (atom.find(aname) == atom.end()) {
			atom[aname] = NULL;
			missing++;
		}
	}
*/

	int missing = SetAtomMap();

	return missing;

}

int Residue::SetAtomMap() {

	// save map for existing atoms
	atom.clear();
	for (auto &A : atoms) {
		atom[A.name] = &A;
		A.residue = this;
	}

	const set<string> &ASET = topology::CANONICAL20_ATOMS.find(name)->second;

	// fill in missing atoms
	int missing = 0;
	for (const auto &aname : ASET) {
		if (atom.find(aname) == atom.end()) {
			atom[aname] = NULL;
			missing++;
		}
	}

	return missing;

}

vector<pair<string, string> > Residue::Flip(int doit) {

	vector<pair<string, string> > flipvec, flipvec_valid;

	if (string(name) == "GLU") {
		flipvec.push_back( { "OE1", "OE2" } );
	} else if (string(name) == "ASP") {
		flipvec.push_back( { "OD1", "OD2" } );
	} else if (string(name) == "PHE" || string(name) == "TYR") {
		flipvec.push_back( { "CD1", "CD2" } );
		flipvec.push_back( { "CE1", "CE2" } );
	} else if (string(name) == "ARG") {
		flipvec.push_back( { "NH1", "NH2" } );
	}

	for (auto &flippair : flipvec) {
		if (SwapAtoms(flippair, doit) > 0) {
			flipvec_valid.push_back(flippair);
		}
	}

	// update pointers, if needed
	if (doit > 0 && flipvec_valid.size() > 0) {
		SetAtomMap();
	}

	return flipvec_valid;

}

int Residue::SwapAtoms(const pair<string, string> &flippair, int doit) {

	map<string, Atom*>::iterator it;

	// is first atom present?
	it = atom.find(flippair.first);
	if (it == atom.end() || it->second == NULL) {
		return 0;
	}
	Atom *A1 = it->second;

	// is second atom present?
	it = atom.find(flippair.second);
	if (it == atom.end() || it->second == NULL) {
		return 0;
	}
	Atom *A2 = it->second;

	// swap atoms
	if (doit > 0) {
		Atom BUF = *A1;
		*A1 = *A2;
		*A2 = BUF;
		strcpy(A1->name, flippair.first.c_str());
		strcpy(A2->name, flippair.second.c_str());
	}

	return 1;

}
