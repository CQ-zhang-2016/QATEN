#include "Chain.h"
#include "AtomRecord.h"
#include "Info.h"
#include "Topology.h"

#include <cstring>
#include <cassert>
#include <cmath>
#include <algorithm>

using namespace std;

Chain::Chain() :
		kd(NULL), nRes(0), nAtoms(0), residues(0), atoms(0) {

	/* */

}

Chain::Chain(const string &name) :
		kd(NULL), nRes(0), nAtoms(0), residues(0), atoms(0) {

	// open file
	FILE *F;
	F = fopen(name.c_str(), "r");
	if (F == NULL) {
		fprintf(stderr, "ERROR: Can't open PDB file for reading '%s'\n",
			name.c_str());
		exit(1);
	}

	// read file line by line
	char buf[81];
	vector<vector<AtomRecord::Atom> > aRecVec;
	while (fgets(buf, 81, F) != NULL) {

		// skip all records except ATOM
		if (!((buf[0] == 'A') && (buf[1] == 'T'))) {
			continue;
		}

		// split into residues on the fly
		AtomRecord::Atom A = AtomRecord::ReadAtomRecord(buf);

		if (aRecVec.size() == 0) {
			// the first atom
			aRecVec.push_back(vector<AtomRecord::Atom>(1, A));
		} else {
			AtomRecord::Atom &B = aRecVec.back().back();
			if (B.resNum == A.resNum && B.insCode == A.insCode) {
				// same residue
				aRecVec.back().push_back(A);
			} else {
				// new residue
				aRecVec.push_back(vector<AtomRecord::Atom>(1, A));
			}
		}

	}
	fclose(F);

	// keep standard AAs only
	for (auto &res_vec : aRecVec) {
		Residue R(res_vec);
		if (R.type < 0 || R.type > 19) {
			continue;
		}
		if (R.atom["N"] == NULL || R.atom["CA"] == NULL || R.atom["C"] == NULL) {
			continue;
		}
		residues.push_back(R);
	}

	// fill in missing info
	Synchronize();
	nRes = residues.size();
	//nAtoms = atoms.size();

	SetKD();

}

Chain::Chain(const Chain &source) :
		kd(NULL), nRes(source.nRes), nAtoms(source.nAtoms), residues(source.residues), atoms(0) {

	Synchronize();
	SetKD();

}

Chain::Chain(const vector<Residue> &rvec) :
		kd(NULL), nRes(rvec.size()), nAtoms(0), residues(rvec), atoms(0) {

	Synchronize();
	SetKD();

}

Chain::~Chain() {

	if (kd != NULL) {
		kd_free(kd);
	}

}

Chain & Chain::operator=(const Chain &source) {

	assert(this != &source);

	residues = source.residues;
	nRes = source.nRes;
	nAtoms = source.nAtoms;

	Synchronize();
	SetKD();

	return *this;

}

void Chain::SetKD() {

	if (kd != NULL) {
		kd_free(kd);
	}

	kd = kd_create(3);

	for (auto A : atoms) {
		kd_insert3(kd, A->x, A->y, A->z, A);
	}

}

void Chain::Synchronize() {

	atoms.clear();
	resmap.clear();

	// 1) set atoms[] vector
	// 2) place correct pointers for Residue->Chain
	// 3) reset residue map
	nAtoms = 0;
	for (auto &R : residues) {
		nAtoms += R.nAtoms;
		R.chain = this;
		for (auto &A : R.atoms) {
			atoms.push_back(&A);
		}
		resmap[ { R.seqNum, R.insCode } ] = &R;
	}

}

void Chain::Save(const string &name) const {

	FILE *F = fopen(name.c_str(), "w");
	if (F == NULL) {
		fprintf(stderr, "ERROR: Cannot write to file '%s'\n", name.c_str());
		return;
	}

	for (auto A : atoms) {
		if (strlen(A->name) == 4 || A->name[0] == '1' || A->name[0] == '2'
			|| A->name[0] == '3' || A->name[0] == '4') {
			fprintf(F, "ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
				A->atomNum, A->name, A->altLoc, A->residue->name,
				A->residue->chainId, A->residue->seqNum,
				A->residue->insCode, A->x, A->y, A->z, A->occup, A->temp,
				A->element, A->charge);
		} else {
			fprintf(F, "ATOM  %5d  %-3s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
				A->atomNum, A->name, A->altLoc, A->residue->name,
				A->residue->chainId, A->residue->seqNum,
				A->residue->insCode, A->x, A->y, A->z, A->occup, A->temp,
				A->element, A->charge);
		}
	}
	fprintf(F, "TER\n");
	fclose(F);

}

string Chain::GetSequence() const {

	string seq(nRes, 'X');

	for (unsigned i = 0; i < nRes; i++) {
		int t = residues[i].type;
		if (t >= 0 && t < 20) {
			seq[i] = AAA1[t];
		}
	}

	return seq;

}

Residue* Chain::GetResidue(int n, char ins) const {

	map<pair<int, char>, Residue*>::const_iterator it;
	it = resmap.find( { n, ins });
	if (it == resmap.end()) {
		return NULL;
	} else {
		return it->second;
	}

}

Atom* Chain::GetAtom(int n, string aname) const {

	Residue *R = GetResidue(n);
	if (R == NULL) {
		//printf("Residue %d not found\n", n);
		return NULL;
	}

	return R->atom[aname];

}

Atom* Chain::GetAtom(const Atom &A) const {

	return GetAtom(A.residue->seqNum, A.name);

}

vector<pair<int, int> > Chain::GetBonds() const {

	vector<pair<int, int> > bonds;

	// intra-residue bonds
	for (auto &R : residues) {
		const set<pair<string, string> > &BSET = topology::CANONICAL20_BONDS.find(R.name)->second;
		for (const auto &ab : BSET) {
			Atom *A = R.atom.find(ab.first)->second;
			Atom *B = R.atom.find(ab.second)->second;
			if (A != NULL && B != NULL) {
				bonds.push_back( { A->atomNum, B->atomNum } );
			}
		}
	
	}

	// bonds between adjacent residues
	for (unsigned i = 1; i < nRes; i++) {
		Atom *A = residues[i - 1].atom.find("C")->second;
		Atom *B = residues[i].atom.find("N")->second;
		if (Atom::Dist(*A, *B) < 4.5) {
			bonds.push_back( { A->atomNum, B->atomNum } );
		}
	}

	return bonds;

}

vector<tuple<int, int, double> > Chain::GetContacts(double dmax) const {

	vector<tuple<int, int, double> > conts;

	for (auto A : atoms) {

		kdres *res = kd_nearest_range3(kd, A->x, A->y, A->z, dmax);
		double pos[3];
		while (!kd_res_end(res)) {
			Atom *B = (Atom*) kd_res_item(res, pos);
			if (B->atomNum > A->atomNum) {
				double d = Atom::Dist(*A, *B);
				conts.push_back(make_tuple(A->atomNum, B->atomNum, d));
			}
			kd_res_next(res);
		}
		kd_res_free(res);
	} 

	return conts;

}


vector<tuple<int, int, double> > Chain::GetContacts(double dmax, unsigned top) const {

	vector<tuple<int, int, double> > conts;

	for (auto &A : atoms) {

		vector<tuple<int, int, double> > contsa;

		kdres *res = kd_nearest_range3(kd, A->x, A->y, A->z, dmax);
		double pos[3];
		while (!kd_res_end(res)) {
			Atom *B = (Atom*) kd_res_item(res, pos);
			if (B->atomNum > A->atomNum) {
				double d = Atom::Dist(*A, *B);
				contsa.push_back(make_tuple(A->atomNum, B->atomNum, d));
			}
			kd_res_next(res);
		}
		kd_res_free(res);

		std::sort(std::begin(contsa), std::end(contsa),
			[](tuple<int, int, double> a, tuple<int, int, double> b)
			{return get<2>(a) < get<2>(b);});

		if (contsa.size() > top) {
			contsa = vector<tuple<int, int, double> >(contsa.begin(), contsa.begin() + top);
		}

		if (contsa.size() > 0) {
			conts.insert(end(conts), begin(contsa), end(contsa));
		}

	}

	return conts;

}

/*
vector<tuple<int, int, double> > Chain::GetContacts(double dmax, unsigned top) const {

	vector<tuple<int, int, double> > conts = GetContacts(dmax);

	std::sort(std::begin(conts), std::end(conts),
		[](tuple<int, int, double> a, tuple<int, int, double> b)
		{return get<2>(a) < get<2>(b);});

	if (conts.size() > top) {
		conts = vector<tuple<int, int, double> >(conts.begin(), conts.begin() + top);
	}

	return conts;

}
*/

int Chain::MaxAtoms() const {

	int n = 0;
	for (auto &R : residues) {
		n += topology::CANONICAL20_ATOMS.find(R.name)->second.size();
	}

	return n;

}


int Chain::MaxBonds() const {

	int n = 0;
	for (auto &R : residues) {
		n += topology::CANONICAL20_BONDS.find(R.name)->second.size();
	}

	n += nRes - 1; // inter-residue bonds

	return n;

}

int Chain::SetAtomsLFR() {

	// set default local frames first
	for (auto A : atoms) {
		A->SetDefaultLFR();
	}

	// fix backbone N,C
	for (unsigned i = 0; i < residues.size(); i++) {
		Residue &R = residues[i];
		Atom *N = R.atom["N"];
		Atom *CA = R.atom["CA"];
		Atom *C = R.atom["C"];
		Atom *O = R.atom["O"];
		if (i > 0) {
			Atom *Cprev = residues[i-1].atom["C"];
			N->SetLFR(Cprev, N, CA);
		}
		if (i < residues.size() - 1) {
			Atom *Nnext = residues[i+1].atom["N"];
			C->SetLFR(CA, C, Nnext);
			O->SetLFR(CA, O, Nnext);
		}
	}

	return 0;

}
