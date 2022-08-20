#include "MyLDDT.h"

#include <cstring>
#include <cmath>
#include <cassert>

MyLDDT::MyLDDT(const Chain &Model, const Chain &Reference, int v) :
		verbosity(v), dmax(15.0), dbins(0), mlabel(""), Mod(Model), Ref(Reference) {

	// select reference sesidue set
	// present in both structures
	vector<Residue> Mod_vec, Ref_vec;
	int mismatch = 0;
	for (auto &ref_res : Ref.residues) {
		Residue *mod_res = Mod.GetResidue(ref_res.seqNum);

		// residue in the Reference
		// is missing from the Model
		if (mod_res == NULL /* || strcmp(ref_res.name, mod_res->name) != 0 */ ) {
			continue;
		}
		mismatch += strcmp(ref_res.name, mod_res->name) != 0;

		// save otherwise
		Ref_vec.push_back(ref_res);
		Mod_vec.push_back(*mod_res);

	}

	Ref_sub = Chain(Ref_vec);
	Mod_sub = Chain(Mod_vec);

	// save unflippable core
	Ref_vec.clear();
	Mod_vec.clear();
	for (unsigned i = 0; i < Ref_sub.nRes; i++) {
		Residue &ref_res = Ref_sub.residues[i];
		Residue &mod_res = Mod_sub.residues[i];
		if (mod_res.Flip(0).size() == 0) {
			Ref_vec.push_back(ref_res);
			Mod_vec.push_back(mod_res);
		}
	}

	Ref_sub_fixed = Chain(Ref_vec);
	Mod_sub_fixed = Chain(Mod_vec);

	if (verbosity > 0) {
		printf("# %s\n", string(78, '-').c_str());
		printf("# Local distance test (LDDT) - unofficial %38s\n", "v0.1");
		printf("# %s\n", string(78, '-').c_str());
		printf("# Number of residues:\n");
		printf("# %15s : %u\n", "model", Mod.nRes);
		printf("# %15s : %u\n", "reference", Ref.nRes);
		printf("# %15s : %u\n", "shared", Mod_sub.nRes);
		printf("# %15s : %d\n", "mismatches", mismatch);
		printf("# %15s : %u\n", "flippable", Mod_sub.nRes - Mod_sub_fixed.nRes);
		printf("# %s\n", string(78, '-').c_str());
	}

}

MyLDDT::~MyLDDT() {

	/* */

}

MyLDDT::Result MyLDDT::GetScore(MODE mode) {

	// set distance bins depending on the mode
	switch (mode) {
		case HA:
			dbins = { 0.5, 1.0, 2.0, 4.0 };
			mlabel = "HA";
			break;
		case TS:
			dbins = { 1.0, 2.0, 4.0, 8.0 };
			mlabel = "TS";
			break;
		default:
			fprintf(stderr, "ERROR: invalid LDDT mode %d\n", mode);
			break;
	}

	// probe every flippabel residue, one at a time:
	// if the score calculated over all unflippable residues
	// improves upon flipping, then flip!
	Chain Mod_sub_flip(Mod_sub); // temp Chain for flipping
	if (verbosity > 1) {
		printf("# Flipped residues (%s):\n", mlabel.c_str());
		printf("# %5s %5s %8s %8s %6s\n", "num", "res", "orig", "flipped", "%imp");
	}
	int nflipped = 0;
	for (auto &mod_res : Mod_sub_flip.residues) {

		double lddt_orig = 0.0;
		double lddt_flip = 0.0;

		bool fl = false;
		// first, check but not flip
		for (auto &flippair : mod_res.Flip(0)) {

			fl = true;

			Atom &A1 = *(mod_res.atom[flippair.first]);
			Atom &A2 = *(mod_res.atom[flippair.second]);

			// skip if atom(s) are missing in the Reference
			Atom *A1ref = Ref_sub.GetAtom(A1);
			Atom *A2ref = Ref_sub.GetAtom(A2);
			if (A1ref == NULL || A2ref == NULL) {
				continue;
			}

			// score original placement
			lddt_orig += lddt(A1, Mod_sub_fixed, *A1ref, Ref_sub_fixed);
			lddt_orig += lddt(A2, Mod_sub_fixed, *A2ref, Ref_sub_fixed);

			// score flipped atoms
			lddt_flip += lddt(A1, Mod_sub_fixed, *A2ref, Ref_sub_fixed);
			lddt_flip += lddt(A2, Mod_sub_fixed, *A1ref, Ref_sub_fixed);
		}

		// flip if score improves
		if (fl && lddt_flip - lddt_orig > 1e-3) {
			mod_res.Flip(1);
			nflipped++;
			if (verbosity > 1) {
				printf("# %5d %5s %8.5lf %8.5lf %5.1lf%%\n",
					mod_res.seqNum, mod_res.name, lddt_orig, lddt_flip,
					(lddt_flip - lddt_orig) / lddt_orig * 100.0);
			}
		}

	}
	if (verbosity > 1 && nflipped == 0) {
		printf("# (none)\n");
	}

	// update Mod_sub_flip
	if (nflipped > 0) {
		Mod_sub_flip.Synchronize();
		Mod_sub_flip.SetKD();
	}

	//
	// calculate final scores
	//
	Result result = {
		0.0,
		vector<double>(Mod.nAtoms, 0.0),
		vector<double>(Mod.nRes, 0.0),
		vector<int>(Mod.nAtoms, 0)
	};

	for (auto Amod : Mod_sub_flip.atoms) {
		Atom *Aref = Ref_sub.GetAtom(*Amod);
		double s = 0.0;
		if (Aref != NULL) {
			s = lddt(*Amod, Mod_sub_flip, *Aref, Ref);
			result.amask[Amod->atomNum] = 1;
		}
		result.ascore[Amod->atomNum] = s;
		result.score += s;
	}

	// normalize by atoms in Reference:
	// atoms missing in Model contribute zero
	// to the total score
	result.score /= Ref.nAtoms;

	// aggregate atomic scores into residue scores
	for (unsigned i = 0; i < Mod.nRes; i++) {
		const Residue &R = Mod.residues[i];
		double &s = result.rscore[i];
		double n = 0;
		for (auto &A : R.atoms) {
			s += result.ascore[A.atomNum];
			n += result.amask[A.atomNum];
		}
		if (n > 0) {
			s /= n;
		}
	}

	if (verbosity > 2) {
		printf("# %s\n", string(78, '-').c_str());
		for (unsigned i = 0; i < Mod.nRes; i++) {
			const Residue &R = Mod.residues[i];
			double s = result.rscore[i];
			int l = (int)round(s * 40);
			int n = 0;
			for (auto &A : R.atoms) {
				n += result.amask[A.atomNum];
			}
			string bar = l > 0 ? string(l, '*') : "";
			printf("# %2s %5d %5s %8.5lf   %2d/%-2lu   . %-40s .\n",
				mlabel.c_str(), R.seqNum, R.name, s,
				n, R.atom.size(), bar.c_str() /* string(l, '*').c_str() */ );
		}
		printf("# %s\n", string(78, '-').c_str());
	}

	if (verbosity > 0) {
		printf("# LDDT(%s)= %.5lf\n", mlabel.c_str(), result.score);
	}

	return result;

}

double MyLDDT::lddt(const Atom &Amod, const Chain &Mod_, const Atom &Aref, const Chain &Ref_) {

	// save neighbors and corresponding Mod/Ref distances
	// in this vector
	map<pair<int, string>, pair<double, double> > neighbors;

	// loop over contacts in the Model
	{
		kdres *res = kd_nearest_range3(Mod_.kd, Amod.x, Amod.y, Amod.z, dmax);
		double pos[3];
		while (!kd_res_end(res)) {

			Atom *Bmod = (Atom*) kd_res_item(res, pos);
			Atom *Bref = Ref_.GetAtom(*Bmod);

			// skip contact if:
			// 1) atoms are from same residue
			// 2) neighboring atom Bref is missing from Ref_
			if (Amod.residue->seqNum != Bmod->residue->seqNum && Bref != NULL) {
				double dmod = Atom::Dist(Amod, *Bmod);
				double dref = Atom::Dist(Aref, *Bref);
				pair<int, string> key = { Bmod->residue->seqNum, Bmod->name };
				neighbors[key] = { dmod, dref };
			}
			kd_res_next(res);
		}
		kd_res_free(res);
	}

	// loop over contacts in the Reference
	{
		kdres *res = kd_nearest_range3(Ref_.kd, Aref.x, Aref.y, Aref.z, dmax);
		double pos[3];
		while (!kd_res_end(res)) {

			Atom *Bref = (Atom*) kd_res_item(res, pos);

			pair<int, string> key = { Bref->residue->seqNum, Bref->name };

			// skip contact if:
			// 1) it is already among neighbors
			// 2) atoms are from same residue
			if (Aref.residue->seqNum != Bref->residue->seqNum && 
				neighbors.find(key) == neighbors.end() ) {

				Atom *Bmod = Mod_.GetAtom(*Bref);

				double dref = Atom::Dist(Aref, *Bref);
				double dmod = Bmod != NULL ? Atom::Dist(Amod, *Bmod) : 1e9;

				neighbors[key] = { dmod, dref };
			}
			kd_res_next(res);
		}
		kd_res_free(res);
	}

	// calculate final score
	double score = 0.0, smax = 0.0;
	for (auto &n : neighbors) {
		double dmin = n.second.first < n.second.second ? n.second.first : n.second.second;
		if (dmin >= dmax) {
			continue;
		}
		double delta = fabs(n.second.first - n.second.second);
		smax += 1.0;

		for (unsigned i = 0; i < dbins.size(); i++) {
			if (delta < dbins[i]) {
				score += 1.0;
			}
		}

	}
	if (smax > 1e-5) {
		score /= smax * dbins.size();
	} else {
		score = 0.0;
	}

	return score;

}
