#include "Options.h"
#include "Chain.h"
#include "MyLDDT.h"

#include <vector>
#include <utility>

int main(int argc, char *argv[]) {

	//
	// (1) Initializations
	//
	OPTS opts = { "", "", "", "", 999.99, 100, 1 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	Chain Model(opts.model);

	Chain *Ref = NULL;
	if (opts.ref != "") {
		Ref = new Chain;
		*Ref = Chain(opts.ref);
	}

	Model.SetAtomsLFR();

	// make atom numbering in the Model sequential
	// and starting from zero
	// !!! very important !!!
	for (unsigned i = 0; i < Model.nAtoms; i++) {
		Atom *A = Model.atoms[i];
		A->atomNum = i;
		//A->SetDefaultLFR();
		/*
		// X
		printf("ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
			A->atomNum*2, "C", ' ', "XXX", 'A', A->residue->seqNum, ' ', A->x, A->y, A->z, 1.0, A->temp, "", "");
		printf("ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
			A->atomNum*2+1, "C", ' ', "XXX", 'A', A->residue->seqNum, ' ', A->x + A->lfr[0][0], A->y + A->lfr[0][1], A->z + A->lfr[0][2], 1.0, A->temp, "", "");
		printf("TER\nEND\n");
		// Y
		printf("ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
			A->atomNum*2, "C", ' ', "YYY", 'A', A->residue->seqNum, ' ', A->x, A->y, A->z, 1.0, A->temp, "", "");
		printf("ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
			A->atomNum*2+1, "C", ' ', "YYY", 'A', A->residue->seqNum, ' ', A->x + A->lfr[1][0], A->y + A->lfr[1][1], A->z + A->lfr[1][2], 1.0, A->temp, "", "");
		printf("TER\nEND\n");
		// Z
		printf("ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
			A->atomNum*2, "C", ' ', "ZZZ", 'A', A->residue->seqNum, ' ', A->x, A->y, A->z, 1.0, A->temp, "", "");
		printf("ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
			A->atomNum*2+1, "C", ' ', "ZZZ", 'A', A->residue->seqNum, ' ', A->x + A->lfr[2][0], A->y + A->lfr[2][1], A->z + A->lfr[2][2], 1.0, A->temp, "", "");
		printf("TER\nEND\n");
		*/
	}

	// additional data on the Model
	std::vector<std::pair<int, int> > bonds = Model.GetBonds();
	std::vector<std::tuple<int, int, double> > contacts = Model.GetContacts(opts.dmax, opts.topn);

/*
	if (opts.verb > 0) {
		printf("# %s\n", std::string(78, '-').c_str());
		printf("# %15s : %s\n", "model", opts.model.c_str());
		printf("# %s\n", std::string(78, '-').c_str());
		printf("# %15s   %9s %9s\n", "", "all", "missing");
		printf("# %15s : %9d %9d\n", "atoms", Model.nAtoms, Model.MaxAtoms() - Model.nAtoms);
		printf("# %15s : %9lu %9lu\n", "bonds", bonds.size(), Model.MaxBonds() - bonds.size());
		printf("# %15s : %9lu %9s   (%.1lfA)\n", "contacts", contacts.size(), "-", opts.dmax);
	}
*/
	//
	// (2) Run LDDT if Ref is set
	//
	MyLDDT::Result lddt_ts, lddt_ha;
	if (Ref != NULL) {
		MyLDDT lddt(Model, *Ref, opts.verb);
		lddt_ha = lddt.GetScore(MyLDDT::HA);
		lddt_ts = lddt.GetScore(MyLDDT::TS);
	}

	//
	// (3) Save results in JSON
	//
	if (opts.json != "") {

		FILE *F = fopen(opts.json.c_str(), "w");
		if (F == NULL) {
			fprintf(stderr, "ERROR: cannot open '%s' file for writing\n", 
				opts.json.c_str());
		}

		// atoms
		fprintf(F, "{\n\"atoms\":[");
		char buf[20];
		for (unsigned i = 0; i < Model.nAtoms; i++) {
			Atom *A = Model.atoms[i];
			sprintf(buf, "%s_%s", A->residue->name, A->name);
			if (i < Model.nAtoms - 1) {
				fprintf(F, "\"%s\",", buf);
			}
		}
		fprintf(F, "\"%s\"],\n", buf);

		// residue indices
		fprintf(F, "\"res_idx\":[");
		for (unsigned i = 0; i < Model.residues.size(); i++) {
			Residue &R = Model.residues[i];
			for (unsigned j = 0; j < R.atoms.size(); j++) {
				if (i == Model.residues.size() - 1 && j == R.atoms.size() - 1) {
					fprintf(F, "%d],\n", i);
				} else {
					fprintf(F, "%d,", i);
				}
			}
		}

		// bonds
		fprintf(F, "\"bonds\":[");
		for (unsigned i = 0; i < bonds.size() - 1; i++) {
			fprintf(F, "[%d,%d],", bonds[i].first, bonds[i].second);
		}
		fprintf(F,"[%d,%d]],\n\"contacts\":[",bonds.back().first, bonds.back().second);

		// contacts in the model
		for (unsigned i = 0; i < contacts.size(); i++) {
			int a = std::get<0>(contacts[i]);
			int b = std::get<1>(contacts[i]);
			double d = std::get<2>(contacts[i]);
			double xyz_ab[3], xyz_ba[3];
			Model.atoms[a]->Project(Model.atoms[b], xyz_ab);
			Model.atoms[b]->Project(Model.atoms[a], xyz_ba);
			if (i < contacts.size() - 1) {
				fprintf(F, "[%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf],",
					a, b, d, xyz_ab[0], xyz_ab[1], xyz_ab[2], xyz_ba[0], xyz_ba[1], xyz_ba[2]);
			} else {
				fprintf(F, "[%d,%d,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf]]",
					a, b, d, xyz_ab[0], xyz_ab[1], xyz_ab[2], xyz_ba[0], xyz_ba[1], xyz_ba[2]);
			}
		}

		// contacts in the reference
		if (Ref != NULL) {
			fprintf(F,",\n\"ref_dist\":[");
			for (unsigned i = 0; i < contacts.size(); i++) {
				int a = std::get<0>(contacts[i]);
				int b = std::get<1>(contacts[i]);
				Atom *A = Ref->GetAtom(*(Model.atoms[a]));
				Atom *B = Ref->GetAtom(*(Model.atoms[b]));
				if (A != NULL && B != NULL) {
					if (i < contacts.size() - 1) {
						fprintf(F,"%.4lf,", Atom::Dist(*A, *B));
					} else {
						fprintf(F,"%.4lf]", Atom::Dist(*A, *B));
					}
				} else {
					if (i < contacts.size() - 1) {
						fprintf(F,"nan,");
					} else {
						fprintf(F,"nan]");
					}
				}
			}
		}

		// scores
		if (Ref != NULL ) {

			// residue-level
			fprintf(F, ",\n\"rscores\":[");
			for (unsigned i = 0; i < lddt_ha.rscore.size(); i++) {
				double ha = lddt_ha.rscore[i];
				double ts = lddt_ts.rscore[i];
				if (i != lddt_ha.rscore.size() - 1) {
					fprintf(F, "[%.5lf,%.5lf],", ha, ts);
				} else {
					fprintf(F, "[%.5lf,%.5lf]]", ha, ts);
				}
			}

			// atom-level
			fprintf(F, ",\n\"ascores\":[");
			for (unsigned i = 0; i < lddt_ha.ascore.size(); i++) {
				double ha = lddt_ha.ascore[i];
				double ts = lddt_ts.ascore[i];
				if (i != lddt_ha.ascore.size() - 1) {
					fprintf(F, "[%.5lf,%.5lf],", ha, ts);
				} else {
					fprintf(F, "[%.5lf,%.5lf]],\n", ha, ts);
				}
			}

			// full structure
			fprintf(F, "\"lddt\":[%.5f,%.5f]\n}\n", lddt_ha.score, lddt_ts.score);
		} else {
			fprintf(F, "\n}\n");
		}

		fclose(F);
	}

	//
	// (4) Save cleaned PDB of the Model
	//
	if (opts.pdb != "") {

		// save LDDT atomic scores as B-factors
		if (Ref != NULL) {
			for (auto A : Model.atoms) {
				A->temp = lddt_ts.ascore[A->atomNum];
			}
		}
		Model.Save(opts.pdb);
	}

	//
	// clean up
	//
	if (Ref != NULL) {
		delete Ref;
	}

	return 0;

}

