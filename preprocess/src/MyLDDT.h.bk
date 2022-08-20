#ifndef MYLDDT_H_
#define MYLDDT_H_

#include "Chain.h"

#include <map>
//#include <unordered_map>

using namespace std;

class MyLDDT {

private:

	int verbosity;
	double dmax;
	std::vector<double> dbins;
	string mlabel;

	const Chain &Mod;
	const Chain &Ref;

	// substructures which can be compared
	Chain Ref_sub;
	Chain Mod_sub;

	// substructures without flippable residues
	Chain Ref_sub_fixed;
	Chain Mod_sub_fixed;

	// leave undefined
	MyLDDT();
	MyLDDT(const MyLDDT&);

	double lddt(const Atom &Amod, const Chain &Mod_, const Atom &Aref, const Chain &Ref_);

public:

	enum MODE {
		HA = 0, // high-res
		TS = 1  // low-res
	};

	struct Result {
		double score;
		vector<double> ascore;
		vector<double> rscore;
		vector<int> amask;
	};

	MyLDDT(const Chain &Model, const Chain &Reference, int verb = 1);
	~MyLDDT();

	Result GetScore(MODE mode = TS);

};

#endif // MYLDDT_H_
