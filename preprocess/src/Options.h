
#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>

struct OPTS {
	std::string model;
	std::string ref;
	std::string json;
	std::string pdb;
	double dmax;
	int topn;
	int verb;
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);

#endif /* OPTIONS_H_ */
