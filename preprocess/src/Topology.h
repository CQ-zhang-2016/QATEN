/*
 * Topology.h
 *
 *  Created on: Oct 24, 2018
 *      Author: aivan
 */

#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

#include <map>
#include <set>
#include <utility>
#include <string>

using namespace std;

namespace topology {

const map<string, set<string> > CANONICAL20_ATOMS =
{
	{ "ALA", { "C", "CA", "CB", "N", "O" } },
	{ "ARG", { "C", "CA", "CB", "CD", "CG", "CZ", "N", "NE", "NH1", "NH2", "O" } },
	{ "ASN", { "C", "CA", "CB", "CG", "N", "ND2", "O", "OD1" } },
	{ "ASP", { "C", "CA", "CB", "CG", "N", "O", "OD1", "OD2" } },
	{ "CYS", { "C", "CA", "CB", "N", "O", "SG" } },
	{ "GLN", { "C", "CA", "CB", "CD", "CG", "N", "NE2", "O", "OE1" } },
	{ "GLU", { "C", "CA", "CB", "CD", "CG", "N", "O", "OE1", "OE2" } },
	{ "GLY", { "C", "CA", "N", "O" } },
	{ "HIS", { "C", "CA", "CB", "CD2", "CE1", "CG", "N", "ND1", "NE2", "O" } },
	{ "ILE", { "C", "CA", "CB", "CD1", "CG1", "CG2", "N", "O" } },
	{ "LEU", { "C", "CA", "CB", "CD1", "CD2", "CG", "N", "O" } },
	{ "LYS", { "C", "CA", "CB", "CD", "CE", "CG", "N", "NZ", "O" } },
	{ "MET", { "C", "CA", "CB", "CE", "CG", "N", "O", "SD" } },
	{ "PHE", { "C", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "CZ", "N", "O" } },
	{ "PRO", { "C", "CA", "CB", "CD", "CG", "N", "O" } },
	{ "SER", { "C", "CA", "CB", "N", "O", "OG" } },
	{ "THR", { "C", "CA", "CB", "CG2", "N", "O", "OG1" } },
	{ "TRP", { "C", "CA", "CB", "CD1", "CD2", "CE2", "CE3", "CG", "CH2", "CZ2", "CZ3", "N", "NE1", "O" } },
	{ "TYR", { "C", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "CZ", "N", "O", "OH" } },
	{ "VAL", { "C", "CA", "CB", "CG1", "CG2", "N", "O" } }
};

const map<pair<string, string>, tuple<string, string, string> > CANONICAL20_LFR = {

	// { "ALA", { "C", "CA", "CB", "N", "O" } }
	{ { "ALA", "CA" }, { "N", "CA", "C" } },
	{ { "ALA", "CB" }, { "CA", "N", "CB" } },
	{ { "ALA", "O" }, { "C", "CA", "O" } },

	//{ "ARG", { "C", "CA", "CB", "CD", "CG", "CZ", "N", "NE", "NH1", "NH2", "O" } },
	{ { "ARG", "CA" }, { "N", "CA", "C" } },
	{ { "ARG", "CB" }, { "CA", "CB", "CG" } },
	{ { "ARG", "CG" }, { "CB", "CG", "CD" } },
	{ { "ARG", "CD" }, { "CG", "CD", "NE" } },
	{ { "ARG", "NE" }, { "CD", "NE", "CZ" } },
	{ { "ARG", "CZ" }, { "NH1", "CZ", "NH2" } },
	{ { "ARG", "NH1" }, { "NE", "CZ", "NH1" } },
	{ { "ARG", "NH2" }, { "NE", "CZ", "NH2" } },
	{ { "ARG", "O" }, { "C", "CA", "O" } },

	// { "ASN", { "C", "CA", "CB", "CG", "N", "ND2", "O", "OD1" } },
	{ { "ASN", "CA" }, { "N", "CA", "C" } },
	{ { "ASN", "CB" }, { "CA", "CB", "CG" } },
	{ { "ASN", "CG" }, { "OD1", "CG", "ND2" } },
	{ { "ASN", "ND2" }, { "CG", "CB", "ND2" } },
	{ { "ASN", "OD1" }, { "CG", "CB", "OD1" } },
	{ { "ASN", "O" }, { "C", "CA", "O" } },

	// { "ASP", { "C", "CA", "CB", "CG", "N", "O", "OD1", "OD2" } },
	{ { "ASP", "CA" }, { "N", "CA", "C" } },
	{ { "ASP", "CB" }, { "CA", "CB", "CG" } },
	{ { "ASP", "CG" }, { "OD1", "CG", "OD2" } },
	{ { "ASP", "OD2" }, { "CG", "CB", "OD2" } },
	{ { "ASP", "OD1" }, { "CG", "CB", "OD1" } },
	{ { "ASP", "O" }, { "C", "CA", "O" } },

	// { "CYS", { "C", "CA", "CB", "N", "O", "SG" } },
	{ { "CYS", "CA" }, { "N", "CA", "C" } },
	{ { "CYS", "CB" }, { "CA", "CB", "SG" } },
	{ { "CYS", "SG" }, { "CB", "CA", "SG" } },
	{ { "CYS", "O" }, { "C", "CA", "O" } },

	// { "GLN", { "C", "CA", "CB", "CD", "CG", "N", "NE2", "O", "OE1" } },
	{ { "GLN", "CA" }, { "N", "CA", "C" } },
	{ { "GLN", "CB" }, { "CA", "CB", "CG" } },
	{ { "GLN", "CG" }, { "CB", "CG", "CD" } },
	{ { "GLN", "CD" }, { "OE1", "CD", "NE2" } },
	{ { "GLN", "OE1" }, { "CD", "CG", "OE1" } },
	{ { "GLN", "NE2" }, { "CD", "CG", "NE2" } },
	{ { "GLN", "O" }, { "C", "CA", "O" } },

	// { "GLU", { "C", "CA", "CB", "CD", "CG", "N", "O", "OE1", "OE2" } },
	{ { "GLU", "CA" }, { "N", "CA", "C" } },
	{ { "GLU", "CB" }, { "CA", "CB", "CG" } },
	{ { "GLU", "CG" }, { "CB", "CG", "CD" } },
	{ { "GLU", "CD" }, { "OE1", "CD", "OE2" } },
	{ { "GLU", "OE1" }, { "CD", "CG", "OE1" } },
	{ { "GLU", "OE2" }, { "CD", "CG", "OE2" } },
	{ { "GLU", "O" }, { "C", "CA", "O" } },

	// { "GLY", { "C", "CA", "N", "O" } },
	{ { "GLY", "CA" }, { "N", "CA", "C" } },
	{ { "GLY", "O" }, { "C", "CA", "O" } },

	// { "HIS", { "C", "CA", "CB", "CD2", "CE1", "CG", "N", "ND1", "NE2", "O" } },
	{ { "HIS", "CA" }, { "N", "CA", "C" } },
	{ { "HIS", "CB" }, { "CA", "CB", "CG" } },
	{ { "HIS", "CG" }, { "ND1", "CG", "CD2" } },
	{ { "HIS", "ND1" }, { "CG", "ND1", "CE1" } },
	{ { "HIS", "CD2" }, { "CG", "CD2", "NE2" } },
	{ { "HIS", "CE1" }, { "ND1", "CE1", "NE2" } },
	{ { "HIS", "NE2" }, { "CD2", "NE2", "CE1" } },
	{ { "HIS", "O" }, { "C", "CA", "O" } },

	// { "ILE", { "C", "CA", "CB", "CD1", "CG1", "CG2", "N", "O" } },
	{ { "ILE", "CA" }, { "N", "CA", "C" } },
	{ { "ILE", "CB" }, { "CA", "CB", "CG1" } },
	{ { "ILE", "CG1" }, { "CB", "CG1", "CD1" } },
	{ { "ILE", "CG2" }, { "CB", "CA", "CG2" } },
	{ { "ILE", "CD1" }, { "CG1", "CB", "CD1" } },
	{ { "ILE", "O" }, { "C", "CA", "O" } },

	// { "LEU", { "C", "CA", "CB", "CD1", "CD2", "CG", "N", "O" } },
	{ { "LEU", "CA" }, { "N", "CA", "C" } },
	{ { "LEU", "CB" }, { "CA", "CB", "CG" } },
	{ { "LEU", "CG" }, { "CD1", "CG", "CD2" } },
	{ { "LEU", "CD1" }, { "CG", "CB", "CD1" } },
	{ { "LEU", "CD2" }, { "CG", "CB", "CD2" } },
	{ { "LEU", "O" }, { "C", "CA", "O" } },

	// { "LYS", { "C", "CA", "CB", "CD", "CE", "CG", "N", "NZ", "O" } },
	{ { "LYS", "CA" }, { "N", "CA", "C" } },
	{ { "LYS", "CB" }, { "CA", "CB", "CG" } },
	{ { "LYS", "CG" }, { "CB", "CG", "CD" } },
	{ { "LYS", "CD" }, { "CG", "CD", "CE" } },
	{ { "LYS", "CE" }, { "CD", "CE", "NZ" } },
	{ { "LYS", "NZ" }, { "CE", "CD", "NZ" } },
	{ { "LYS", "O" }, { "C", "CA", "O" } },

	// { "MET", { "C", "CA", "CB", "CE", "CG", "N", "O", "SD" } },
	{ { "MET", "CA" }, { "N", "CA", "C" } },
	{ { "MET", "CB" }, { "CA", "CB", "CG" } },
	{ { "MET", "CG" }, { "CB", "CG", "SD" } },
	{ { "MET", "SD" }, { "CG", "SD", "CE" } },
	{ { "MET", "CE" }, { "SD", "CG", "CE" } },
	{ { "MET", "O" }, { "C", "CA", "O" } },

	// { "PHE", { "C", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "CZ", "N", "O" } },
	{ { "PHE", "CA" }, { "N", "CA", "C" } },
	{ { "PHE", "CB" }, { "CA", "CB", "CG" } },
	{ { "PHE", "CG" }, { "CD1", "CG", "CD2" } },
	{ { "PHE", "CD1" }, { "CG", "CD1", "CE1" } },
	{ { "PHE", "CE1" }, { "CD1", "CE1", "CZ" } },
	{ { "PHE", "CD2" }, { "CG", "CD2", "CE2" } },
	{ { "PHE", "CE2" }, { "CD1", "CE2", "CZ" } },
	{ { "PHE", "CZ" }, { "CE1", "CZ", "CE2" } },
	{ { "PHE", "O" }, { "C", "CA", "O" } },

	// { "PRO", { "C", "CA", "CB", "CD", "CG", "N", "O" } },
	{ { "PRO", "CA" }, { "N", "CA", "C" } },
	{ { "PRO", "CB" }, { "CA", "CB", "CG" } },
	{ { "PRO", "CG" }, { "CB", "CG", "CD" } },
	{ { "PRO", "CD" }, { "CG", "CD", "N" } },
	{ { "PRO", "O" }, { "C", "CA", "O" } },

	// { "SER", { "C", "CA", "CB", "N", "O", "OG" } },
	{ { "SER", "CA" }, { "N", "CA", "C" } },
	{ { "SER", "CB" }, { "CA", "CB", "OG" } },
	{ { "SER", "OG" }, { "CB", "CA", "OG" } },
	{ { "SER", "O" }, { "C", "CA", "O" } },

	// { "THR", { "C", "CA", "CB", "CG2", "N", "O", "OG1" } },
	{ { "THR", "CA" }, { "N", "CA", "C" } },
	{ { "THR", "CB" }, { "CA", "CB", "OG1" } },
	{ { "THR", "OG1" }, { "CB", "CA", "OG1" } },
	{ { "THR", "CG2" }, { "CB", "CA", "CG2" } },
	{ { "THR", "O" }, { "C", "CA", "O" } },

	// { "TRP", { "C", "CA", "CB", "CD1", "CD2", "CE2", "CE3", "CG", "CH2", "CZ2", "CZ3", "N", "NE1", "O" } },
	{ { "TRP", "CA" }, { "N", "CA", "C" } },
	{ { "TRP", "CB" }, { "CA", "CB", "CG" } },
	{ { "TRP", "CG" }, { "CD1", "CG", "CD2" } },
	{ { "TRP", "CD1" }, { "CG", "CD1", "NE1" } },
	{ { "TRP", "CD2" }, { "CE2", "CD2", "CE3" } },
	{ { "TRP", "NE1" }, { "CD1", "NE1", "CE2" } },
	{ { "TRP", "CE2" }, { "CD2", "CE2", "CZ2" } },
	{ { "TRP", "CE3" }, { "CD2", "CE3", "CZ3" } },
	{ { "TRP", "CZ2" }, { "CE2", "CZ2", "CH2" } },
	{ { "TRP", "CZ3" }, { "CE3", "CZ3", "CH2" } },
	{ { "TRP", "CH2" }, { "CZ2", "CH2", "CZ3" } },
	{ { "TRP", "O" }, { "C", "CA", "O" } },

	// { "TYR", { "C", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "CZ", "N", "O", "OH" } },
	{ { "TYR", "CA" }, { "N", "CA", "C" } },
	{ { "TYR", "CB" }, { "CA", "CB", "CG" } },
	{ { "TYR", "CG" }, { "CD1", "CG", "CD2" } },
	{ { "TYR", "CD1" }, { "CG", "CD1", "CE1" } },
	{ { "TYR", "CE1" }, { "CD1", "CE1", "CZ" } },
	{ { "TYR", "CD2" }, { "CG", "CD2", "CE2" } },
	{ { "TYR", "CE2" }, { "CD1", "CE2", "CZ" } },
	{ { "TYR", "CZ" }, { "CE1", "CZ", "CE2" } },
	{ { "TYR", "OH" }, { "CE2", "CZ", "CE1" } },
	{ { "TYR", "O" }, { "C", "CA", "O" } },

	// { "VAL", { "C", "CA", "CB", "CG1", "CG2", "N", "O" } }
	{ { "VAL", "CA" }, { "N", "CA", "C" } },
	{ { "VAL", "CB" }, { "CA", "CB", "CG1" } },
	{ { "VAL", "CG1" }, { "CB", "CA", "CG1" } },
	{ { "VAL", "CG2" }, { "CB", "CA", "CG2" } },
	{ { "VAL", "O" }, { "C", "CA", "O" } }

};

const map<string, set<pair<string, string> > > CANONICAL20_BONDS =
{
	{ "ALA", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" } } },
	{ "ARG", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD" }, { "CD", "NE" }, { "NE", "CZ" }, { "CZ", "NH1" }, { "CZ", "NH2" } } },
	{ "ASN", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "OD1" }, { "CG", "ND2" } } },
	{ "ASP", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "OD1" }, { "CG", "OD2" } } },
	{ "CYS", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "SG" } } },
	{ "GLN", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD" }, { "CD", "OE1" }, { "CD", "NE2" } } },
	{ "GLU", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD" }, { "CD", "OE1" }, { "CD", "OE2" } } },
	{ "GLY", { { "N", "CA" }, { "CA", "C" }, { "C", "O" } } },
	{ "HIS", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "ND1" }, { "ND1", "CE1" }, { "CG", "CD2" }, { "CD2", "NE2" }, { "CE1", "NE2" } } },
	{ "ILE", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG1" }, { "CG1", "CD1" }, { "CB", "CG2" } } },
	{ "LEU", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD1" }, { "CG", "CD2" } } },
	{ "LYS", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD" }, { "CD", "CE" }, { "CE", "NZ" } } },
	{ "MET", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "SD" }, { "SD", "CE" } } },
	{ "PHE", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD1" }, { "CD1", "CE1" }, { "CE1", "CZ" }, { "CG", "CD2" }, { "CD2", "CE2" }, { "CE2", "CZ" } } },
	{ "PRO", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD" }, { "CD", "N" } } },
	{ "SER", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "OG" } } },
	{ "THR", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "OG1" }, { "CB", "CG2" } } },
	{ "TRP", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD1" }, { "CG", "CD2" }, { "CD1", "NE1" }, { "CD2", "CE2" }, { "CD2", "CE3" }, { "NE1", "CE2" }, { "CE2", "CZ2" }, { "CZ2", "CH2" }, { "CH2", "CZ3" }, { "CZ3", "CE3" } } },
	{ "TYR", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG" }, { "CG", "CD1" }, { "CG", "CD2" }, { "CD1", "CE1" }, { "CD2", "CE2" }, { "CE1", "CZ" }, { "CE2", "CZ" }, { "CZ", "OH" } } },
	{ "VAL", { { "N", "CA" }, { "CA", "C" }, { "C", "O" }, { "CA", "CB" }, { "CB", "CG1" }, { "CB", "CG2" } } }
};

const map<string, char> CANONICAL20_MAP = 
{
	{ "CYS", 'C' }, { "ASP", 'D' }, { "SER", 'S' }, { "GLN", 'Q' },
	{ "LYS", 'K' }, { "ILE", 'I' }, { "PRO", 'P' }, { "THR", 'T' }, 
	{ "PHE", 'F' }, { "ASN", 'N' }, { "GLY", 'G' }, { "HIS", 'H' }, 
	{ "LEU", 'L' }, { "ARG", 'R' }, { "TRP", 'W' }, { "ALA", 'A' }, 
	{ "VAL", 'V' }, { "GLU", 'E' }, { "TYR", 'Y' }, { "MET", 'M' }, 
	{ "MSE", 'M' }, { "HYP", 'P' }, { "MLY", 'K' }, { "SEP", 'S' }, 
	{ "TPO", 'T' }, { "CSO", 'C' }, { "PTR", 'Y' }, { "KCX", 'K' }, 
	{ "CME", 'C' }, { "CSD", 'A' }, { "CAS", 'C' }, { "MLE", 'L' }, 
	{ "DAL", 'A' }, { "CGU", 'E' }, { "DLE", 'L' }, { "FME", 'M' }, 
	{ "DVA", 'V' }, { "OCS", 'C' }, { "DPR", 'P' }, { "MVA", 'V' }, 
	{ "TYS", 'Y' }, { "M3L", 'K' }, { "SMC", 'C' }, { "ALY", 'K' }, 
	{ "CSX", 'C' }, { "DCY", 'C' }, { "NLE", 'L' }, { "DGL", 'E' }, 
	{ "DSN", 'S' }, { "CSS", 'C' }, { "DLY", 'K' }, { "MLZ", 'K' }, 
	{ "DPN", 'F' }, { "DAR", 'R' }, { "PHI", 'F' }, { "IAS", 'D' }, 
	{ "DAS", 'D' }, { "HIC", 'H' }, { "MP8", 'P' }, { "DTH", 'T' }, 
	{ "DIL", 'I' }, { "MEN", 'N' }, { "DTY", 'Y' }, { "CXM", 'M' }, 
	{ "DGN", 'G' }, { "DTR", 'W' }, { "SAC", 'S' }, { "DSG", 'N' }, 
	{ "MME", 'M' }, { "MAA", 'A' }, { "YOF", 'Y' }, { "FP9", 'P' }, 
	{ "FVA", 'V' }, { "MLU", 'L' }, { "OMY", 'Y' }, { "FGA", 'E' }, 
	{ "MEA", 'F' }, { "CMH", 'C' }, { "DHI", 'H' }, { "SEC", 'C' }, 
	{ "OMZ", 'Y' }, { "SCY", 'C' }, { "MHO", 'M' }, { "MED", 'M' }, 
	{ "CAF", 'C' }, { "NIY", 'Y' }, { "OAS", 'S' }, { "SCH", 'C' }, 
	{ "MK8", 'L' }, { "SME", 'M' }, { "LYZ", 'K' }
};

}

#endif /* TOPOLOGY_H_ */
