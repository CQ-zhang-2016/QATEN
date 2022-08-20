#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "AtomRecord.h"

namespace AtomRecord {

Atom ReadAtomRecord(char *str) {

	Atom A;

	/* complete string to 80 characters */
	char *pch = strpbrk(str, "\n");
	if (pch != NULL && pch - str < 80) {
		while (pch != str + 80) {
			*pch++ = ' ';
		}
		*pch = '\0';
	}

	/*
	 * reading PDB ATOM string in the reverse direction
	 */

	/* Charge on the atom */
	str[80] = '\0';
	strncpy(A.charge, str + 78, 3);

	/* Element symbol, right-justified */
	str[78] = '\0';
	strncpy(A.element, str + 76, 3);

	/* Segment identifier, left-justified */
	str[76] = '\0';
	strncpy(A.segId, str + 72, 5);

	/* Temperature factor */
	str[67] = '\0';
	A.temp = atof(str + 60);

	/* Occupancy */
	str[60] = '\0';
	A.occup = atof(str + 54);

	/* Z coordinate */
	str[54] = '\0';
	A.z = atof(str + 46);

	/* Y coordinate */
	str[46] = '\0';
	A.y = atof(str + 38);

	/* X coordinate */
	str[38] = '\0';
	A.x = atof(str + 30);

	/* Code for insertion of residues */
	A.insCode = str[26];

	/* Residue sequence number */
	str[26] = '\0';
	A.resNum = atoi(str + 22);

	/* Chain identifier */
	A.chainId = str[21];

	/* Residue name */
	str[20] = '\0';
	sscanf(str + 17, "%s", A.resName); /* no spaces in residue name */

	/* Alternate location indicator */
	A.altLoc = str[16];

	/* Atom name */
	str[16] = '\0';
	sscanf(str + 12, "%s", A.atomName);

	/* Atom serial number */
	str[11] = '\0';
	A.atomNum = atoi(str + 6);

	return A;

}

}
