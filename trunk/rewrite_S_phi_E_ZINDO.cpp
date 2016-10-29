#include <iostream>
#include <cmath>
#include <iomanip>

#include "components.h"
#define ZINDO 1
#define HOMO_mon 180 //Counting from zero!

using namespace std;


int main(int argc, char * argv [])
{
    if (argc!=5) {
	cerr << "***Error*** : Should be called with (in order!):\n"
	     << "                  name of monomer_1 calculations, name of monomer_2 calculations, name of pair calculations and either zindo or dft\n"
	     << "              For each monomer and dimer should have a 'log' and a 'pun' file\n"; 
	     return (-1);
    }
    if (ZINDO) { cerr << "Assuming ZINDO calculations\n";}
    else       { cerr << "Assuming not using ZINDO\n";}
    components mon1, mon2, pair;

    mon1.readS(argv[1]);
    mon1.readMOs(argv[1]);
    //mon1.normaliseMOs();
    //mon1.checkOrthNorm();
    mon2.readS(argv[2]);
    mon2.readMOs(argv[2]);
    //mon2.normaliseMOs();
    pair.readS(argv[3]);
    pair.readMOs(argv[3]);
    if (mon1.Nbasis!=mon2.Nbasis || pair.Nbasis!=2*mon1.Nbasis) {
	cerr << "***Error***: Unexpected basis sets\n";
	return (-2);
    }
    if(!ZINDO) {
	//pair.normaliseMOs();
	mon1.justPrintS("S_1.txt");
	mon1.justPrintMOs("MOs_1.txt");
	mon2.justPrintS("S_2.txt");
	mon2.justPrintMOs("MOs_2.txt");
	pair.justPrintS("S_pair.txt");
	pair.justPrintMOs("MOs_pair.txt");
	pair.justPrintEvls("Evls_pair.txt");
    }
    if (ZINDO) {
	/********************************************************************************
	* If using zindo orbitals are already orthogonal and can simply calculate F_local
	********************************************************************************/
	double F[pair.Nbasis][pair.Nbasis];
	double B[pair.Nbasis][pair.Nbasis]; 	//B[i][k] = MOs_mons[i][j] * MOs_pair_[j][k]
						//     = mon1[i][j] * MOs_pair[j][k] 	 	for 0   < i,j < N/2
						//     = mon2[i-N/2][j-N/2] * MOs_pair[j][k]    for N/2 < i,j < N
	for (int k=0; k<pair.Nbasis; k++) {
	    for (int i=0; i<mon1.Nbasis; i++) {
		B[i][k]=0.0;
		for (int j=0; j<mon1.Nbasis; j++) {
		    B[i][k]+=mon1.MOs[j][i] * pair.MOs[j][k];
		}
	    }
	    for (int i=mon2.Nbasis; i<pair.Nbasis; i++) {
		B[i][k]=0.0;
		for (int j=mon2.Nbasis; j<pair.Nbasis; j++) {
		    B[i][k]+=mon2.MOs[j-mon2.Nbasis][i-mon2.Nbasis] * pair.MOs[j][k];	
		}
	    }
	}
	for (int i=0; i<pair.Nbasis; i++) {
	    for (int j=0; j<pair.Nbasis; j++) {
		F[i][j]=0.0;
		for (int k=0; k<pair.Nbasis; k++) {
		    F[i][j]+=B[i][k] * pair.Evls[k] * B[j][k];
		}
	    }
	}
	cout << mon1.Nbasis << endl;
	cout << F[HOMO_mon][HOMO_mon+mon1.Nbasis] << endl << F[HOMO_mon+mon1.Nbasis][HOMO_mon] << endl;
    }//if(ZINDO)
	return 0;
	
}
