#ifndef FILE_COMP
#define FILE_COMP

#include <iostream>
#include <fstream>

using namespace std;

class components
{
    public: 
	int Nat;
	int Nbasis;
	int NBsUse;	//The number of MOs produced (may differ from AOs if basis has been collapsed)
	double ** MOs;
	double ** S;
	double * Evls;
	string * Atypes;

	components(){
	    Nat=0;
	    Nbasis=0;
	}

	~components(){
	    if (Nbasis != 0 ){
		for (int i=0; i <Nbasis;i++) {
		    delete [] MOs[i];
		    delete [] S[i];
		}
		delete [] S;
		delete [] MOs;
		delete [] Evls;
		//delete [] Atypes;
	    }
	    MOs=0;
	    S=0;
	    Atypes=0;
	    Evls=0;
	    Nbasis=0;
	}

	void readS(char * namefile){
	    //read all the info from a punch file
		char infile[256];
		strcpy(infile,namefile);
		strcat(infile,".log");
		ifstream in(infile);	 
		if (!in)       {cerr << "***Error*** : " << infile << " doesn't exist.\n";}
	        string word;	
		while( in ){
		    getline(in, word);
		    if (   word == " Two-electron integral symmetry is turned off." 
		        || word == " Two-electron integral symmetry is turned on.") {
			in >> word;
			sscanf(word.c_str(), "%i", &Nbasis);
			cerr << "NBasis = " << Nbasis << endl;
		    }
		    if (word.find("*** Overlap")!=string::npos) {
			S = new double * [Nbasis];
			for (int m=0;m<Nbasis;m++) {S[m] = new double[Nbasis];}
			for (int j=0; j<Nbasis; j+=5) {
    			    getline(in,word);				//Get labels 
			    int line=0;
			    for (int i=j; i<Nbasis; i++) {
				in >> word;        			//Get label
				if (line<4) {      			//Don't have all fields
				    for (int l=0; l<=line; l++) {
					in >> word;			
					size_t pos=word.find("D");
					word=word.substr(0,pos)+"e"+word.substr(pos+1);
					S[i][j+l]=S[j+l][i]=atof(word.c_str());
				    }
				}
				else{					//Have 5 fields
				    for (int l=0; l<=4; l++) {
					in >> word;
					size_t pos=word.find("D");
					word=word.substr(0,pos)+"e"+word.substr(pos+1);
					S[i][j+l]=S[j+l][i]=atof(word.c_str());
				    }
				}
				line++;
			    }//for(i)
			    getline(in,word);				//Appears to be necessary to finish line
			}//for(j)
		   }//if(*** Overlap)
		}//while(in)	
		cerr << "Read in S...  ";
	}//components(char * namefile)
	void readMOs(char * namefile) {
	    char infile[256];
		    //TIDY
	    char logfile[256];
	    strcpy(logfile,namefile);
	    strcat(logfile,".log");
	    ifstream logIn(logfile);
	    string word;
	    while (logfile) {
		getline(logIn,word);						
	    	if (word.find("NBsUse=")!=string::npos) {
		    NBsUse=atoi(word.substr(8,6).c_str());
		    cerr << "NBsUse = " << NBsUse << endl;
		    break;
		}
	    }
	    if (Nbasis==0 || NBsUse==0) {cerr << "***Error*** : Trying to read MOs without knowing the number of basis sets\n";}
	    strcpy(infile,namefile);
	    strcat(infile,".pun");
	    ifstream in(infile);
	    if (!in)       {cerr << "***Error*** : " << infile << " doesn't exist.\n";}
	    getline(in,word);						//First line of rubbish
	    MOs = new double * [Nbasis];
	    for (int i=0; i<Nbasis; i++) {MOs[i] = new double [NBsUse];}
	    Evls = new double [NBsUse];
	    //N.B. MOs[i][j] = <Basis_i | MO_j>
	    for (int i=0; i<NBsUse; i++) {
		getline(in,word);					//MO label and energy
	    	size_t pos=word.find("OE=");
    		word=word.substr(pos+3);
		pos=word.find("D");
		word=word.substr(0,pos)+"e"+word.substr(pos+1);
		Evls[i]=27.211396*atof(word.c_str());			//eV
		for (int j=0; j<Nbasis; j+=5) {
		    getline(in,word);
		    if (Nbasis-j>=5) {
			for (int k=0; k<5; k++) {
			    string f=word.substr(15*k, 15*(k+1));
			    f=f.substr(0,(f.find("D")))+"e"+f.substr(f.find("D")+1);
			    MOs[j+k][i]=atof(f.c_str());
			}//for(k)
		    }//if(Nbasis-j>5)
		    else {
			for (int k=0; k<Nbasis-j; k++) {
			    string f=word.substr(15*k, 15*(k+1));
			    f=f.substr(0,(f.find("D")))+"e"+f.substr(f.find("D")+1);
			    MOs[j+k][i]=atof(f.c_str());
			}//for(k)
		    }//else
		}//for(Nbasis)
	    }//for(Nbasis) 
	    cerr << "Read in MOs\n";
	}//readMOs(char * namefile)
	void checkOrthNorm() {
	    for (int psi1=0; psi1<NBsUse; psi1++) {
		for (int psi2=psi1; psi2<NBsUse; psi2++) {
		    double dotProd=0.0;
		    for (int basis=0; basis<Nbasis; basis++) {
			dotProd+=MOs[basis][psi1]*MOs[basis][psi2];
		    }
		    if (dotProd>1e-2) cerr << "MO_" << psi1 << ".MO_" << psi2 << " = " << dotProd << endl;
		}
	    }
	}
	void normaliseMOs() {
	    for (int psi=0; psi<Nbasis; psi++) {
		double N_sqd=0.0;
		for (int basis=0; basis<Nbasis; basis++) {
		    N_sqd+=MOs[basis][psi]*MOs[basis][psi];
		}
		double N=sqrt(N_sqd); cout << N << '\t' ;
		for (int basis=0; basis<Nbasis; basis++) {
		    MOs[basis][psi]/=N;
		}
	    }
	}//normaliseMOs()

	void set (char * namefile){
	    //read all the info from a punch file
		ifstream in(namefile);	 
		if (!in) {cerr << "Error file non existant" <<endl;}
	        string word;	
		while( in ){
		    in >> word;
		    if ( word == "(Angstroms)"){
			getline(in, word);
			getline(in, word);
			in >> word;
			Atypes = new string[Nat];
			in >> word;
			in >> word;
			in >> word;
			for (int i=0; i< Nat;i++)
			{
			    in >> Atypes[i];
			    in >> word;
			    in >> word;
			    in >> word;
			}
		    }
		    if ( word == "Overlap" ){
				in >> word;
				in >> word;
				in >> word;
				in >> word;
				in >> word;
				in >> word;
				sscanf(word.c_str(), "%i", &Nbasis);	   
				S = new double * [Nbasis];
				for (int m=0;m<Nbasis;m++)
				{
				    S[m] = new double[Nbasis];
				}
				int i=0;
				int j=0;
				while ( i< Nbasis){
				    in >>word;
				    sscanf(word.c_str(), "%lf", &S[i][j]);
				    j++;
				    if(j==Nbasis) { 
					j=0;
					i++;
			   		}
				}
		    }

		    if ( word == "vectors" ){
			in >> word;
			in >> word;
			in >> word;
			in >> word;
			in >> word;
			in >> word;
			in >> word;
			in >> word;
			in >> word;
			sscanf(word.c_str(), "%i", &Nbasis);	   
			MOs = new double * [Nbasis];
			for (int m=0;m<Nbasis;m++)
			{
			    MOs[m] = new double[Nbasis];
			}
			getline(in, word); //skip the rest of the line
			int i=0;
			int j=0;
			while ( i< Nbasis){
			    in >>word;
			    sscanf(word.c_str(), "%lf", &MOs[i][j]);
			    j++;
			    if(j==Nbasis) { 
				j=0;
				i++;
			    }
			}
		    }
		}
	}


	void justPrintS(char * namefile) {
	    ofstream out(namefile);
	    if (!out) {
		cerr << "***ERROR*** : Can't open " << namefile << endl;
		exit(-1);
	    }
	    for (int i=0; i<Nbasis; i++) {
		for (int j=0; j<Nbasis; j++) {
		    out << scientific << setprecision(8) << S[i][j] << "  "; 
		}
		out << endl;
	    }
	}
	void justPrintMOs(char * namefile) {
	    ofstream out(namefile);
	    if (!out) {
		cerr << "***ERROR*** : Can't open " << namefile << endl;
		exit(-1);
	    }
	    for (int i=0; i<Nbasis; i++) {
		for (int j=0; j<Nbasis; j++) {
		    out << scientific << setprecision(8) << MOs[i][j] << "  "; 
		}
		out << endl;
	    }
	}
	void justPrintEvls(char * namefile) {
	    ofstream out(namefile);
	    if (!out) {
		cerr << "***ERROR*** : Can't open " << namefile << endl;
		exit(-1);
	    }
	    for (int i=0; i<Nbasis; i++ ) 
		out << scientific << setprecision(8) << Evls[i] << endl;
	}

	void print_S(){
	    int i,j,k,l;
	    double test_S;
	    for (i=0; i<Nbasis;i++){
		for( j=0;j<Nbasis;j++){
		    test_S=0.0;
		    for (k=0;k<Nbasis;k++){
			for (l=0;l<Nbasis;l++){
				test_S+= S[k][l] * MOs[i][l] * MOs[j][k] ;
			}
		    }
		    cout << "test_S" << '\t' <<  test_S <<'\t' ;
		}
		cout <<endl;
	    }
	}
	
	
};

#endif // FILE_COMP
