/* Header file for the single particle class (SingleParticleState.h)*/

#ifndef INFINITEMATTERSP_H
#define INFINITEMATTERSP_H

#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using std::string;
using std::vector;
using std::array;
using std::cout;
using std::endl;

const double PI = 3.14159265357;
const double KAPPAR = 1.487, KAPPAS = 0.465, KAPPAT = 0.639;
const double VR = 200., VS = -91.85, VT = -178.;
const double hbarc = 197.326; // MeV * fm
const double nmass = 938.92; //939.565; // MeV / c^2

struct qState
{
	short nx, ny, nz, nshell;
	bool isUp, isHole;

	short spin() const
	{
		return 2 * isUp - 1;
	}

	double r(){
		return sqrt(nx*nx + ny*ny + nz*nz);
	}

	void print() const{
		cout << "(nx,ny,nz): " << "(" << nx << ",";
		cout << ny << "," << nz <<  ")";
		cout << " \tisUp: " << isUp;
		cout << " \tnshell: " << nshell << " \tisHole: " << isHole << endl;
	}
};

struct pair_t
{
	int i, j;

	int P(int ii, const vector<qState> &s) const
	{
		if(ii == 0) return s[i].nx + s[j].nx;
		if(ii == 1) return s[i].ny + s[j].ny;
		if(ii == 2) return s[i].nz + s[j].nz;
		return 0.;
	}

	double rp(int ii, const vector<qState> &s) const
	{
		if(ii == 0) return s[i].nx - s[j].nx;
		if(ii == 1) return s[i].ny - s[j].ny;
		if(ii == 2) return s[i].nz - s[j].nz;
		return 0.;
	}

	void print(const vector<qState> &s) const
	{
		cout << "PAIRPAIR -- i: " << i << " \tj: " << j << endl;
		s[i].print();
		s[j].print();
	}
};

struct P_group_t;

class InfiniteMatterSP{
	
	public:
		InfiniteMatterSP(int Nmax, int nShell, double rho);
		~InfiniteMatterSP();
		void GenerateSP();
		double GetfEnUn(){return fEnUn;}
		vector<qState> & GetfStates(){return fStates;}
		void ConstructPairs();
		// get the P group by (nX,nY,nZ)
		P_group_t *P(int nX, int nY, int nZ);
		double CorrelationEnergy(double deltaS);
		void Print() const;
		double Minnesota(const pair_t * t, const pair_t * s);
		double Minnesota(const qState &ti, const qState &tj,
		const qState &si, const qState &sj);

	private:
		vector<qState> fStates;
		const int fNmax;
		const int fNShell; // zero -> 2 singleSP states, itself included
		int fA, fNSP; // number of particles & single SP states
		double fRho; // nucleon density      unit: fm^-3
		P_group_t *fP; // pair grouped by total P
		double fL;
		double fEnUn;
};

struct corr_en_t{
	int P[3];
	double corr_en;

	corr_en_t(){
		P[0] = -1; P[1] = -1; P[2] = -1;
		corr_en = -99.;
	}
	corr_en_t(int *p, double cor){
		P[0] = p[0]; P[1] = p[1]; P[2] = p[2];
		corr_en = cor;
	}
};


struct P_group_t
{
	array<pair_t, 10000> pr; // pairs
	int npp, nph, nhh, np; // np: n of pairs
	int P[3];
	static vector<corr_en_t> co_enVec;

	P_group_t()
	{
		npp = 0; nph = 0; nhh = 0; np = 0;
		P[0] = -99; P[1] = -99; P[2] = -99;
	}
	double A()
	{
		return 0.5 + sqrt(0.25 + 2*npp);
	}
	double Nu()
	{
		return 0.5+sqrt(0.25+2*nhh);
	}
	double N()
	{
		return A() + Nu();
	}
	void print(const vector<qState> &s) const{
		cout << "(Px,Py,Pz): " << "(" << P[0] << ",";
		cout << P[1] << "," << P[2] << ")" << endl;
		cout << endl;
//		return;
		cout << "npp: " << npp;
		cout << "\t nph: " << nph << "\tnhh: " << nhh << endl;
		for(int i = 0; i < np; i++){
			cout << i << ": ";
			pr[i].print(s);
		} // end for over i
	}
	double Psqrt(){
		return sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
	}
	double CorrelationEnergy(InfiniteMatterSP *infSP); // for CCD calculation
	double CorrelationEnergyIMSRG(InfiniteMatterSP *infSP, double deltaS); // for CCD calculation
};



#endif // INFINITEMATTERSP_H