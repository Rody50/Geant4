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
const double KAPPAR = 1.5, KAPPAS = 0.5, KAPPAT = 0.6;
const double VR = 200, VS = -90, VT = -180;

struct qState
{
	short nx, ny, nz, nshell;
	bool isUp, isHole;

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

	double rp(int ii, const vector<qState> &s)
	{
		if(ii == 0) return s[i].nx - s[j].nx;
		if(ii == 1) return s[i].ny - s[j].ny;
		if(ii == 2) return s[i].nz - s[j].nz;
	}

	void print(const vector<qState> &s) const
	{
		cout << "PAIRPAIR -- i: " << i << " \tj: " << j << endl;
		s[i].print();
		s[j].print();
	}
};

struct P_group_t
{
	array<pair_t, 10000> pr; // pairs
	int npp, nph, nhh, np; // np: n of pairs
	int P[3];

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
		return A() + H();
	}
	void print(const vector<qState> &s) const{
		cout << "(Px,Py,Pz): " << "(" << P[0] << ",";
		cout << P[1] << "," << P[2] << ")" << endl;
		cout << "npp: " << npp;
		cout << "\t nph: " << nph << "\tnhh: " << nhh << endl;
		for(int i = 0; i < np; i++){
			cout << i << ": ";
			pr[i].print(s);
		} // end for over i
	}
	void CorrelationEnergy(); // for CCD calculation
};

class InfiniteMatterSP{
	
	public:
		InfiniteMatterSP(int Nmax, int nShell, double rho);
		~InfiniteMatterSP();
		void GenerateSP();
		void ConstructPairs();
		// get the P group by (nX,nY,nZ)
		P_group_t *P(int nX, int nY, int nZ);
		double CorrelationEnergy();
		void Print() const;
		double Minnesota(const pair_t * t, const pair_t * s);

	private:
		vector<qState> fStates;
		const int fNmax;
		const int fNShell; // zero -> 2 singleSP states, itself included
		int fA, fNSP; // number of particles & single SP states
		double fRho; // nucleon density      unit: fm^-3
		P_group_t *fP; // pair grouped by total P
		const double fL;
};

#endif // INFINITEMATTERSP_H
