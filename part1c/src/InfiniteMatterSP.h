/* Header file for the single particle class (SingleParticleState.h)*/

#ifndef INFINITEMATTERSP_H
#define INFINITEMATTERSP_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using std::string;
using std::vector;
using std::cout;
using std::endl;


struct qState
{
	int kx;
	int ky;
	int kz;
	int s;
};


class InfiniteMatterSP
{

	private:
		vector<qState> fStates;
		int fNmax;
		int fNmagic;
		double fRho;

	public:
		InfiniteMatterSP(int Nmax, int Nmagic, double rho);

		void GenerateSP();

		void Print();
};

#endif // INFINITEMATTERSP_H