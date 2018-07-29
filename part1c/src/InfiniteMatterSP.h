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

const double PI = 3.14159265357;

struct qState
{
	double nx;
	double ny;
	double nz;
	int two_s;
};

struct conState
{
	int part1;
	int part2;
};


class InfiniteMatterSP
{

	private:
		vector<qState> fStates;
		vector<conState> fConStates;
		int fNmax;
		int fNShell;
		double fRho;

	public:
		InfiniteMatterSP(int Nmax, int nShell, double rho);

		void GenerateSP();

		void Connected();

		void Print();
};

#endif // INFINITEMATTERSP_H