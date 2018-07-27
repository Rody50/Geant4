/* Header file for the single particle class (SingleParticleState.h)*/

#ifndef INFINITEMATTERSP_H
#define INFINITEMATTERSP_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using std::string;
using std::vector;
using std::cout;
using std::endl;


struct quantumNumber
{
	int kx;
	int ky;
	int kz;
	int s;
};


class InfiniteMatterSP
{

	private:
		vector<quantumNumber> fStates;
		int fN;

	public:
		InfiniteMatterSP(int N);

		// void setfN(int N);

		// void GenerateSPS();

		// int getn(int i);

		// int getl(int i);

		// int getTwo_j(int i);

		// int getTwo_m(int i);

		// void print();
};

#endif // INFINITEMATTERSP_H