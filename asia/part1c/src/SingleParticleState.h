/* Header file for the single particle class (SingleParticleState.h)*/

#ifndef SINGLEPARTICLESTATE_H
#define SINGLEPARTICLESTATE_H

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
	int n;
	int l;
	int two_j;
	int two_m;
};


class SingleParticleState
{

	private:
		vector<quantumNumber> fStates;
		int fN;

	public:
		SingleParticleState(int N);

		void setfN(int N);

		void GenerateSPS();

		int getn(int i);

		int getl(int i);

		int getTwo_j(int i);

		int getTwo_m(int i);

		void print();
};

#endif // SINGLEPARTICLESTATE_H