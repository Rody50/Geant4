/* Header file for the single particle class (PairingConfiguration.h)*/

#ifndef PAIRINGCONFIGURATION_H
#define PAIRINGCONFIGURATION_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <bitset>
#include <cstdlib>

using std::string;
using std::vector;
using std::cout;
using std::endl;

typedef long long LL;

class PairingConfiguration
{

	private:
		const int fPairs;
		const int fN;
		vector<LL> fConfigurations;

	public:
		PairingConfiguration(int N, int pairs);

		// void setfN(int N);

		// void setfPairs(int pairs);

		int getfN();

		int getfPairs();

		const vector<LL> & getConfigVec();

		int CountBits1(LL config, int option);

		void GenerateConfig();

		void Print();
};

#endif // PAIRINGCONFIGURATION_H