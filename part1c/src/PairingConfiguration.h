/* Header file for the single particle class (PairingConfiguration.h)*/

#ifndef PairingConfiguration
#define PairingConfiguration

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>

using std::string;
using std::vector;
using std::cout;
using std::endl;

typedef long long LL;

class PairingConfiguration
{

	private:
		int fPairs;
		int fN;
		vector<LL> fConfigurations;

	public:
		PairingConfiguration(int N, int pairs);

		void setfN(int N);

		void setfPairs(int pairs);

		int getfN();

		int getfPairs();

		int CountBits1(LL config, int option);

		void GenerateConfig();

		void Print();
};

#endif // PairingConfiguration