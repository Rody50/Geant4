/* Header file for the single particle class (PairingConfiguration.h)*/

#ifndef PairingConfiguration
#define PairingConfiguration

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using std::string;
using std::vector;
using std::cout;
using std::endl;

class PairingConfiguration
{

	private:
		int fPairs;
		int fN;
		vector<bool*> fConfigurations;

	public:
		PairingConfiguration(int N);

		void setfN(int N);

		void GenerateSPS();

		int getn(int i);

		int getTwo_m(int i);

		void print();
};

#endif // PairingConfiguration