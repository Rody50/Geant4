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


struct configurationLabel
{
	int N;
	int two_m;
};


class PairingConfiguration
{

	private:
		vector<configurationLabel> fConfiguration;
		int fN;

	public:
		PairingConfiguration(int N);

		void setfN(int N);

		void GenerateSPS();

		int getn(int i);

		int getTwo_m(int i);

		void print();
};

#endif // PairingConfiguration