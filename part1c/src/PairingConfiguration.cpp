/* Implementation of the PairingConfiguration class (PairingConfiguration.cpp)*/

#include "PairingConfiguration.h"

PairingConfiguration::PairingConfiguration(int N, int pairs) : fN(N),
 fPairs(pairs)
{}

void PairingConfiguration::setfN(int N)
{
	fN = N;
}

void PairingConfiguration::setfPairs(int pairs)
{
	fPairs = pairs;
}

int PairingConfiguration::getfN()
{
	return fN;
}

int PairingConfiguration::getfPairs()
{
	return fPairs;
}

int PairingConfiguration::CountBits1(LL config, int option) //Calculate the "1" bits of an integer.
{
	int n = 0;

	if (option == 1)
	{
		do
		{
			config &= config-1;
			n++;
		} while(config!=0);

		return n;
	}

	for (int state = 0; state < fN; state++)
	{
		LL ni = pow(2, state);
		short ii = config & ni;
		if(ii == 1) n++;
	}

	return n;
}


void PairingConfiguration::GenerateConfig()
{
	LL N_max = pow(2, fN) - pow(2, fN - fPairs);

	for (int config = N_max; config = 0; config--)
	{
		if (CountBits1(config, 0) == fPairs)
			fConfigurations.push_back(config);
	}

}

void PairingConfiguration::Print()
{
	cout << "The different configurations are: " << endl;
	
	for (auto t : fConfigurations)
	{
		std::bitset<fN> y(t);
		cout << y << endl;
	}

}