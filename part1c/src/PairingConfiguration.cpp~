/* Implementation of the PairingConfiguration class (PairingConfiguration.cpp)*/

#include "PairingConfiguration.h"

LL InverseBit(LL num, int N); // inverse the bit ordering of an integer, N: number of total bits

PairingConfiguration::PairingConfiguration(int N, int pairs) : fN(N),
 fPairs(pairs)
{}

// void PairingConfiguration::setfN(int N)
// {
// 	fN = N;
// }

// void PairingConfiguration::setfPairs(int pairs)
// {
// 	fPairs = pairs;
// }

int PairingConfiguration::getfN()
{
	return fN;
}

int PairingConfiguration::getfPairs()
{
	return fPairs;
}

const vector<LL> & PairingConfiguration::getConfigVec()
{
	return fConfigurations;
}

int PairingConfiguration::CountBits1(LL config, int option) //Calculate the "1" bits of an integer.
{
	if (config == 0) return 0;

	int n = 0;

	if (option == 1)
	{
		do
		{
			config &= config - 1;
			n++;
		} while(config!=0);

		return n;
	}

	for (int state = 0; state < fN; state++)
	{
		short ii = (config>>state) & 0x1;
		if(ii == 1) n++;
	}

	return n;
}


void PairingConfiguration::GenerateConfig()
{
	LL N_max = pow(2, fN) - pow(2, fN - fPairs);

	for (LL config = N_max; config >= 0; config--)
	{
		if (CountBits1(config, 0) == fPairs)
			fConfigurations.push_back(InverseBit(config, fN)); // 
	}

}

void PairingConfiguration::Print()
{
	cout << "The different configurations are: " << endl;
	cout << "The size of the vector is: " << fConfigurations.size() << endl;
	
	for (auto t : fConfigurations)
	{
		for (int state = 0; state < fN; state++)
		{
			short ii = (t>>state) & 0x1;
			cout << ii;
		}	
		cout << endl;
	}

}

LL InverseBit(LL num, int N){ // inverse the bit ordering of an integer, N: number of total bits
	if(0 == num) return 0;

	LL num_in = 0;
	short *bitVec = new short[N], size = 0;
	for(int i = N;i--;) bitVec[i] = 0;
	while(num > 0){
		bitVec[size++] = num & 0x1;
		num >>= 1;
	}
	for(int i = N - 1; i >= 0; i--){
		num_in += bitVec[i] << N - i - 1;
	}
		
	delete bitVec;
	return num_in;
}

