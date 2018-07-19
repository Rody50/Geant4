/* Implementation of the PairingConfiguration class (PairingConfiguration.cpp)*/

#include "PairingConfiguration.h"

PairingConfiguration::PairingConfiguration(int N) : fN(N)
{}

void PairingConfiguration::setfN(int N)
{
	fN = N;
}

void PairingConfiguration::GenerateSPS()
{
	quantumNumber temp;

	for (int n = 0; fN - 2 * n >= 0; n++)
	{
		temp.n = n;

		int l = fN - 2 * n; 
		temp.l = l;

		int temp_two_j = 2 * l - 1;

		for (int twoj = temp_two_j; twoj <= temp_two_j + 2; twoj+=2)
		{
			temp.two_j = twoj;
			for (int two_mj = -twoj; two_mj <= twoj; two_mj+=2)
			{
				temp.two_m = two_mj;
				fConfigurations.push_back(temp);
			}
		}	
	}
}

int PairingConfiguration::getn(int i)
{
	return fConfigurations[i].n;
}

int PairingConfiguration::getl(int i)
{
	return fConfigurations[i].l;
}

int PairingConfiguration::getTwo_j(int i)
{
	return fConfigurations[i].two_j;
}

int PairingConfiguration::getTwo_m(int i)
{
	return fConfigurations[i].two_m;
}

void PairingConfiguration::print()
{
	cout << "The quantum numbers are: " << endl;
	cout << "n \t" << "l \t" << "twoj \t" << "two_m" << endl;
	for (auto t : fConfigurations)
		cout << t.n << " \t" << t.l << " \t" << t.two_j << " \t" << t.two_m << endl;
	cout << "Done!" << endl;
}