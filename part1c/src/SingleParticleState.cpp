/* Implementation of the SingleParticleState class (SingleParticleState.cpp)*/

#include "SingleParticleState.h"

SingleParticleState::SingleParticleState(int N) : fN(N)
{}

void SingleParticleState::setfN(int N)
{
	fN = N;
}

void SingleParticleState::GenerateSPS()
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
				fStates.push_back(temp);
			}
		}	
	}
}

int SingleParticleState::getn(int i)
{
	return fStates[i].n;
}

int SingleParticleState::getl(int i)
{
	return fStates[i].l;
}

int SingleParticleState::getTwo_j(int i)
{
	return fStates[i].two_j;
}

int SingleParticleState::getTwo_m(int i)
{
	return fStates[i].two_m;
}

void SingleParticleState::print()
{
	cout << "The quantum numbers are: " << endl;
	cout << "n \t" << "l \t" << "twoj \t" << "two_m" << endl;
	for (auto t : fStates)
		cout << t.n << " \t" << t.l << " \t" << t.two_j << " \t" << t.two_m << endl;
	cout << "Done!" << endl;
}