/* Implementation of the SingleParticleState class (SingleParticleState.cpp)*/

#include "InfiniteMatterSP.h"

InfiniteMatterSP::InfiniteMatterSP(int Nmax, int Nmagic, double rho) : 
	fNmax(Nmax), fNmagic(Nmagic), fRho(rho)
{}

void InfiniteMatterSP::GenerateSP()
{
	qState temp;
	double L = pow(fNmagic / fRho, 1 / 3.);

	// for (int nx = -fNmax; nx <= fNmax; nx++)
	// 	for (int ny = -fNmax; ny <= fNmax; ny++)
	// 		for (int nz = -fNmax; nz <= fNmax; nz++)

	// 	{
	// 		temp.two_j = twoj;
	// 		for (int two_mj = -twoj; two_mj <= twoj; two_mj+=2)
	// 		{
	// 			temp.two_m = two_mj;
	// 			fStates.push_back(temp);
	// 		}
	// 	}	

}

void InfiniteMatterSP::Print()
{
	// cout << "The quantum numbers are: " << endl;
	// cout << "n \t" << "l \t" << "twoj \t" << "two_m" << endl;
	// for (auto t : fStates)
	// 	cout << t.n << " \t" << t.l << " \t" << t.two_j << " \t" << t.two_m << endl;
	// cout << "Done!" << endl;
}