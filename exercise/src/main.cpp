#include "SingleParticleState.h"
#include "PairingConfiguration.h"
#include "InfiniteMatterSP.h"
#include "FCI.h"
#include "CCM.h"
#include "CCMInf.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using std::ofstream;
using std::to_string;
using std::flush;
using std::vector;
using std::cin;

int main()
{

// <- <- <- // HO // -> -> -> //
// Setting up the harmonic oscillator single particle basis
  SingleParticleState sPS(2);
	sPS.GenerateSPS();
//	sPS.print();
// <- <- <- // end HO // -> -> -> //


// <- <- <- // FCI // -> -> -> //
	// int N = 20;
	// int no_pairs = 2;

	// FCI hamiltonian(N, no_pairs);

	// hamiltonian.GenerateMatrix(1.0, 1.0);
	// hamiltonian.Print();
	
	// double g = -1.;
	// double eigenVec00, eigenVal0;
	// std::ofstream file0("results/FCI/eigVal0_20.txt");
	// std::ofstream file1("results/FCI/eigVec00_20.txt");
	
	// for(int i = 0 ; i < 21; i++)
	// {
	// 	g += 0.1;
	// 	hamiltonian.GenerateMatrix(g, 1.0);
	// 	hamiltonian.eigenValue(eigenVal0, eigenVec00);
	// 	file0 << g << "\t" << eigenVal0 - (2*1.0-g) << endl;
	// 	file1 << g << "\t" << eigenVec00 << endl;
	// }
// <- <- <- // end FCI // -> -> -> //


// <- <- <- // CCM // -> -> -> //
	int nLevels = 8;
	int particles = 4;

	CCM coupled(particles, nLevels, .5, 1.);
	cout << "The Hamiltonian: " << endl;
	coupled.SolveT();

	double g = -1.1;
	double corr_en;
	std::ofstream file0("results/CCM/CCM_N_" + to_string(nLevels)
		+ "_A_" + to_string(particles) + "_0.txt");
	
	for(int i = 0 ; i < 21; i++)
	{
		g += 0.1;
		CCM coupled(particles, nLevels, g, 1.);
		corr_en = coupled.SolveT();
		file0 << g << "\t" << corr_en << endl;
	}
// <- <- <- // end CCM // -> -> -> //

// <- <- <- // CCM INF // -> -> -> //
// 	int nshell = 2, NMax = 2;
// 	double rho = 0.001; // nucleon density      unit: fm^-3

// 	InfiniteMatterSP spMatter(NMax, nshell, rho);
// 	spMatter.GenerateSP();
// 	spMatter.ConstructPairs();
// 	double corr_en = spMatter.CorrelationEnergy();
// 	cout << "The total correlation energy is: " << corr_en << endl;
// 	return 0;
	
// 	double step = 0.15 / 10.;
// 	cout << "nshell: "; cin >> nshell; cout << "\ninput nshell: " << nshell << endl;
// //	cout << "NMax: "; cin >> NMax; cout << "\ninput NMax: " << NMax << endl;
// 	ofstream file(("results/corr_en_inf_nshell" + to_string(nshell) 
// 		+ "_Nmax_" + to_string(NMax) + ".txt").c_str());
// 	for(int i = 0; i < 10; i++){
// 		cout << "---------------- rho: " << rho << " ----------------------" << endl;
// 		InfiniteMatterSP *spMatter = new InfiniteMatterSP(NMax, nshell, rho);
// 		spMatter->GenerateSP();
// 	//	spMatter.Print();
// 		spMatter->ConstructPairs();
// 		double corr_en = spMatter->CorrelationEnergy();
// 		cout << "The total correlation energy is: " << corr_en << endl;
// 		file << rho << " " << corr_en << endl;
// 		delete spMatter;
// 		rho += step;
// 	} // end loop over 1
// <- <- <- // end CCM INF // -> -> -> //


// <- <- <- // IMSRG INF // -> -> -> //
	// int nshell = 1, NMax = 1;
	// double rho = 0.08; // nucleon density unit: fm^-3
	// double deltaS = 0.1;

	// InfiniteMatterSP spMatter(NMax, nshell, rho);
	// spMatter.GenerateSP();
	// spMatter.ConstructPairs();
	// double result = spMatter.CorrelationEnergy(deltaS);

// <- <- <- // end IMSRG INF // -> -> -> //

	return 0;
}

