#include "SingleParticleState.h"
#include "PairingConfiguration.h"
#include "FCI.h"

int main()
{
    SingleParticleState sPS(2);

	sPS.GenerateSPS();
	sPS.print();

	int N = 10;
	int no_pairs = 2;

	FCI hamiltonian(N, no_pairs);

	hamiltonian.GenerateMatrix(1.0, 1.0);
	hamiltonian.Print();
	
	double g = -1.;
	double eigenVec00, eigenVal0;
	std::ofstream file0("results/eigVal0_10.txt");
	std::ofstream file1("results/eigVec00_10.txt");
	
	double Nvec[4] = {4, 6, 8, 10};
	for(int i = 0 ; i < 21; i++){
		g += 0.1;
		hamiltonian.GenerateMatrix(g, 1.0);
		hamiltonian.eigenValue(eigenVal0, eigenVec00);
		file0 << g << "\t" << eigenVal0 - (2*1.0-g) << endl;
		file1 << g << "\t" << eigenVec00 << endl;
	}

}
