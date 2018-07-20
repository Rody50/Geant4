#include "SingleParticleState.h"
#include "PairingConfiguration.h"
#include "FCI.h"

int main()
{
    SingleParticleState sPS(2);

	sPS.GenerateSPS();
	sPS.print();

	int N = 4;
	int no_pairs = 2;

	FCI hamiltonian(N, no_pairs);

	hamiltonian.GenerateMatrix(1.0, 1.0);
	hamiltonian.Print();

}
