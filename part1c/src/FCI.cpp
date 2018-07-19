
#include "FCI.h"

FCI::FCI(const vector<LL> & configuration, int N) : 
fConfiguration(configuration), fN(N), 
fHamiltonian(configuration.size(), configuration.size())
{}

double FCI::MatrixElement(int config, int config_prime, double g, double d)
{
	int energy = 0;
	int n = 0;

	for (int state = 0; state < fN; state++)
	{
		short i = (config>>state) & 0x1;
		short i_prime = (config_prime>>state) & 0x1; 

		if (i == i_prime && i == 1)
			energy += 2 * state;

		if (i != i_prime) n++;
	}

	if (n ==  0) return energy * d - g;
	else if (n == 2) return - g / 2.0;
	else return 0.0;
}

void FCI::GenerateMatrix(double g, double d)
{
	int dim = fConfiguration.size();

	for (int i = 0; i < dim; i++)
		for (int j = i; j < dim; j++)
		{
			double ME = 
				MatrixElement(fConfiguration[i], fConfiguration[j], g, d);
			fHamiltonian(i + j * dim) = ME;
			
			if (i != j)
				fHamiltonian(j + i * dim) = ME;
		}
}

double FCI::eigenValue()
{
	int length = fConfiguration.size();

	vec eigval(length);

	eig_sym(eigval, fHamiltonian);

	return eigval[0];
}

void FCI::Print()
{
	fHamiltonian.print();
}

void FCI::WriteToFile()
{

}