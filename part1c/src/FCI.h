#ifndef FCI_H
#define FCI_H

#include <fstream>
#include "PairingConfiguration.h"
#include "armadillo"

using namespace arma;

class FCI
{
	private:
		int fN, fN_pairs;
		PairingConfiguration fPairConfig;
		vector<LL> fConfiguration;
		mat *fHamiltonian;		

	public:
		FCI(int N, int n_pairs);

		double MatrixElement(int config, int config_p, double g, double d);

		void GenerateMatrix(double g, double d);

		void eigenValue(double &eigenValue, double &eigenVec00);

		void Print();

		void WriteToFile();

};

#endif // FCI_H
