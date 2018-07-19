#ifndef FCI_H
#define FCI_H


#include "PairingConfiguration.h"
#include "armadillo"

using namespace arma;

class FCI
{
	private:
		int fN;
		vector<LL> fConfiguration;
		mat fHamiltonian;

	public:
		FCI(const vector<LL> & configuration, int N);

		double MatrixElement(int config, int config_p, double g, double d);

		void GenerateMatrix(double g, double d);

		double eigenValue();

		void Print();

		void WriteToFile();

};

#endif // FCI_H