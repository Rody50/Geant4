/* Header file for the FillMatrix class (FillMatrix.h) */
#ifndef FILLMATRIX_H
#define FILLMATRIX_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <bitset>
#include <cstdlib>

using std::string;
using std::vector;
using std::cout;
using std::endl;

typedef long long LL;

class FillMatrix
{
	private:
		int fNu;
		int fA;
		int fN;

		Tensor4 fV;
		Tensor4 fT;
		Tensor2 fF;

	public:
		FillMatrix(int N, int particles);

		void FillV();

		void FillF();

		void FillT();

		bool IsPair(int a, int b);


}
#endif  // FILLMATRIX_H