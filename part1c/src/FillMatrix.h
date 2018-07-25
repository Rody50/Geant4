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

class FillMatrix
{

	public:
		FillMatrix();

		static void FillV(Tensor4 & v);

		static void FillF(Tensor2 & f, bool isHole);

		static void FillT(Tensor4 & t, Tensor4 & v, Tensor2 & fp,
			Tensor2 & fh);

		static bool IsPair(int a, int b);
}
#endif  // FILLMATRIX_H