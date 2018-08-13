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
#include <armadillo>
#include "InfiniteMatterSP.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using namespace arma;

class Tensor2;
class Tensor4;
class CCMInf;
class IMSRGInf;

class FillMatrix
{

	public:
		FillMatrix();

		static void FillV(Tensor4 & v, double g);

		static void FillV(mat & vpppp, mat & vpphh,
			mat & vhhhh, P_group_t * p, InfiniteMatterSP *InfSP_p);

		static void FillV(mat & v,
			P_group_t *p, InfiniteMatterSP *InfSP_p);
		
		static void FillF(Tensor2 & f, int A, double g, double d);

		static void FillF(mat & f,
			P_group_t *p, IMSRGInf *imsrg_p);

		static void FillT(Tensor4 & t, Tensor4 & v, Tensor2 & fp,
			Tensor2 & fh);

		static void FillT(mat & tpphh,
			P_group_t * p, CCMInf *obj);

		static void FillEta(mat & etapphh, mat & matpphh, 
			P_group_t * p, IMSRGInf * imsrg_p);

		static bool IsPair(int a, int b);
};

#endif  // FILLMATRIX_H