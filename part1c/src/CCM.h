/* Header file for the CCM class (CCM.h) */
#ifndef CCM_H
#define CCM_H

#include "Tensor2.h"
#include "Tensor4.h"
#include "FillMatrix.h"

class CCM
{
	private:
		Tensor4 fHamil;
		Tensor4 fVpppp, fVpphh, fVhhhh, fVhhpp;
		Tensor4 fTpphh;
		Tensor2 fFpp, fFhh;
		int fNu, fA, fN;
		double fG, fD;

	public:
		CCM(int A, int N, double g, double d);

		void ComputeH();

		void SolveT();

		void Print();

};

#endif  // CCM_H
