/* Header file for the CCM class (CCM.h) */
#ifndef CCM_H
#define CCM_H

#include "Tensor2.h"
#include "Tensor4.h"
#include "FillMatrix.h"

class CCM
{
	public:
		Tensor4 fHamil;
		Tensor4 fVpppp, fVpphh, fVhhhh;
		Tensor4 fTpppp, fTpphh, fThhhh;
		Tensor2 fFpp, fFhh;
		int fNu;
		int fA;
		int fN;

	private:
		CCM(int A, int N);

		void ComputeH();

		void SolveT();

}

#endif  // CCM_H