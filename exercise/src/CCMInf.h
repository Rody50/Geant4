/* Header file for the CCM class (CCM.h) */
#ifndef CCMINF_H
#define CCMINF_H

#include <armadillo>
#include "FillMatrix.h"
#include "InfiniteMatterSP.h"

class CCMInf
{
	private:
		mat fHamil;
		mat fVpppp, fVpphh, fVhhhh, fVhhpp;
		mat fTpphh;
		mat fInterm;
		int fNu, fA, fN;
		InfiniteMatterSP * fInfSP;
		P_group_t * fP;

	public:		
		CCMInf(P_group_t *p, InfiniteMatterSP *infSP);

		double FSum(pair_t *t, pair_t *s);

		void Interm();

		InfiniteMatterSP * GetfInfSP(){return fInfSP;}

		void ComputeH();

		double SolveT(); // return correlation energy

		void Print();

};

#endif  // CCMINF_H