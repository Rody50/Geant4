/* Header file for the CCM class (CCM.h) */
#ifndef IMSRGINF_H
#define IMSRGINF_H

#include <armadillo>

#include "FillMatrix.h"
#include "InfiniteMatterSP.h"

class IMSRGInf
{
	private:
		mat fHamil0, fHamilS;
		mat fF, fV;
		mat fEta, fOmega;
		int fNu, fA;
		double fDeltaS;
		InfiniteMatterSP * fInfSP;
		P_group_t * fP;

	public:		
		IMSRGInf(P_group_t *p, InfiniteMatterSP *infSP, double deltaS);

		double F(qState *s);

		double FSum(pair_t *t, pair_t *s);

		InfiniteMatterSP * GetfInfSP(){return fInfSP;}
		
		mat Commute(const mat &a, const mat &b);

		void UpdateOm();

		void UpdateEta();

		double Iterate(); // return correlation energy

		//void Print();

};

#endif  // IMSRGINF_H