
#include "CCM.h"

CCM::CCM(int A, int N) : fA(A), fN(N), fNu(N - A), fHamil(fNu, fNu, fA, fA)
{
	fVpppp(fNu, fNu, fNu, fNu);
	fVpphh(fNu, fNu, fA, fA);
	fVhhhh(fA, fA, fA, fA);

	fTpppp(fNu, fNu, fNu, fNu);
	fTpphh(fNu, fNu, fA, fA);
	fThhhh(fA, fA, fA, fA);

	fFpp(fNu, fNu);
	fFhh(fA, fA);

	FillV(fVpppp); FillV(fVpphh); FillV(fVhhhh);
	FillF(fFpp, 0); FillF(fFhh, 1);
	FillT(fTpppp, fVpppp, fFpp, fFhh);
	FillT(fTpphh, fVpphh, fFpp, fFhh);
	FillT(fThhhh, fVhhhh, fFpp, fFhh);
}

inline void Permute(int * i1, int * i2, int ind0, int ind1)
{
	for(int i = 0; i < 4; i++)
	{
		if(i1[i] == ind0) i1[i] = ind1;
		else if(i1[i] == ind1) i1[i] = ind0;
		
		if(i2[i] == ind0) i2[i] = ind1;
		else if(i2[i] == ind1) i2[i] = ind0;
	}
}

void CCM::ComputeH()
{
	int i1[4], i2[4], i3[4], i4[4], i5[2];

	for (int a = 0; a < fNu; a++)
		for (int b = 0; b < fNu; b++)
			for (int i = 0; i < fA; i++)
				for (int j = 0; j < fA; j++)
					// term 1, 2 and 3
					fHamil[a][b][i][j] += fVpphh[a][b][i][j]
						+ fFpp[b][b] * fTpphh[a][b][i][j] - fFpp[a][a] * fTpphh[b][a][i][j] 
						- fFpp[j][j] * fTpphh[a][b][i][j] + fFpp[i][i] * fTpphh[a][b][j][i];

					// term 4 and 5
					i1[0] = 0; i1[1] = 1; i1[2] = 4; i1[3] = 5;
					i2[0] = 4; i2[1] = 5; i2[2] = 2; i2[3] = 3; 
					fHamil[a][b][i][j] +=  0.5 * fVpppp.Prod4x4(fTpphh, i1 , i2)
						+ 0.5 * fVhhhh.Prod4x4(fTpphh, i2, i1)
					
					// term 6 = 0
					// term 7
					i1[0] = 4; i1[1] = 1; i1[2] = 5; i1[3] = 3;
					i2[0] = 0; i2[1] = 5; i2[2] = 2; i2[3] = 4;
					i3[0] = 0; i3[1] = 4; i3[2] = 2; i3[3] = 5;
					i4[0] = 5; i4[1] = 1; i4[2] = 4; i4[3] = 3; // a b i j
					fHamil[a][b][i][j] +=
						 0.5 * fVpppp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4);
					
					Permute(i1, i2, 0, 1); // - b a i j
					fHamil[a][b][i][j] +=
						- 0.5 * fVpppp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4);
					
					Permute(i3, i4, 2, 3); // b a j i
					fHamil[a][b][i][j] +=
						 0.5 * fVpppp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4);
					
					Permute(i1, i2, 0, 1); // - a b j i
					fHamil[a][b][i][j] +=
						- 0.5 * fVpppp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4);

					// term 8
					i1[0] = 2; i1[1] = 0; i1[2] = 3; i1[3] = 4;
					i2[0] = 3; i2[1] = 4; i2[2] = 1; i2[3] = 2;
					i3[0] = 0; i3[1] = 1; i3[2] = 4; i3[3] = 3;
					i5[0] = 4; i5[1] = 2; // a b i j
					fHamil[a][b][i][j] +=
						0.5 * fTpphh.Prod4x2(fVpppp.Prod4x4(fTpphh, i1 , i2), i3, i5);

					// term 9

					// term 10
}

void CCM::SolveT()
{

}