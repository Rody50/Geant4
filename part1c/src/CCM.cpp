
#include "CCM.h"

CCM::CCM(int A, int N, double g, double d) : fA(A), fN(N),
	fNu(N - A), fG(g), fD(d), fHamil(N - A, N - A, A, A),
	fVpppp(N - A, N - A, N - A, N - A),	
	fVpphh(N - A, N - A, A, A), 
	fVhhpp(A, A, N - A, N - A), 
	fVhhhh(A, A, A, A),
	fTpphh(N - A, N - A, A, A),
	fFpp(N - A, N - A),
	fFhh(A, A)
{
	FillMatrix::FillV(fVpppp, fG); FillMatrix::FillV(fVpphh, fG);
	FillMatrix::FillV(fVhhhh, fG); FillMatrix::FillV(fVhhpp, fG);
	FillMatrix::FillF(fFpp, fA, fG, fD); FillMatrix::FillF(fFhh, 0, fG, fD);
	FillMatrix::FillT(fTpphh, fVpphh, fFpp, fFhh);

	cout << "The Vpphh matrix: " << endl;
	fVpphh.Print();
	cout << "The Fpp matrix: " << endl;
	fFpp.Print();
	cout << "The Fhh matrix: " << endl;
	fFhh.Print();
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
	int i1[4], i2[4], i3[4], i4[4], i5[4];

	for (int a = 0; a < fNu; a++)
		for (int b = 0; b < fNu; b++)
			for (int i = 0; i < fA; i++)
				for (int j = 0; j < fA; j++)
				{	
					cout << "Did we enter the first calcuation? " << endl;
					// term 1, 2 and 3
					fHamil(a, b, i, j) += fVpphh(a, b, i, j)
						+ fFpp(b, b) * fTpphh(a, b, i, j) - fFpp(a, a) * fTpphh(b, a, i, j) 
						- fFpp(j, j) * fTpphh(a, b, i, j) + fFpp(i, i) * fTpphh(a, b, j, i);
					//getchar();
					//cout << "#1" << endl;
					// term 4 and 5
					i1[0] = 0; i1[1] = 1; i1[2] = 4; i1[3] = 5;
					i2[0] = 4; i2[1] = 5; i2[2] = 2; i2[3] = 3; 
					fHamil(a, b, i, j) +=  0.5 * fVpppp.Prod4x4(fTpphh, i1 , i2)(a, b, i, j)
						+ 0.5 * fVhhhh.Prod4x4(fTpphh, i2, i1)(a, b, i, j);
					//getchar();
					//cout << "#2" << endl;
					// term 6 = 0

					// term 7
					i1[0] = 4; i1[1] = 1; i1[2] = 5; i1[3] = 3;
					i2[0] = 0; i2[1] = 5; i2[2] = 2; i2[3] = 4;
					i3[0] = 0; i3[1] = 4; i3[2] = 2; i3[3] = 5;
					i4[0] = 5; i4[1] = 1; i4[2] = 4; i4[3] = 3; // a b i j
					fHamil(a, b, i, j) +=
						 0.5 * fVhhpp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4)(a, b, i, j);
					
					Permute(i1, i2, 0, 1); // - b a i j
					fHamil(a, b, i, j) +=
						- 0.5 * fVhhpp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4)(a, b, i, j);
					
					Permute(i3, i4, 2, 3); // b a j i
					fHamil(a, b, i, j) +=
						 0.5 * fVhhpp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4)(a, b, i, j);
					
					Permute(i1, i2, 0, 1); // - a b j i
					fHamil(a, b, i, j) +=
						- 0.5 * fVhhpp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4)(a, b, i, j);
					
					getchar();
					cout << "#3" << endl;
					// term 8
					i1[0] = 2; i1[1] = 0; i1[2] = 3; i1[3] = 4;
					i2[0] = 3; i2[1] = 4; i2[2] = 1; i2[3] = 2;
					i3[0] = 0; i3[1] = 1; i3[2] = 4; i3[3] = 3;
					i5[0] = 4; i5[1] = 2; i5[2] = 9; i5[3] = 9; // a b i j
					
					
					Tensor2 intermediate = fVhhpp.Prod4x4s3(fTpphh, i1 , i2);
					 fHamil(a, b, i, j) +=
					 	0.5 * fTpphh.Prod4x2(intermediate, i3, i5)(a, b, i, j);
					 getchar();
					cout << "#4" << endl;
					// Permute(i2, i3, 2, 3); // a b j i
					// fHamil(a, b, i, j) +=
					// 	- 0.5 * fTpphh.Prod4x2(fVhhpp.Prod4x4s3(fTpphh, i1 , i2), i3, i5)(a, b, i, j);
					// getchar();
					// cout << "#4" << endl;
					// // term 9
					// i1[0] = 2; i1[1] = 3; i1[2] = 4; i1[3] = 1;
					// i2[0] = 0; i2[1] = 4; i2[2] = 2; i2[3] = 3;
					// i3[0] = 4; i3[1] = 1; i3[2] = 2; i3[3] = 3;
					// i5[0] = 0; i5[1] = 4; i5[2] = 9; i5[3] = 9; // a b i j
					// fHamil(a, b, i, j) +=
					// 	0.5 * fTpphh.Prod4x2(fVhhpp.Prod4x4s3(fTpphh, i1 , i2), i3, i5)(a, b, i, j);

					// Permute(i2, i3, 0, 1); // b a i j
					// fHamil(a, b, i, j) +=
					// 	- 0.5 * fTpphh.Prod4x2(fVhhpp.Prod4x4s3(fTpphh, i1 , i2), i3, i5)(a, b, i, j);					
					// getchar();
					// cout << "#5" << endl;
					// // term 10
					// i1[0] = 0; i1[1] = 1; i1[2] = 4; i1[3] = 5;
					// i2[0] = 4; i2[1] = 5; i2[2] = 2; i2[3] = 3;
					// i3[0] = 4; i3[1] = 5; i3[2] = 2; i3[3] = 3;
					// i4[0] = 0; i4[1] = 1; i4[2] = 4; i4[3] = 5; // a b i j
					// fHamil(a, b, i, j) +=
					// 	 0.25 * fVhhpp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4)(a, b, i, j);
					// getchar();
					// cout << "#6" << endl;
				}
}

void CCM::SolveT()
{
	double delta = 1.;
	double delta_pre = 0.;

	while(fabs(delta - delta_pre) > 1E-2)
	{
		cout << "Did we enter the first iteration?" << endl;
		ComputeH();
		cout << "Could we compute H?" << endl;

		for (int a = 0; a < fNu; a++)
			for (int b = 0; b < fNu; b++)
				for (int i = 0; i < fA; i++)
					for (int j = 0; j < fA; j++)
					{	
						double deno = fFhh(i, i) + fFhh(j, j) - fFpp(a, a) - fFpp(b, b);
						double tmp = fHamil(a, b, i, j) / deno;
						delta += fabs(fTpphh(a, b, i, j) - tmp);
						fTpphh(a, b, i, j) += tmp;
					}

					delta_pre = delta;
					cout << "The new T matrix is: " << endl;
					fTpphh.Print();
					cout << endl;
					cout << "The new value of delta-delta is: " 
							 << fabs(delta - delta_pre) << endl;
	}
}

void CCM::Print()
{
	fHamil.Print();
}
