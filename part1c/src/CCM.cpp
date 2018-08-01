
#include "CCM.h"
#include <fstream>
#include <string>

using std::ofstream;
using std::to_string;
using std::flush;

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

	// cout << "The Vpphh matrix: " << endl;
	// fVpphh.Print();
	// cout << "The Fpp matrix: " << endl;
	// fFpp.Print();
	// cout << "The Fhh matrix: " << endl;
	// fFhh.Print();
	// cout << "The Tpphh matrix " << endl;
	// fTpphh.Print();
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
		for (int i = 0; i < fA; i++)
		{	
			int b = a + 1;
			int j = i + 1;

			if (!(a / 2 == b / 2 && i / 2 ==  j / 2)) continue;
			// term 1, 2 and 3
			fHamil(a, b, i, j) = fVpphh(a, b, i, j)
				+ fFpp(b, b) * fTpphh(a, b, i, j) - fFpp(a, a) * fTpphh(b, a, i, j)
				- fFhh(j, j) * fTpphh(a, b, i, j) + fFhh(i, i) * fTpphh(a, b, j, i);

			// term 4 and 5
			i1[0] = 0; i1[1] = 1; i1[2] = 4; i1[3] = 5;
			i2[0] = 4; i2[1] = 5; i2[2] = 2; i2[3] = 3; 
			fHamil(a, b, i, j) +=  0.5 * fVpppp.Prod4x4(fTpphh, i1 , i2)(a, b, i, j)
				+ 0.5 * fVhhhh.Prod4x4(fTpphh, i2, i1)(a, b, i, j);

			// // term 7
			i1[0] = 4; i1[1] = 1; i1[2] = 5; i1[3] = 3;
			i2[0] = 0; i2[1] = 5; i2[2] = 2; i2[3] = 4;
			i3[0] = 0; i3[1] = 4; i3[2] = 2; i3[3] = 5;
			i4[0] = 5; i4[1] = 1; i4[2] = 4; i4[3] = 3; 
			Tensor4 term7(fNu, fNu, fA, fA);
			term7 = fVhhpp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4);
			// a b i j
			fHamil(a, b, i, j) += 0.5 * term7(a, b, i, j);
			
			// - b a i j
			fHamil(a, b, i, j) += - 0.5 * term7(b, a, i, j);

			// b a j i
			fHamil(a, b, i, j) += 0.5 * term7(b, a, j, i);
			
			// - a b j i
			fHamil(a, b, i, j) += - 0.5 * term7(a, b, j, i);
			
			// term 8
			i1[0] = 2; i1[1] = 0; i1[2] = 3; i1[3] = 4;
			i2[0] = 3; i2[1] = 4; i2[2] = 1; i2[3] = 2;
			i3[0] = 0; i3[1] = 1; i3[2] = 4; i3[3] = 3;
			i5[0] = 4; i5[1] = 2; i5[2] = 9; i5[3] = 9; 
			Tensor4 term8(fNu, fNu, fA, fA);
			term8 = fTpphh.Prod4x2(fVhhpp.Prod4x4s3(fTpphh, i1 , i2), i3, i5);
			
			// a b i j
			fHamil(a, b, i, j) += 0.5 * term8(a, b, i, j);

			// a b j i
			fHamil(a, b, i, j) += - 0.5 * term8(a, b, j, i);

			// term 9
			i1[0] = 2; i1[1] = 3; i1[2] = 4; i1[3] = 1;
			i2[0] = 0; i2[1] = 4; i2[2] = 2; i2[3] = 3;
			i3[0] = 4; i3[1] = 1; i3[2] = 2; i3[3] = 3;
			i5[0] = 0; i5[1] = 4; i5[2] = 9; i5[3] = 9; 
			Tensor4 term9(fNu, fNu, fA, fA);
			term9 = fTpphh.Prod4x2(fVhhpp.Prod4x4s3(fTpphh, i1 , i2), i3, i5);
			
			// a b i j
			fHamil(a, b, i, j) += 0.5 * term9(a, b, i, j);

			// b a i j
			fHamil(a, b, i, j) += - 0.5 * term9(b, a, i, j);					

			// term 10
			i1[0] = 0; i1[1] = 1; i1[2] = 4; i1[3] = 5;
			i2[0] = 4; i2[1] = 5; i2[2] = 2; i2[3] = 3;
			i3[0] = 4; i3[1] = 5; i3[2] = 2; i3[3] = 3;
			i4[0] = 0; i4[1] = 1; i4[2] = 4; i4[3] = 5; // a b i j
			fHamil(a, b, i, j) +=
				 0.25 * fVhhpp.Prod4x4(fTpphh, i1 , i2).Prod4x4(fTpphh, i3, i4)(a, b, i, j);

			// antisymmetry
			fHamil(b, a, i, j) = - fHamil(a, b, i, j);
			fHamil(a, b, j, i) = - fHamil(a, b, i, j);
			fHamil(b, a, j, i) = fHamil(a, b, i, j);
		}
}

double CCM::SolveT()
{
	double factor = 1.;

	double corr_en_pre = 1.;
	double corr_en = 0.;
	vector<double> correlationEn;

	int NStepsMax = 1E5;
	int NSteps = 0;

	ofstream corr_en_file(("results/corr_en-steps" + to_string(fN) 
		+ "_" + to_string(fA) + ".txt").c_str());


	while (fabs(corr_en - corr_en_pre) > 1E-5 || NSteps > NStepsMax)
	{	
		ComputeH();
		
		corr_en_pre = corr_en;
		corr_en = 0.;

		for (int a = 0; a < fNu; a++)
			for (int i = 0; i < fA; i++)
			{	
				int b = a + 1;
				int j = i + 1;
				if (!(a / 2 == b / 2 && i / 2 == j /2)) continue;
				
				// a b i j
				corr_en += 0.25 * fVhhpp(i, j, a, b)
					* fTpphh(a, b, i, j);

				// b a i j 
				corr_en += 0.25 * fVhhpp(i, j, b, a) 
					* fTpphh(b, a, i, j);

				// a b j i 
				corr_en += 0.25 * fVhhpp(j, i, a, b)
					* fTpphh(a, b, j, i);

				// b a j i
				corr_en += 0.25 * fVhhpp(j, i, b, a)
					* fTpphh(b, a, j, i);

				double deno = fFhh(i, i) + fFhh(j, j) - fFpp(a, a) - fFpp(b, b);
				double temp = fHamil(a, b, i, j) / deno;	

					fTpphh(a, b, i, j) += factor * temp;
			
				cout << "a: " << a;
				cout << "\tcorr_en: " << corr_en;
				cout << " \tNSteps: " << NSteps;
				cout << "\r" << flush;

				fTpphh(b, a, i, j) = - fTpphh(a, b, i, j);
				fTpphh(a, b, j, i) = - fTpphh(a, b, i, j);
				fTpphh(b, a, j, i) = fTpphh(a, b, i, j);

			}

		correlationEn.push_back(corr_en);
		NSteps++;
		corr_en_file << corr_en << endl;
	}
	corr_en_file.close();

	cout << "The converged correlation energy is: " << corr_en << endl;
	
	return corr_en;
}

void CCM::Print()
{
	fHamil.Print();
}
