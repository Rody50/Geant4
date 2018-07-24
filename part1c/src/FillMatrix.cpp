
#include "FillMatrix.h"

FillMatrix::FillMatrix(int N, int particles) : fA(particles), 
	fN(N), fNu(N-particles)
{

}

void FillMatrix::FillV()
{

	for (int i = 0; i < fN; i++)
		for (int j = 0; j < fN; j++)
			for (int a = 0; a < fN; a++)
				for (int b = 0; b < fN; b++)
				{
					if (IsPair(i, j) && IsPair(a, b) && i < j && a < b )
						{
							fV[i, j, a, b] = - g / 2.;
							fV[i, j, b, a] =   g / 2.;
							fV[j, i, a, b] =   g / 2.;
							fV[j, i, b, a] = - g / 2.;
						}	
				}
}

void FillMatrix::FillF()
{

	for (int i = 0; i < fN; i++)
	{
			fF[i][i] = d * floor(i / 2.);
		if (i <= fA)
			fF[i][i] += - g / 2.;
	}

}

void FillMatrix::FillT()
{

	for (int i = 0; i < fN; i++)
		for (int j = 0; j < fN; j++)
			for (int a = 0; a < fN; a++)
				for (int b = 0; b < fN; b++)
				{
					if (IsPair(i, j) && IsPair(a, b) && i < j && a < b )
						{
							fT[i, j, a, b] = fV[i, j, a, b] / 
								(fF[a][a] + fF[b][b] - fF[i][i] - fF[j][j]);
							fT[i, j, b, a] = - fT[i, j, a, b];
							fT[j, i, a, b] = - fT[i, j, a, b];
							fT[j, i, b, a] =   fT[i, j, a, b];
						}

				}
}

void FillMatrix::Print()
{
	fConfigObj.Print();
	fV->print();
}

bool FillMatrix::IsPair(int a, int b)
{
	return floor(a / 2.) == floor(b / 2.);
}