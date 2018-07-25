
#include "FillMatrix.h"
#include "Tensor2.h"
#include "Tensor4.h"


FillMatrix::FillMatrix()
{}

void FillMatrix::FillV(Tensor4 & v)
{
	int * iMax[4];

	for (int i = 0; i < 4; i++)
		iMax[i] = v.GetL(i);

	for (int a = 0; a < iMax[0]; a++)
		for (int b = 0; b < iMax[1]; b++)
			for (int i = 0; i < iMax[2]; i++)
				for (int j = 0; j < iMax[3]; j++)
				{
					if (IsPair(i, j) && IsPair(a, b) && i < j && a < b )
						{
							v[a][b][i][j] = - g / 2.;
							v[a][b][j][i] =   g / 2.;
							v[b][a][i][j] =   g / 2.;
							v[b][a][j][i] = - g / 2.;
						}	
				}
}

void FillMatrix::FillF(Tensor2 & f, bool isHole)
{

	for (int i = 0; i < f.GetL(1); i++)
	{
		f[i][i] = d * floor(i / 2.);
		if (isHole) f[i][i] += - g / 2.;
	}
}

void FillMatrix::FillT(Tensor4 & t, Tensor4 & v, Tensor2 & fp,
	Tensor2 & fh)
{
	int * iMax[4];

	for (int i = 0; i < 4; i++)
		iMax[i] = t.GetL(i);

	for (int a = 0; a < iMax[0]; a++)
		for (int b = 0; b < iMax[1]; b++)
			for (int i = 0; i < iMax[2]; i++)
				for (int j = 0; j < iMax[3]; j++)
				{
					if (IsPair(i, j) && IsPair(a, b) && i < j && a < b)
						{
							t[a][b][i][j] = v[a][b][i][j] / 
								(fh[i][i] + fh[j][j] - fp[a][a] - fp[b][b]);
							t[a][b][j][i] = - t[a][b][i][j];
							t[b][a][i][j] = - t[a][b][i][j];
							t[b][a][j][i] =   t[a][b][i][j];
						}

				}
}

void FillMatrix::Print()
{
	fConfigObj.Print();
	v->print();
}

bool FillMatrix::IsPair(int a, int b)
{
	return floor(a / 2.) == floor(b / 2.);
}