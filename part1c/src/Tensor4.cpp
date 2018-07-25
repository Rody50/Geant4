
#include "Tensor4.h"

Tensor4::Tensor4(int l1, int l2, int l3, int l4) : fL1(l1), fL2(l2),
 fL3(l3), fL4(l4), fTensor(l1 * l1 * l2 * l2, 0.0)
{}

double Tensor4::operator[][][][](int i, int j, int a, int b)
{
	return fTensor[((i - 1) * l1 + j - 1) * l2 * l2 + (a - 1) * l2 + b];
}

int Tensor4::GetL(int i)
{
	switch(i)
	{
		case 0: return fL1;
		case 1: return fL2;
		case 2: return fL3;
		case 3: return fL4;
		default: break;
	}
	return -1;
}

void Tensor4::Multiplication4x4(const Tensor4 & t, Tensor4 & result, 
	int * iO1, int * iO2, int option)
{
	int iMax[6];
	for (int i = 0; i < 6; i++)
		for(int j = 0; j < 6; j++)
		{
			if(iO1[i] == j) iMax[i] = GetL(i);
			if(iO2[i] == j) iMax[i] = t.GetL(i);
		}

	if (option == -1)
	{
		for (int i1 = 0; i1 < iMax[0]; i1++)
			for (int i2 = 0; i2 < iMax[1]; i2++)
				for (int i3 = 0; i3 < iMax[2]; i3++)
					for (int i4 = 0; i4 < iMax[3]; i4++)
					{
						double result_temp = 0.; 
						for (int i = 0; i < iMax[5]; i++)
							for (int j = 0; j < iMax[6]; j++)
							{
								int iv[4], it[4];
								int id[6] = {i1, i2, i3, i4, i, j};

								for (int ii = 0; ii < 4; ii++)
								{	
									iv[ii] = id[iO1[ii]];
									it[ii] = id[iO2[ii]];
								}
								 
								result_temp +=
									fTensor[iv[0]][iv[1]][iv[2]][iv[3]] 
									* t[it[0]][it[1]][it[2]][it[3]];
							}
						result[i1][i2][i3][i4] = result_temp
					}
	}
}

void Tensor4::Multiplication4x2(const Tensor2 & t, int i, Tensor2 & result)
{

}

void Tensor4::Print()
{

}