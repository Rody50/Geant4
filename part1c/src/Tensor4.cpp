
#include <iostream>
#include "Tensor4.h"
#include "Tensor2.h"


Tensor4::Tensor4(int l1, int l2, int l3, int l4)
	 : Tensor(l1 * l2 * l3 * l4, 0.0),
	 fL1(l1), fL2(l2), fL3(l3), fL4(l4)
{}

double & Tensor4::operator()(int a, int b, int i, int j)
{
	return fTensor[ fL1 * fL2 * (fL3 * a + b) + fL1 * i + j];
}

int Tensor4::GetL(int i) const
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

Tensor4 Tensor4::Prod4x4(Tensor4 & t, 
	int * iO1, int * iO2)
{
	int iMax[6];
	for (int i = 0; i < 4; i++)
		for(int j = 0; j < 6; j++)
		{
			if(iO1[i] == j) iMax[j] = GetL(i);
			if(iO2[i] == j) iMax[j] = t.GetL(i);
		}

	Tensor4 result(iMax[0], iMax[1], iMax[2], iMax[3]);

	for (int i1 = 0; i1 < iMax[0]; i1++)
		for (int i2 = 0; i2 < iMax[1]; i2++)
			for (int i3 = 0; i3 < iMax[2]; i3++)
				for (int i4 = 0; i4 < iMax[3]; i4++)
				{
					double result_temp = 0.; 

					for (int i = 0; i < iMax[4]; i++)
						for (int j = 0; j < iMax[5]; j++)
						{								
							int iv[4], it[4];
							int id[6] = {i1, i2, i3, i4, i, j};

							for (int ii = 0; ii < 4; ii++)
							{	
								iv[ii] = id[iO1[ii]];
								it[ii] = id[iO2[ii]];
							}
							 
							result_temp +=
								(*this)(iv[0], iv[1], iv[2], iv[3]) 
								* t(it[0], it[1], it[2], it[3]);
						}
						result(i1, i2, i3, i4) = result_temp;
				}

	return result;

}

Tensor2 Tensor4::Prod4x4s3(Tensor4 & t,
 int * iO1, int * iO2)
{
	int iMax[5];
	for (int i = 0; i < 4; i++)
		for(int j = 0; j < 5; j++)
		{
			if(iO1[i] == j) iMax[j] = GetL(i);
			if(iO2[i] == j) iMax[j] = t.GetL(i);
		}

	Tensor2 result(iMax[0], iMax[1]);

	for (int i1 = 0; i1 < iMax[0]; i1++)
		for (int i2 = 0; i2 < iMax[1]; i2++)
		{
			double result_temp = 0.; 

			for (int i = 0; i < iMax[2]; i++)
				for (int j = 0; j < iMax[3]; j++)
					for (int k = 0; k < iMax[4]; k++)
					{
						int iv[4], it[4];
						int id[5] = {i1, i2, i, j, k};

						for (int ii = 0; ii < 4; ii++)
						{	
							iv[ii] = id[iO1[ii]];
							it[ii] = id[iO2[ii]];
						}
						 
						result_temp +=
							(*this)(iv[0], iv[1], iv[2], iv[3]) 
							* t(it[0], it[1], it[2], it[3]);
					}

			result(i1, i2) = result_temp;
		}
	return result;
}

Tensor4 Tensor4::Prod4x2(const Tensor2 & t, 
	int * iO1, int * iO2)
{
	int iMax[5];
	for (int i = 0; i < 4; i++)
		for(int j = 0; j < 5; j++)
		{
			if(iO1[i] == j) iMax[j] = GetL(i);
			if(iO2[i] == j && i < 2) iMax[j] = t.GetL(i);
		}

	Tensor4 result(iMax[0], iMax[1], iMax[2], iMax[3]);

	for (int i1 = 0; i1 < iMax[0]; i1++)
		for (int i2 = 0; i2 < iMax[1]; i2++)
			for (int i3 = 0; i3 < iMax[2]; i3++)
				for (int i4 = 0; i4 < iMax[3]; i4++)
				{
					double result_temp = 0.; 
					for (int i = 0; i < iMax[4]; i++)
					{								
						int iv[4], it[2];
						int id[5] = {i1, i2, i3, i4, i};

						for (int ii = 0; ii < 4; ii++)
						{	
							iv[ii] = id[iO1[ii]];
							if (ii < 2) it[ii] = id[iO2[ii]];
						}
						result_temp +=
							(*this)(iv[0], iv[1], iv[2], iv[3]) 
							* t(it[0], it[1]);
					}

					result(i1, i2, i3, i4) = result_temp;
				}

	return result;
}

void Tensor4::Print()
{
	//for (int a = 0; a < fL1; a++)
	//	for (int b = 0; b < fL2; b++)
		{
			for (int i = 0; i < fL3; i++)
				for (int j = 0; j < fL4; j++)
				{
					std::cout << (*this)(0, 1, i, j) << " \t";

				}

			std::cout << std::endl;	
		}		
}
