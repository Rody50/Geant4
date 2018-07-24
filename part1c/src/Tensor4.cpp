
#include "Tensor4.h"

Tensor4::Tensor4(int l1, int l2) : fL1(l1), fL2(l2),
 fTensor(l1 * l1 * l2 * l2, 0.0)
{}

double Tensor4::operator[][][][](int i, int j, int a, int b)
{
	return fTensor[((i - 1) * l1 + j - 1) * l2 * l2 + (a - 1) * l2 + b];
}
void Tensor4::Multiplication4x4(const Tensor4 & t, Tensor4 & result, 
	int sum_index1, int sum_index2, int sum_index3, int option)
{
	if (option == -1)
	{
		for (int i1 = 0; i1 < ; i1++)
			for (int i2 = 0; i2 < ; i2++)
				for (int i3 = 0; i3 < ; i3++)
					for (int i4 = 0; i4 < ; i4++)
					{
						double result_temp = 0.; 
							for (int i = 0; i < ; i++)
								for (int j = 0; j < ; j++)
									result_temp += fTensor[][][][] * t[][][][];
						
						result[i1][i2][i3][i4] = result_temp
					}
}

void Tensor4::Multiplication4x2(const Tensor2 & t, int i, Tensor2 & result)
{

}

void Tensor4::Print()
{

}