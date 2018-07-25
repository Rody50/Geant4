
#include "Tensor2.h"

Tensor2::Tensor2(int l1, int l2) : fL1(l1), fL2(l2), fTensor(l1 * l2, 0.0)
{}

double Tensor2::operator[][](int i, int j)
{
	return fTensor[i * fL1 + j];
}

int Tensor2::GetL(int i)
{
	switch(i)
	{
		case 0: return fL1;
		case 1: return fL2;
		default: break;
	}
	return -1;
}

void Tensor2::Print()
{

}