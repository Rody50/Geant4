
#include "Tensor2.h"

Tensor2::Tensor2(int length) : fLength(length), fTensor(length, 0.0)
{

}

double Tensor2::operator[][](int i, int j)
{
	if (i == j) return fTensor[i];

	return 0;
}

void Tensor2::Print()
{

}