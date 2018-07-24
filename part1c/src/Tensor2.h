/* Header file of the tensor class */
#ifndef TENSOR2_H
#define TENSOR2_H

#include "Tensor.h"

class Tensor2 : Tensor
{
	private:
		int fLength;

	public:
		Tensor2(int length);

		double operator[][](int i, int j);

		virtual void Print();

}

#endif // TENSOR2_H