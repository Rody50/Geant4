/* Header file of the tensor class */
#ifndef TENSOR4_H
#define TENSOR4_H

#include "Tensor.h"

class Tensor2;

class Tensor4 : Tensor
{
	private: 
		int fL1, fL2, fL3, fL4;

	public:
		Tensor4(int l1, int l2, int l3, int l4);

		double operator[][][][](int i, int j, int a, int b);

		int GetL(int i);

		void Multiplication4x4(const Tensor4 & t, Tensor4 & result);

		void Multiplication4x2(const Tensor2 & t, Tensor2 & result);

		virtual void Print();


}

#endif // TENSOR4_H