/* Header file of the tensor class */
#ifndef TENSOR4_H
#define TENSOR4_H

#include "Tensor.h"

class Tensor2;

class Tensor4 : public Tensor
{
	private: 
		int fL1, fL2, fL3, fL4;

	public:
		Tensor4(int l1, int l2, int l3, int l4);

		double & operator()(int i, int j, int a, int b);
		const double & operator()(int i, int j, int a, int b) const{
			return operator()(i, j, a, b);
		}

		int GetL(int i) const;

		Tensor4 Prod4x4(Tensor4 & t, int * iO1, int * iO2);

		Tensor2 Prod4x4s3(Tensor4 & t, int * iO1, int * iO2);

		Tensor4 Prod4x2(const Tensor2 & t, int * iO1, int * iO2);

		virtual void Print();


};

#endif // TENSOR4_H
