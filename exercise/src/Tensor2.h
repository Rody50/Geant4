/* Header file of the tensor class */
#ifndef TENSOR2_H
#define TENSOR2_H

#include "Tensor.h"
#include <iostream>

class Tensor2 : public Tensor
{
	private:
		int fL1, fL2;

	public:
		Tensor2(int l1, int l2);

		double & operator()(int i, int j);
		
		const double & operator()(int i, int j) const
		{	
			return fTensor[i * fL1 + j];
		}

		int GetL(int i) const;

		virtual void Print();

};

#endif // TENSOR2_H
