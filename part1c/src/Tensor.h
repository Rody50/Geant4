/* Header file of the tensor class */
#ifndef TENSOR_H
#define TENSOR_H

#include <vector>

using std::vector;

class Tensor
{
	protected:
		vector<double> fTensor;

	public:
		Tensor(int l, double val = 0.);

		virtual void Print() = 0;

};

#endif // TENSOR_H
