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
		Tensor();

		virtual void Print();

}

#endif // TENSOR_H