/* Code that Diagonalizes the Hamiltonian (Diagonalization.cpp) */
#include "Diagonalization.h"

double eigenValue(const double g, const double d, const double * matrix, const int length)
{
	mat hamiltonian(length, length);

	for (int index = 0; index < length * length; index++)
		hamiltonian(index) = matrix[index];

	vec eigval(length);

	eig_sym(eigval, hamiltonian);

	return eigval[0];
}