/* Code that Diagonalizes the Hamiltonian (Diagonalization.cpp) */

#include "Hamiltonian.h"


int main()
{

	const double d = 1.0;
	const double g = -1.0;

 double matrix[n * n] = 
 			{2 * d - g/2.0, - g/2.0, - g/2.0, - g/2.0, - g/2.0, 0.0,
			- g/2.0, 4 * d - g/2.0, - g/2.0, - g/2.0, 0.0, - g/2.0,
			- g/2.0, - g/2.0, 6 * d - g/2.0, 0.0, - g/2.0, - g/2.0,
			- g/2.0, - g/2.0, 0.0, 6 * d - g/2.0, - g/2.0, - g/2.0,
			- g/2.0, 0.0, - g/2.0, - g/2.0, 8 * d - g/2.0, - g/2.0,
			0.0, - g/2.0, - g/2.0, - g/2.0, - g/2.0, 10 * d - g/2.0};


for (int index = 0; index < n * n; index++)
	hamiltonian(index) = matrix[index];

hamiltonian.print();

vec Eigval(n);
eig_sym(Eigval, hamiltonian);
cout << "The first eigenvalue is: " << Eigval[0] << std::endl;
}