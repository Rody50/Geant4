/* Code that Diagonalizes the Hamiltonian -- the main function (main.cpp) */

#include "Hamiltonian.h"

extern const int n;

int main(int argc, char *argv []){
	const double d = 1.0;
	const int nn = 21; // the number of the correlation energy sampling
	double g_vec[nn], cor_e_vec[nn];
	
	char txtfile[] = "results/correlation_value.txt";
	ofstream file(txtfile);

	for(int i = 0; i < nn; i++){
		g_vec[i] = i / 10. - 1;
		cor_e_vec[i] = eigenvalue(g_vec[i], d) - (2 * d - g_vec[i]);
		cout << "The first correlation energy for g = ";
		cout << g_vec[i] << " is:" <<  cor_e_vec[i] << endl;
		file << g_vec[i] << "\t" << cor_e_vec[i] << endl;
	} // end loop over i
	file.close();

	Draw(txtfile);
	system((string("evince ")+txtfile+".pdf").c_str());
} // end of the main function
