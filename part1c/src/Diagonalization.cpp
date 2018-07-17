/* Code that Diagonalizes the Hamiltonian (Diagonalization.cpp) */

#include "Hamiltonian.h"

extern const int n;

double eigenvalue(double g, double d){
	double matrix[n * n] = {
	 2 * d - g,	  - g/2.0,   - g/2.0,   - g/2.0,   - g/2.0,        0.0,
	   - g/2.0, 4 * d - g,   - g/2.0,   - g/2.0,       0.0,    - g/2.0,
	   - g/2.0,   - g/2.0, 6 * d - g,       0.0,   - g/2.0,    - g/2.0,
	   - g/2.0,   - g/2.0,       0.0, 6 * d - g,   - g/2.0,    - g/2.0,
	   - g/2.0,       0.0,   - g/2.0,   - g/2.0, 8 * d - g,    - g/2.0,
	       0.0,   - g/2.0,   - g/2.0,   - g/2.0,   - g/2.0, 10 * d - g
	};

	mat hamiltonian(n, n);
	for (int index = 0; index < n * n; index++){
		hamiltonian(index) = matrix[index];
	}

	vec Eigval(n);
	eig_sym(Eigval, hamiltonian);
//	cout << "The first eigenvalue is: " << Eigval[0] << endl;

	return Eigval[0];
}

void Draw(const char *txtfile){
	TGraph *g = new TGraph(txtfile, "%lg %lg");
	g->SetNameTitle("g", "Correlation Energy in Paring Model;g/d;binding energy");
	TCanvas *c = new TCanvas("c", "Binding Energy in Paring Model", 800, 600);
	g->SetLineWidth(2); g->SetLineColor(2); g->SetMarkerStyle(20);
	g->GetYaxis()->SetTitleOffset(1.3);
	g->Draw("APL");
}


