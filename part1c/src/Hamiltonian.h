/* Header file for exercise 1c (Hamiltonian.h) */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <armadillo>
#include "TGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TAxis.h"

using namespace arma;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;

const int n = 6;

double eigenvalue(double g, double d);
void Draw(const char *txtfile); // show the final result in graphical form
