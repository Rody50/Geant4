#include "Diagonalization.h"
#include "SingleParticleState.h"

int main()
{
	const int n = 6;

	const double d = 1.0;
	double g_vec[21];
	double corr_vec[21];

	for (int i = 0; i < 21; i++)
	{
		g_vec[i] = (i - 10) / 10.;
		double g = g_vec[i];
		double matrix[n * n] = 
		 	{2 * d - g, - g/2.0, - g/2.0, - g/2.0, - g/2.0, 0.0,
			- g/2.0, 4 * d - g, - g/2.0, - g/2.0, 0.0, - g/2.0,
			- g/2.0, - g/2.0, 6 * d - g, 0.0, - g/2.0, - g/2.0,
			- g/2.0, - g/2.0, 0.0, 6 * d - g, - g/2.0, - g/2.0,
			- g/2.0, 0.0, - g/2.0, - g/2.0, 8 * d - g, - g/2.0,
			0.0, - g/2.0, - g/2.0, - g/2.0, - g/2.0, 10 * d - g};
		corr_vec[i] = eigenValue(g_vec[i], d, matrix, n) - (2 * d - g_vec[i]);
		cout << "The correlation energy is: " << corr_vec[i] << std::endl; 

	}

//Write data to file 
    string fileName = "./results/correlation_energy.txt";

    FILE *formfactorFILE;
    formfactorFILE = fopen(fileName.c_str(), "w");

    for (int count = 0; count <= 21; count ++) 
        fprintf(formfactorFILE, "%lf \t %lf \n", g_vec[count], corr_vec[count]);

    fclose(formfactorFILE);

  SingleParticleState sPS(2);

	sPS.GenerateSPS();
	sPS.print();

}