/* Implementation of the SingleParticleState class (SingleParticleState.cpp)*/

#include "InfiniteMatterSP.h"

InfiniteMatterSP::InfiniteMatterSP(int Nmax, int nShell, double rho) : 
	fNmax(Nmax), fNShell(nShell), fRho(rho)
{}

void InfiniteMatterSP::GenerateSP()
{
	qState temp;
	// double L = pow(fNmagic / fRho, 1 / 3.);

	for (int nx = -fNmax; nx <= fNmax; nx++)
	{
		temp.nx = nx;
		for (int ny = -fNmax; ny <= fNmax; ny++)
		{
			temp.ny = ny;
			for (int nz = -fNmax; nz <= fNmax; nz++)
			{
				temp.nz = nz;
				for (int twoS = -1; twoS < 2; twoS+=2)
				{
					temp.two_s = twoS;
					fStates.push_back(temp);
				}
			}
		}
	}	
}

double Radius(qState t)
{
	int nx = t.nx;
	int ny = t.ny;
	int nz = t.nz;
	return sqrt(nx*nx + ny*ny + nz*nz);
}

void InfiniteMatterSP::Connected()
{
	conState temp;

	int rp0, rp = 0, nshell, nf = 0;
	for(auto &t : fStates){
		nf++;
		rp0 = rp;
		rp = Radius(t);
		cout << rp0 << "\t" << rp << endl; getchar();
		if(rp - rp0 > 0) nshell++;
		if(nshell >= fNShell + 1) break;
	}

	cout << "nf: " << nf - 1 << endl;

	for (auto &t : fStates)
	//	for (int j = 0; j < size; j++)
		{
			//if(i == j) continue;

			// checking for hole states
			if (Radius(t) < fNShell - 1)// && Radius(fStates[j]) < fNShell)
			{
				temp.part2 = 1;
				fConStates.push_back(temp);
			}

			// checking for particle states
		//	else if (Radius(fStates[i]) <= lambda && Radius(fStates[j]) <= lambda)
			{

			}
		}
		cout << "The size of the hole states is: " << fConStates.size() << endl;
}



				// temp.Pijx = fStates[i].kx + fStates[j].kx;
				// temp.Pijy = fStates[i].ky + fStates[j].ky;
				// temp.Pijz = fStates[i].kz + fStates[j].kz;

				// temp.kRelx = (fStates[i].kx - fStates[j].kx) / 2;
				// temp.kRely = (fStates[i].ky - fStates[j].ky) / 2;
				// temp.kRelz = (fStates[i].kz - fStates[j].kz) / 2;

				// temp.two_Sij = fStates[i].two_s + fStates[j].two_s;

void InfiniteMatterSP::Print()
{	
	cout << "The quantum numbers are: " << endl;
	cout << "nx \t" << "ny \t" << "nz \t" <<  "two_s" << endl;
	for (auto t : fStates)
		cout << t.nx << " \t" << t.ny << " \t" << t.nz << " \t" << t.two_s << endl;
	cout << "The number of single particle states is: " << fStates.size() << endl;
	cout << "Done!" << endl;
}