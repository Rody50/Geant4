#include <iostream>
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

using namespace std;

//////////////// Coupled Cluster ///////////////////

const int particle_no = 10;
const int n_mu = 100; // Available positions above Fermion
const double g = 0.5;

/////// Fuction <rs|V|pq> ////////
double V_pqrs(int p, int q, int r, int s)
{
	///// If we use V_2body in form of Matrix, we should save pp-pp, hh-hh, ph-ph ... to save the space.
	if(p/2 == q/2 && r/2 == s/2) return (-g/2);
	else return 0.;
}

double f_pq(int p, int q)
{
	/////// f_pq = <q|v^(1)|p> + Sum_i(<iq|v^(2)|ip>)////////
	////// In pairing model, the 1st term <q||p> would be 0.///////
	double sum = 0.;
	for(int i = 0; i < particle_no; i++)
	{
		sum += V_pqrs(i,p,i,q);
	}
	return sum;
}

void Fill_HN()
{
	//////// HN_bar(ij,ab) 4D-Matrix ////////
	double ****HN_bar = (double ****)malloc(particle_no * sizeof(double***));
	for(int i=0; i<particle_no; i++)
	{
		HN_bar[i] = (double***)malloc(particle_no * sizeof(double**));
		for(int j = 0; j<particle_no; j++)
		{
			HN_bar[i][j] = (double**)malloc(n_mu * sizeof(double*));
			for(int k = 0; k<particle_no; k++)
			{
				HN_bar[i][j][k] = (double*)malloc(particle_no * sizeof(double));
			}
		}
	}

	//////// t_(ij,ab) 4D-Matrix ////////
	double ****t_ijab = (double ****)malloc(particle_no * sizeof(double***));
	for(int i=0; i<particle_no; i++)
	{
		t_ijab[i] = (double***)malloc(particle_no * sizeof(double**));
		for(int j = 0; j<particle_no; j++)
		{
			t_ijab[i][j] = (double**)malloc(n_mu * sizeof(double*));
			for(int k = 0; k<particle_no; k++)
			{
				t_ijab[i][j][k] = (double*)malloc(particle_no * sizeof(double));
			}
		}
	}


//////////// Initialized the t_ijab. Trial t_ijab_0. /////////
	for(int i = 0; i<particle_no; i++) 
	{
		for(int j = i+1; j<particle_no; j++)
		{
			for(int a = particle_no; a < (particle_no + n_mu) ; a++)
			{
				for(int b = a+1; b < (particle_no + n_mu); b++)
				{
					t_ijab[i][j][a][b] = V_pqrs(i,j,a,b) / ( f_pq(a,a) + f_pq(b,b) - f_pq(i,i) - f_pq(j,j) );
				}
			}
		}
	}




	////////// Free HN /////////
	for(int i = 0; i<particle_no; i++)
	{
		for(int j =0; j<particle_no; j++)
		{
			for(int k = 0; k<n_mu; k++)
			{
				free(HN_bar[i][j][k]);
			}
		}
	}
	for(int i = 0; i<particle_no; i++)
	{
		for(int j = 0; j<particle_no; j++)
		{
			free(HN_bar[i][j]);
		}
	}
	for (int i = 0; i < particle_no; i++)
	{
		free(HN_bar[i]);
	}
	free(HN_bar);

	////////// Free t_ijab /////////
	for(int i = 0; i<particle_no; i++)
	{
		for(int j =0; j<particle_no; j++)
		{
			for(int k = 0; k<n_mu; k++)
			{
				free(t_ijab[i][j][k]);
			}
		}
	}
	for(int i = 0; i<particle_no; i++)
	{
		for(int j = 0; j<particle_no; j++)
		{
			free(t_ijab[i][j]);
		}
	}
	for (int i = 0; i < particle_no; i++)
	{
		free(t_ijab[i]);
	}
	free(t_ijab);
}


int main()
{
	cout<<"hello"<<endl;


	return 0;

}