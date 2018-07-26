#include <iostream>
#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

//////////////// Coupled Cluster ///////////////////

const int particle_no = 8;
const int n_mu = 8; // Available positions above Fermion
const double g = 0.5;
const double epsilon = 1.0;

/////// Fuction <rs|V|pq> ////////
double V_pqrs(int p, int q, int r, int s)
{
	///// If we use V_2body in form of Matrix, we should save pp-pp, hh-hh, ph-ph ... to save the space.
	if((p/2 == q/2) && (r/2 == s/2)) return (-g/2);
	else return 0.;
}

double f_pq(int p, int q)
{
	/////// f_pq = <q|v^(1)|p> + Sum_i(<iq|v^(2)|ip>)////////
	////// In pairing model, f_pq = sp_energy_p + V(pp,pp).///////
	double sum = 0.;
	double sp_energy = (p/2) * epsilon;
	if((p==q) && (p<particle_no))
	{
		for (int i = 0; i < particle_no; i++)
		{
			if((i!=p) && ((i/2)==(p/2)))
			{
				sum += sp_energy + V_pqrs(i,p,i,q);
			}
		}

	}
	else if((p==q) && (p >= particle_no))
	{
		sum = sp_energy;
	}
	else sum = 0.;
	return sum;
}


double Fill_HN()
{
	//////// apply for HN_bar(ij,ab) 4D-Matrix ////////
	//////// apply for t_(ij,ab) 4D-Matrix ////////
	double ****HN_bar = (double ****)malloc(particle_no * sizeof(double***));
	double ****t_ijab = (double ****)malloc(particle_no * sizeof(double***));
	for(int i=0; i<particle_no; i++)
	{
		HN_bar[i] = (double***)malloc(particle_no * sizeof(double**));
		t_ijab[i] = (double***)malloc(particle_no * sizeof(double**));
		for(int j = 0; j<particle_no; j++)
		{
			HN_bar[i][j] = (double**)malloc(n_mu * sizeof(double*));
			t_ijab[i][j] = (double**)malloc(n_mu * sizeof(double*));
			for(int k = 0; k<n_mu; k++)
			{
				HN_bar[i][j][k] = (double*)malloc(n_mu * sizeof(double));
				t_ijab[i][j][k] = (double*)malloc(n_mu * sizeof(double));
			}
		}
	}


	/////// Initialized the t_ijab. Trial t_ijab_0. and HN_bar matrix. /////////
	for(int i = 0; i<particle_no; i++) 
	{
		for(int j = 0; j<particle_no; j++)
		{
			for(int a = particle_no; a < (particle_no + n_mu); a++)
			{
				for(int b = particle_no; b < (particle_no + n_mu); b++)
				{
					t_ijab[i][j][a - particle_no][b - particle_no] = V_pqrs(i,j,a,b) / ( f_pq(a,a) + f_pq(b,b) - f_pq(i,i) - f_pq(j,j) ); //trial solution t_ijab_0.
//					HN_bar[i][j][a-particle_no][b-particle_no] = 0.;
				}
			}
		}
	}

	//////// FILL HAMILTONIAN HN /////////
	double Ec = 0.; // Correlation energy, Ec = Ecc - Eref.
	int i=0;
	double temp1, temp2, delta;
	do
	{
		temp1 = Ec;
		cout<<i<<" iteration, Ec_"<<i<<" = "<<temp1;
		for(int i = 0; i<particle_no; i++) 
		{
			for(int j = 0; j<particle_no; j++)
			{
				for(int a = particle_no; a < (particle_no + n_mu); a++)
				{
					for(int b = particle_no; b < (particle_no + n_mu); b++)
					{
						Ec = 0.25 * V_pqrs(a,b,i,j) * t_ijab[i][j][a-particle_no][b-particle_no];
						double h[10] = {0.};

						h[0] = V_pqrs(i,j,a,b);

						for(int k = 0; k<particle_no; k++)
						{
							h[1] += (-1) * ( f_pq(j,k) * t_ijab[i][k][a-particle_no][b-particle_no] - f_pq(i,k) * t_ijab[j][k][a-particle_no][b-particle_no]);
						}
						for(int c = particle_no; c < (particle_no + n_mu); c++)
						{
							h[2] += f_pq(c,b) * t_ijab[i][j][a-particle_no][c-particle_no] - f_pq(c,a) * t_ijab[i][j][b-particle_no][c-particle_no];
						}
						for(int k = 0; k<particle_no; k++)
						{
							for(int c = particle_no; c < (particle_no+n_mu); c++)
							{
								h[3] += (t_ijab[i][k][a-particle_no][c-particle_no] * V_pqrs(c,j,k,b) -t_ijab[i][k][b-particle_no][c-particle_no]*V_pqrs(c,j,k,a)) - (t_ijab[j][k][a-particle_no][c-particle_no]*V_pqrs(c,i,k,b) + t_ijab[j][k][b-particle_no][c-particle_no]*V_pqrs(c,i,k,a));
								for(int l = 0; l < particle_no; l++)
								{
									for(int d = particle_no; d<(particle_no+n_mu); d++)
									{
										h[6] += (-1) * 0.5 * (t_ijab[i][j][a-particle_no][c-particle_no]*V_pqrs(c,d,k,l)*t_ijab[k][l][b-particle_no][d-particle_no] - t_ijab[i][j][b-particle_no][c-particle_no]*V_pqrs(c,d,k,l)*t_ijab[k][l][a-particle_no][d-particle_no]);
										h[7] += (-1) * 0.5 * (t_ijab[i][k][a-particle_no][b-particle_no]*V_pqrs(c,d,k,l)*t_ijab[j][l][c-particle_no][d-particle_no] - t_ijab[j][k][a-particle_no][b-particle_no]*V_pqrs(c,d,k,l)*t_ijab[i][l][c-particle_no][d-particle_no]);
										h[8] += 0.5 * (t_ijab[i][k][a-particle_no][c-particle_no]*V_pqrs(c,d,k,l)*t_ijab[l][j][d-particle_no][b-particle_no] - t_ijab[i][k][b-particle_no][c-particle_no]*V_pqrs(c,d,k,l)*t_ijab[l][j][d-particle_no][a-particle_no] - t_ijab[j][k][a-particle_no][c-particle_no]*V_pqrs(c,d,k,l)*t_ijab[l][i][d-particle_no][b-particle_no] + t_ijab[j][k][b-particle_no][c-particle_no]*V_pqrs(c,d,k,l)*t_ijab[l][i][d-particle_no][a-particle_no]);
										h[9] += 0.25 * (t_ijab[i][j][a-particle_no][b-particle_no]*V_pqrs(c,d,k,l)*t_ijab[k][l][a-particle_no][b-particle_no]);
									}
								}
							}						
						}
						for(int k =0; k<particle_no; k++)
						{
							for (int l = 0; l < particle_no; l++)
							{
								h[4] += 0.5 * t_ijab[k][l][a-particle_no][b-particle_no] * V_pqrs(i,j,k,l);
							}
						}
						for(int c = particle_no; c<(particle_no+n_mu); c++)
						{
							for(int d = particle_no; d<(particle_no+n_mu); d++)
							{
								h[5] += 0.5 * t_ijab[i][j][c-particle_no][d-particle_no]*V_pqrs(c,d,a,b);
							}
						}
						HN_bar[i][j][a-particle_no][b-particle_no] = 0.;
						for(int m = 0; m<1; m++)
						{
							HN_bar[i][j][a - particle_no][b - particle_no] += h[m];
						}
					}
				}
			}
			
			for(int i = 0; i<particle_no; i++) 
			{
				for(int j = 0; j<particle_no; j++)
				{
					for(int a = particle_no; a < (particle_no + n_mu); a++)
					{
						for(int b = particle_no; b < (particle_no + n_mu); b++)
						{
						//	if(t_ijab[i][j][a-particle_no][b-particle_no]!=0)cout<<"attention!"<<t_ijab[i][j][a-particle_no][b-particle_no]<<endl;
							t_ijab[i][j][a-particle_no][b-particle_no] -= HN_bar[i][j][a-particle_no][b-particle_no] / (f_pq(a,a) + f_pq(b,b) - f_pq(i,i) - f_pq(j,j));
						}
					}
				}
			}

		}
		
		temp2 = Ec;
		delta = fabs(temp2 - temp1);
		cout<<", Ec_"<<i<<" = "<<temp2<<", d_Ec = "<<delta<<endl;
		
		i++;

		if(i>10) break;

	}while(delta>0.001);


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

	return Ec;
}


int main()
{
//	cout<<"hello"<<endl;
	Fill_HN();
//cout<<V_pqrs(0,1,10,11)<<endl;
	return 0;

}