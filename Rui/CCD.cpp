////////********** Coupled Cluster **********/////////

#include <iostream>
#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

const int particle_no = 4;
const int n_mu = 4; // Available positions above Fermion
//const double g = -0.5;
const double epsilon = 1.0;

/////// Fuction <rs|V|pq> ////////
///// If we use V_2body in form of Matrix, we should save pp-pp, hh-hh, ph-ph ... to save the space.
double V_pqrs(int p, int q, int r, int s, double g)
{
	if((p/2 == q/2) && (r/2 == s/2) && (p!=q) && (r!=s))
	{
		if((p<q)&&(r<s))return (-g/2);
		else if((p<q)&&(r>s))return (g/2);
		else if((p>q)&&(r<s))return (g/2);
		else return (-g/2); // means <12|12>, <12|21>,<21|12>,<21|21>
	}
	else return 0.;
}

/////// f_pq = <q|v^(1)|p> + Sum_i(<iq|v^(2)|ip>)////////
////// In pairing model, below Fermion -> f_pq = sp_energy_p + V(ip,ip)(only if i and p are in the same level)///////
//////////////////////// above Fermion -> f_pq = sp_energy_p /////////////////
double f_pq(int p, int q, double g)
{
	double sum = 0.;
	double sp_energy = (p/2) * epsilon;
	if((p==q) && (p < particle_no))
	{
		for (int i = 0; i < particle_no; i++)
		{
			if((i!=p) && ((i/2)==(p/2)))
			{
				sum += sp_energy + V_pqrs(i,p,i,q,g);
			}
		}
		return sum;
	}
	else if((p==q) && (p >= particle_no))
	{
		sum = sp_energy;
		return sum;
	}
	else 
	{
		sum = 0.;
		return sum;
	}
}


double Fill_HN(double g)
{
	//////// Allocate for HN_bar(ij,ab) and t_(ij,ab) 4D-Matrix ////////
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

	//Ec is the correlation energy, Ec = Ecc - Eref.
	double Ec = 0.;
//	printf("Ec_0 = -%.7f\n",Ec);


	/////// Initialized the t_ijab. Trial solution t_ijab_0 /////////
	for(int i = 0; i<particle_no; i++) 
	{
		for(int j = 0; j<particle_no; j++)
		{
			for(int a = particle_no; a < (particle_no + n_mu); a++)
			{
				for(int b = particle_no; b < (particle_no + n_mu); b++)
				{
					t_ijab[i][j][a - particle_no][b - particle_no] = V_pqrs(i,j,a,b,g) / ( -f_pq(a,a,g) - f_pq(b,b,g) + f_pq(i,i,g) + f_pq(j,j,g) );
				}
			}
		}
	}

	//////// FILL HAMILTONIAN HN /////////
	int i=0;
	double h[10];
	double temp1, temp2, delta;
	do
	{
		temp1 = Ec;
		Ec = 0.; 
		for(int i = 0; i<particle_no; i++) 
		{
			for(int j = 0; j<particle_no; j++)
			{
				for(int a = particle_no; a < (particle_no + n_mu); a++)
				{
					for(int b = particle_no; b < (particle_no + n_mu); b++)
					{
						Ec += 0.25 * V_pqrs(a,b,i,j,g) * t_ijab[i][j][a-particle_no][b-particle_no];
					
						for(int loop1 = 0; loop1 <10; loop1++)
						{
							h[loop1] = 0.;	
						}

						h[0] = V_pqrs(i,j,a,b,g);

					//	h[1] = (a/2+b/2-i/2-j/2) * epsilon * t_ijab[i][j][a-particle_no][b-particle_no];
							
						for(int k = 0; k< particle_no; k++)
						{
							h[1] += (-1) * ( f_pq(j,k,g) * t_ijab[i][k][a-particle_no][b-particle_no] - f_pq(i,k,g) * t_ijab[j][k][a-particle_no][b-particle_no]);
						}

						for(int c = particle_no; c < (particle_no + n_mu); c++)
						{
							h[2] += f_pq(c,b,g) * t_ijab[i][j][a-particle_no][c-particle_no] - f_pq(c,a,g) * t_ijab[i][j][b-particle_no][c-particle_no];
						}
						/*for(int c = particle_no; c<(particle_no+n_mu); c++)
						{
							for(int d = particle_no; d<(particle_no+n_mu); d++)
							{
								h[2] += 0.5 * V_pqrs(c,d,a,b,g) * t_ijab[i][j][c-particle_no][d-particle_no];
							}
						}*/
					/*	for(int k = 0; k<particle_no; k++)
						{
							for(int l=0; l<particle_no; l++)
							{
								h[3] += 0.5 * V_pqrs(i,j,k,l,g) * t_ijab[k][l][a-particle_no][b-particle_no];
							}
						}
						for(int k = 0; k<particle_no; k++)
						{
							for(int c = particle_no; c<(particle_no+n_mu); c++)
							{
								h[4] += V_pqrs(c,j,k,b,g) * t_ijab[i][k][a-particle_no][c-particle_no] - V_pqrs(c,j,k,a,g) * t_ijab[i][k][b-particle_no][c-particle_no] - V_pqrs(c,i,k,b,g) * t_ijab[j][k][a-particle_no][c-particle_no] + V_pqrs(c,i,k,a,g) * t_ijab[j][k][b-particle_no][c-particle_no];
							}
						}

						for(int k = 0; k<particle_no; k++)
						{
							for(int c = particle_no; c < (particle_no+n_mu); c++)
							{
								for(int l = 0; l < particle_no; l++)
								{
									for(int d = particle_no; d<(particle_no+n_mu); d++)
									{
										h[5] += 0.25 * V_pqrs(c,d,k,l,g) * t_ijab[i][j][c-particle_no][d-particle_no] * t_ijab[k][l][a-particle_no][b-particle_no];
										h[6] += V_pqrs(c,d,k,l,g) * t_ijab[i][k][a-particle_no][c-particle_no] * t_ijab[j][l][b-particle_no][d-particle_no] - V_pqrs(c,d,k,l,g) * t_ijab[j][k][a-particle_no][c-particle_no] * t_ijab[i][l][b-particle_no][d-particle_no];
										h[7] += (-0.5) * V_pqrs(c,d,k,l,g) * t_ijab[i][k][d-particle_no][c-particle_no] * t_ijab[l][j][a-particle_no][b-particle_no] - (-0.5) * V_pqrs(c,d,k,l,g) * t_ijab[j][k][d-particle_no][c-particle_no] * t_ijab[l][i][a-particle_no][b-particle_no];
										h[8] += (-0.5) * V_pqrs(c,d,k,l,g) * t_ijab[l][k][a-particle_no][c-particle_no] * t_ijab[i][j][d-particle_no][b-particle_no] - (-0.5) * V_pqrs(c,d,k,l,g) * t_ijab[l][k][b-particle_no][c-particle_no] * t_ijab[i][j][d-particle_no][a-particle_no];
									}
								}
							}
						}*/
						
						for(int k = 0; k<particle_no; k++)
						{
							for(int c = particle_no; c < (particle_no+n_mu); c++)
							{
								h[3] += (t_ijab[i][k][a-particle_no][c-particle_no] * V_pqrs(c,j,k,b,g) -t_ijab[i][k][b-particle_no][c-particle_no]*V_pqrs(c,j,k,a,g)) - (t_ijab[j][k][a-particle_no][c-particle_no]*V_pqrs(c,i,k,b,g) + t_ijab[j][k][b-particle_no][c-particle_no]*V_pqrs(c,i,k,a,g));
								for(int l = 0; l < particle_no; l++)
								{
									for(int d = particle_no; d<(particle_no+n_mu); d++)
									{
										h[6] += (-1) * 0.5 * (t_ijab[i][j][a-particle_no][c-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[k][l][b-particle_no][d-particle_no] - t_ijab[i][j][b-particle_no][c-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[k][l][a-particle_no][d-particle_no]);
										h[7] += (-1) * 0.5 * (t_ijab[i][k][a-particle_no][b-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[j][l][c-particle_no][d-particle_no] - t_ijab[j][k][a-particle_no][b-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[i][l][c-particle_no][d-particle_no]);
										h[8] += 0.5 * (t_ijab[i][k][a-particle_no][c-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[l][j][d-particle_no][b-particle_no] - t_ijab[i][k][b-particle_no][c-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[l][j][d-particle_no][a-particle_no] - t_ijab[j][k][a-particle_no][c-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[l][i][d-particle_no][b-particle_no] + t_ijab[j][k][b-particle_no][c-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[l][i][d-particle_no][a-particle_no]);
										h[9] += 0.25 * (t_ijab[i][j][c-particle_no][d-particle_no]*V_pqrs(c,d,k,l,g)*t_ijab[k][l][a-particle_no][b-particle_no]);
									}
								}
							}						
						}

						for(int k =0; k<particle_no; k++)
						{
							for (int l = 0; l < particle_no; l++)
							{
								h[4] += 0.5 * t_ijab[k][l][a-particle_no][b-particle_no] * V_pqrs(i,j,k,l,g);
							}
						}

						for(int c = particle_no; c<(particle_no+n_mu); c++)
						{
							for(int d = particle_no; d<(particle_no+n_mu); d++)
							{
								h[5] += 0.5 * t_ijab[i][j][c-particle_no][d-particle_no]*V_pqrs(c,d,a,b,g);
							}
						}

						HN_bar[i][j][a-particle_no][b-particle_no] = 0.;
						for(int m = 0; m<10; m++)
						{
							HN_bar[i][j][a - particle_no][b - particle_no] += h[m];
						}
					}
				}
			}
		}

		double t_no_this, t_last;
		double alpha = 0.5; //*** Here we could try different alpha!!! ****//
		
		///// Renew t_ijab matrix elements /////
		for(int i = 0; i<particle_no; i++) 
		{
			for(int j = 0; j<particle_no; j++)
			{
				for(int a = particle_no; a < (particle_no + n_mu); a++)
				{
					for(int b = particle_no; b < (particle_no + n_mu); b++)
					{
				//		t_ijab[i][j][a-particle_no][b-particle_no] -= HN_bar[i][j][a-particle_no][b-particle_no] / (f_pq(a,a,g) + f_pq(b,b,g) - f_pq(i,i,g) - f_pq(j,j,g));
					//////**** HERE we should use MIXING !!!****///////
						t_last = t_ijab[i][j][a-particle_no][b-particle_no];
						t_no_this = t_ijab[i][j][a-particle_no][b-particle_no] - HN_bar[i][j][a-particle_no][b-particle_no] / (f_pq(a,a,g) + f_pq(b,b,g) - f_pq(i,i,g) - f_pq(j,j,g));
						t_ijab[i][j][a-particle_no][b-particle_no] = alpha * t_no_this + (1-alpha) * t_last;
					}
				}
			}
		}
		
		
		temp2 = Ec;
		delta = fabs(temp1 - temp2);
		cout<<"Ec_"<<i<<" = "<<temp2<<", d_Ec = "<<delta<<endl;
		
		i++;

		if(i>1000) break;

	}while(delta>10e-9);


	////////// Free HN and t_ijab/////////
	for(int i = 0; i<particle_no; i++)
	{
		for(int j =0; j<particle_no; j++)
		{
			for(int k = 0; k<n_mu; k++)
			{
				free(HN_bar[i][j][k]);
				free(t_ijab[i][j][k]);
	
			}
		}
	}
	for(int i = 0; i<particle_no; i++)
	{
		for(int j = 0; j<particle_no; j++)
		{
			free(HN_bar[i][j]);
			free(t_ijab[i][j]);
		}
	}
	for (int i = 0; i < particle_no; i++)
	{
		free(HN_bar[i]);
		free(t_ijab[i]);
	}
	free(HN_bar);
	free(t_ijab);

	return Ec;
}


int main()
{
	double g;
//	for(int i =-19; i<20; i++)
//	{	
//		g = ((double)i)/20.;
//		cout<<g<<" "<<Fill_HN(g)<<endl;
//	}
	Fill_HN(-0.5);
	return 0;

}