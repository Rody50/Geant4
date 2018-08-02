#include <iostream>
#include <iostream>
#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <time.h>

using namespace std;

#define pi 3.1415926

const int magic_no = 66;
const int Nmax = 4;
const double mass = 938.92; // in unit of MeV
const double h_bar = 197.33;


//********************************************************************//
//****************************  STRUCTURE  ***************************//
//********************************************************************//
struct spstate
{
	int kx; //Nmax = 4, so kx, ky ,kz have (2*Nmax+1).
	int ky;
	int kz;
	int sz; //spin_up = +1, spin_down = -1;
	int occupied; //occupied = 1, unoccupied = 0;
};

struct tbwf//two-body wave function in relative momentum presentation
{
	int Q_kx;
	int Q_ky;
	int Q_kz;
	double q_kx;
	double q_ky;
	double q_kz;
	int sz1;
	int sz2;
	int flag;
};

struct CHANNEL
{
	int Qx; //record the same total momentum
	int Qy;
	int Qz;
	int times; //record how many times this Q occurs.
	tbwf *wf;
};

struct MATRIX_CHANNEL
{
	int Qx;
	int Qy;
	int Qz;
	double **matrix; 
};


//********************************************************************//
//****************************  FUNCTIONS  ***************************//
//********************************************************************//
double factorial(double n) // calculate the factorial
{
	double result = n;
	do
	{
		result *= (n-1);
		n--;
	} while(n != 1);

	return result;
}

double combinator(double n, double m) //calculate the combinator
{
	double result;
	result = factorial(n) / (factorial(m) * factorial(n-m));
	return result;
}

double kronecker(int s1, int s2)
{
	if(s1 == s2) return 1.;
	else return 0;
}

///*** DEFINE THE V FUNCTION ***///
double V_ks(double L, double q_rs_x, double q_rs_y, double q_rs_z, int sz_r, int sz_s, double q_pq_x, double q_pq_y, double q_pq_z, int sz_p, int sz_q)
{
	double result=0;
	double qx, qy, qz, qq;
	qx = (q_pq_x - q_rs_x);
	qy = (q_pq_y - q_rs_y);
	qz = (q_pq_z - q_rs_z);
	qq = (pow(qx,2) + pow(qy,2) + pow(qz,2)) * pow((2*pi/L),2);

	result = 200./pow(L,3) * pow((pi/1.487), 1.5) * exp(-qq/(4*1.487));

	result += -91.85/pow(L,3) * pow((pi/0.465), 1.5) * exp(-qq/(4*0.465));

	result *= kronecker(sz_p, sz_r) * kronecker(sz_q, sz_s) - kronecker(sz_p, sz_s) * kronecker(sz_q, sz_r);
	return result*0.5;
}

///*** DEFINE THE V FUNCTION FOR ANTI-SYMMETRY***///
double V_ks_AS(double L, double q_rs_x, double q_rs_y, double q_rs_z, int sz_r, int sz_s, double q_pq_x, double q_pq_y, double q_pq_z, int sz_p, int sz_q)
{
	double result;
	result = V_ks(L, q_rs_x, q_rs_y, q_rs_z, sz_r, sz_s, q_pq_x, q_pq_y, q_pq_z, sz_p, sz_q) \
			 - V_ks(L, -q_rs_x, -q_rs_y, -q_rs_z, sz_s, sz_r, q_pq_x, q_pq_y, q_pq_z, sz_p, sz_q);
	return result;
}


///*** DEFINE THE f FUNCTION***///
double f_ks(struct spstate *config, int spstates_no, double L, int kp_x, int kp_y, int kp_z, int sp, int kq_x, int kq_y, int kq_z, int sq)
{
	double result = 0;
	double kk;
	double ip_kx, ip_ky, ip_kz;
	if(kp_x == kq_x && kp_y == kq_y && kp_z == kq_z && sp == sq)
	{
		kk = (pow(kp_x,2) + pow(kp_y,2) + pow(kp_z,2)) * pow(2*pi/L,2);
		result = pow(h_bar,2) * kk / (2. * mass);
		for(int i = 0; i<spstates_no; i++)
		{
			if(config[i].occupied == 1)
			{
				if(config[i].kx == kp_x && config[i].ky == kp_y && config[i].kz == kp_z && config[i].sz == sp) continue;
				ip_kx = (double)(config[i].kx - kp_x)/2.;
				ip_ky = (double)(config[i].ky - kp_y)/2.;
				ip_kz = (double)(config[i].kz - kp_z)/2.;

				result += V_ks_AS(L, ip_kx, ip_ky, ip_kz, config[i].sz, sp,ip_kx, ip_ky, ip_kz, config[i].sz, sp);
			}
		}
	}
	return result;
}




void CCD_neutronmatter(double rho) //generate single particle state
{
	double L;
	L = pow((double)magic_no/rho, 1./3.);

	int spstates_no = 2 * pow((2*Nmax+1),3);
	spstate *config = (spstate *)malloc(spstates_no * sizeof(spstate));


	////////////////////////////////////////////////////////////////////////
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	//*******************  GENERATE ALL THE S.P. STATE  ******************//
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////
	//*****************************************************************************//
	//*** Actually it is a MAP more than fullfill the s.p. configuration space. ***//
	//*****************************************************************************//
	/////////////////////////////////////////////////////////////////////////////////

	/*	for(int i = 0; i<spstates_no; i++)
	{
		config[i].sz = (i%2) - 0.5 ;
		config[i].kz = (i/2) % (2*(2*Nmax+1)) - Nmax;
		config[i].ky = ((i/2) / (2*(2*Nmax+1))) % (2*(2*Nmax+1)) - Nmax;
		config[i].kx = (((i/2) / (2*(2*Nmax+1))) / (2*(2*Nmax+1))) % (2*(2*Nmax+1)) - Nmax;
	}*/



	////////////////////////////////////////////////////////////////////////
	//********************************************************************//
	//*** Alternative choice to generate single-particle configuration ***//
	//*******  More straightforward but might be less efficient.  ********//
	//********************************************************************//
	////////////////////////////////////////////////////////////////////////

	int loop1 = 0;
	double mod;
	do
	{
		for(int kx = -Nmax; kx <= Nmax; kx++)
		{
			for(int ky = -Nmax; ky <= Nmax; ky++)
			{
				for (int kz = -Nmax; kz <= Nmax; kz++)
				{
					for(int sz = -1 ; sz <=1; sz+=2)
					{
						config[loop1].kx = kx;
						config[loop1].ky = ky;
						config[loop1].kz = kz;
						config[loop1].sz = sz;
						mod = pow(pow(config[loop1].kx,2) + pow(config[loop1].ky,2)+pow(config[loop1].kz,2),0.5);
						if(mod == 0 || mod == 1) config[loop1].occupied = 1;
						else config[loop1].occupied = 0;
						loop1++;
					}
				}
			}
		}
	}while(loop1<spstates_no);


	/////////////////////////////////////////////////////////////////////////
	//*********************************************************************//
	//*** GENERATE THE CONFIGURATIONS IN RELATIVE MOMENTUM PRESENTATION ***//
	//********************************************************************///
	/////////////////////////////////////////////////////////////////////////

//	double pp_no = (spstates_no - magic_no) * (spstates_no - magic_no - 1);//combinator((spstates_no - magic_no) , 2); // number of pp-configurations in relative presentation
	double hh_no_old = 2 * combinator(magic_no, 2);
//	double ph_no = 2 * (spstates_no - magic_no) * magic_no;

//	tbwf *ph_config = (tbwf *)malloc(ph_no * sizeof(tbwf));



	int loop2 = 0;
	int loop3 = 0;
	int loop4 = 0;


	////////  pp-configuration-number /////
        double pp_no = 0;
	for(int k1 = 0; k1<spstates_no; k1++)
	{
		for(int k2 = 0; k2<spstates_no; k2++)
		{
			if((k1 != k2) && (config[k1].occupied == 0) && (config[k2].occupied == 0) &&(config[k1].sz*config[k2].sz == -1))
			{
				pp_no ++;
			}
		}
	}
//	cout <<  "pp_no_new = "<<pp_no<<endl;
        ////////  pp-configuration-setup /////
	tbwf *pp_config = (tbwf *)malloc(pp_no * sizeof(tbwf)); //configurations in relative presentation
	
        	for(int k1 = 0; k1<spstates_no; k1++)
		{
			for(int k2 = 0; k2<spstates_no; k2++)
			{
				if((k1 != k2) && (config[k1].occupied == 0) && (config[k2].occupied == 0) && (config[k1].sz*config[k2].sz == -1)) // S12 =0 for this V
				{
				//	cout<<"!!"<<endl;
					pp_config[loop2].Q_kx = config[k1].kx + config[k2].kx;
					pp_config[loop2].Q_ky = config[k1].ky + config[k2].ky;
					pp_config[loop2].Q_kz = config[k1].kz + config[k2].kz;
					
					pp_config[loop2].q_kx = (double)(config[k1].kx - config[k2].kx) / 2.;
					pp_config[loop2].q_ky = (double)(config[k1].ky - config[k2].ky) / 2.;
					pp_config[loop2].q_kz = (double)(config[k1].kz - config[k2].kz) / 2.;
				
					pp_config[loop2].sz1 = config[k1].sz;
					pp_config[loop2].sz2 = config[k2].sz;

					pp_config[loop2].flag = 1;
					
					loop2++;
				}
			}
		}

//	        cout <<  "loop2 = "<<loop2<<endl;
	


	////////  hh-configuration-number /////
        double hh_no = 0;
	for(int k1 = 0; k1<spstates_no; k1++)
	{
		for(int k2 = 0; k2<spstates_no; k2++)
		{
			if((k1 != k2) && (config[k1].occupied == 1) && (config[k2].occupied == 1) &&  (config[k1].sz*config[k2].sz == -1))
			{
				hh_no ++;
			}
		}
	}

	////////  pp-configuration-setup /////
 	tbwf *hh_config = (tbwf *)malloc(hh_no * sizeof(tbwf));
		for(int k1 = 0; k1<spstates_no; k1++)
		{
			for(int k2 = 0; k2<spstates_no; k2++)
			{
				if((k1 != k2) && (config[k1].occupied == 1) && (config[k2].occupied == 1) &&  (config[k1].sz*config[k2].sz == -1))
				{
					hh_config[loop3].Q_kx = config[k1].kx + config[k2].kx;
					hh_config[loop3].Q_ky = config[k1].ky + config[k2].ky;
					hh_config[loop3].Q_kz = config[k1].kz + config[k2].kz;
					
					hh_config[loop3].q_kx = (double)(config[k1].kx - config[k2].kx) / 2.;
					hh_config[loop3].q_ky = (double)(config[k1].ky - config[k2].ky) / 2.;
					hh_config[loop3].q_kz = (double)(config[k1].kz - config[k2].kz) / 2.;
				
					hh_config[loop3].sz1 = config[k1].sz;
					hh_config[loop3].sz2 = config[k2].sz;
					hh_config[loop3].flag = 1;
					loop3++;
				}
			}
		}





	/////////////////////////////////////////////////////////////////////////
	//*********************************************************************//
	//***************************** BENCHMARK *****************************//
	//************** For Kinetic energy and HATREE-FOCK energy ************//
	//********************************************************************///
	/////////////////////////////////////////////////////////////////////////

		double Ek = 0.;
		double kk = 0.;
	//	cout<<"L = "<<L<<endl;
		for(int i = 0; i<spstates_no; i++)
		{
			if(config[i].occupied == 1)
			{
	//			cout<<config[i].kx<<" "<<config[i].ky<<" "<<config[i].kz<<endl;
				kk = (double)(pow(config[i].kx,2) + pow(config[i].ky,2) + pow(config[i].kz,2)) * pow((2*pi)/L,2);
	//			cout<<"kk = "<<kk<<endl;
				Ek += pow(h_bar,2) * kk / (2. * mass);
			}
		}
//		cout<<"h^2/2m="<<pow(h_bar,2)/2./mass;
//		cout<<"kinetic energy is "<<Ek/14.<<endl;
		double hf_energy;
		hf_energy = Ek;
		for(int ij = 0; ij < hh_no; ij++)
		{
			hf_energy += 0.5 * V_ks_AS(L, hh_config[ij].q_kx, hh_config[ij].q_ky, hh_config[ij].q_kz, hh_config[ij].sz1, hh_config[ij].sz2, \
			hh_config[ij].q_kx, hh_config[ij].q_ky, hh_config[ij].q_kz, hh_config[ij].sz1, hh_config[ij].sz2);
		}
		cout<<"density = "<<setw(5)<<rho<<" , (Ehf/N) = "  <<setw(7)<<hf_energy/magic_no<<"  ";

/*
		for(int k1 = 0; k1<spstates_no; k1++)
		{
			for(int k2 = 0; k2<spstates_no; k2++)
			{
				if((k1 != k2) && ((config[k1].occupied == 0 && config[k2].occupied == 1) || (config[k1].occupied == 1 && config[k2].occupied == 0)) )
				{
					ph_config[loop4].Q_kx = config[k1].kx + config[k2].kx;
					ph_config[loop4].Q_ky = config[k1].ky + config[k2].ky;
					ph_config[loop4].Q_kz = config[k1].kz + config[k2].kz;
					
					ph_config[loop4].q_kx = (double)(config[k1].kx - config[k2].kx) / 2.;
					ph_config[loop4].q_ky = (double)(config[k1].ky - config[k2].ky) / 2.;
					ph_config[loop4].q_kz = (double)(config[k1].kz - config[k2].kz) / 2.;
				
					ph_config[loop4].sz1 = config[k1].sz;
					ph_config[loop4].sz2 = config[k2].sz;
					ph_config[loop4].flag = 1;
					
					loop4++;
				}
			}
		}
*/

	////////////////////////////////////////////////////////////////////////
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	//***************  GENERATE THE CHANNEL CONFIGURATIONS  **************//
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	//*************************** HH-HH - CHANNEL  ***********************//
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	////////////////////////////////////////////////////////////////////////

	int k,loop5;
	int hh_temp_x, hh_temp_y, hh_temp_z;
	int hh_dimension = 0;

	double **temp_pointer1_V;
    double *temp_pointer2_V;

    double **temp_pointer1_t;
    double *temp_pointer2_t;

    double **temp_pointer1_H;
    double *temp_pointer2_H;

	for(int i = 0; i < hh_no; i++)
	{
		if(hh_config[i].flag == 1)
		{
			hh_temp_x = hh_config[i].Q_kx;
			hh_temp_y = hh_config[i].Q_ky;
			hh_temp_z = hh_config[i].Q_kz;
			for(int j = i; j<hh_no; j++)
			{
				if(hh_config[j].Q_kx == hh_temp_x && hh_config[j].Q_ky == hh_temp_y && hh_config[j].Q_kz == hh_temp_z)
				{
					hh_config[j].flag = 0;
				}
			}
			hh_dimension++; // 这记录了hh里有多少个不同的总动量
		}	
	}

	for(int i = 0; i<hh_no; i++)
	{
		hh_config[i].flag = 1;
	}



	CHANNEL hh_channel[hh_dimension]; //先声明！！！！！！！

	int hh_temp_dimension = 0;
	for(int i = 0; i< hh_no; i++)
	{
		if(hh_config[i].flag == 1)
		{
			hh_channel[hh_temp_dimension].Qx = hh_config[i].Q_kx;
			hh_channel[hh_temp_dimension].Qy = hh_config[i].Q_ky;
			hh_channel[hh_temp_dimension].Qz = hh_config[i].Q_kz;
			hh_channel[hh_temp_dimension].times = 0;

			for(int j = i; j<hh_no; j++)
			{
				if(hh_config[j].Q_kx == hh_channel[hh_temp_dimension].Qx && hh_config[j].Q_ky == hh_channel[hh_temp_dimension].Qy && hh_config[j].Q_kz == hh_channel[hh_temp_dimension].Qz)
				{
					hh_config[j].flag = 0;
					hh_channel[hh_temp_dimension].times++;
				}
			}

			//为每一个hh_channel[hh_temp_dimension]动态分配
			tbwf *temp_star_hh;
			temp_star_hh = (tbwf *)malloc(hh_channel[hh_temp_dimension].times * sizeof(tbwf));
			hh_channel[hh_temp_dimension].wf = temp_star_hh;
			
	//		hh_channel[hh_temp_dimension].wf = (tbwf *)malloc(sizeof(tbwf) * (p));

			hh_temp_dimension++;
		}
	}

	for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
	{
		hh_config[i].flag = 1;
	}


	loop5 = 0;
	k = 0;
	for(int i = 0; i< hh_no; i++)
	{
		if(hh_config[i].flag == 1)
		{
			k = 0;
			for(int j = i; j<hh_no; j++)
			{
				if(hh_config[i].Q_kx == hh_config[j].Q_kx && hh_config[i].Q_ky == hh_config[j].Q_ky && hh_config[i].Q_kz == hh_config[j].Q_kz)
				{
					hh_channel[loop5].wf[k].q_kx = hh_config[j].q_kx;
					hh_channel[loop5].wf[k].q_ky = hh_config[j].q_ky;
					hh_channel[loop5].wf[k].q_kz = hh_config[j].q_kz;
					hh_channel[loop5].wf[k].sz1 = hh_config[j].sz1;
					hh_channel[loop5].wf[k].sz2 = hh_config[j].sz2;
					hh_config[j].flag = 0;                  
					k++;
				}
			}
//			if(k != pp_channel[loop5].times) 
//			{
//					cout<<"k = "<<k<<" , times = "<<pp_channel[loop5].times<<endl;
//			}
			loop5++;
		}
	}


	//// Assignment for the V_ijkl ////
	//为每一个pp_channel[pp_temp_dimension]动态分配
	MATRIX_CHANNEL V_ijkl[hh_dimension];


	for(int i = 0; i< hh_dimension; i++)
	{
		temp_pointer1_V = (double **)malloc(hh_channel[i].times * sizeof(double*));
		V_ijkl[i].matrix = temp_pointer1_V;

		for(int cd = 0; cd < hh_channel[i].times; cd++)
		{
			temp_pointer2_V = (double *)malloc(hh_channel[i].times * sizeof(double));
			V_ijkl[i].matrix[cd] = temp_pointer2_V;
		}
	}
	
	///////////////////////////////////////////////////////////
	//*******************************************************//
	//************** FULFILL THE V_ij_kl MATRIX *************//
	//*******************************************************//
	//*******************************************************//

	for(int i = 0; i< hh_dimension; i++)
	{
//		cout<<"Q = "<<pp_channel[i].Qx<<" "<<pp_channel[i].Qy<<" "<<pp_channel[i].Qz;
//		cout<<" times = "<<pp_channel[i].times<<endl;
		V_ijkl[i].Qx = hh_channel[i].Qx;
		V_ijkl[i].Qy = hh_channel[i].Qy;
		V_ijkl[i].Qz = hh_channel[i].Qz;

		for(int ab = 0; ab < hh_channel[i].times; ab++)
		{
			for(int cd = 0; cd < hh_channel[i].times; cd++)
			{
				V_ijkl[i].matrix[cd][ab]= V_ks_AS(L, hh_channel[i].wf[cd].q_kx, hh_channel[i].wf[cd].q_ky, hh_channel[i].wf[cd].q_kz, hh_channel[i].wf[cd].sz1, hh_channel[i].wf[cd].sz2,\
				hh_channel[i].wf[ab].q_kx,hh_channel[i].wf[ab].q_ky, hh_channel[i].wf[ab].q_kz, hh_channel[i].wf[ab].sz1, hh_channel[i].wf[ab].sz2);
//				cout<<V_ijkl[i].matrix[cd][ab]<<" ";
			}
//			cout<<"\n";
		}
//		cout<<endl;
//		cout<<"///////////////"<<endl;
//		cout<<endl;
	}




	////////////////////////////////////////////////////////////////////////
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	//*************************  HH-PP - CHANNEL  ************************//
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	////////////////////////////////////////////////////////////////////////

	int hhpp_temp_x, hhpp_temp_y, hhpp_temp_z;
	int hhpp_dimension = 0;
	for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
	{
		hh_config[i].flag = 1;
	}
	for(int i = 0; i < hh_no; i++)
	{
		if(hh_config[i].flag == 1)
		{
			hhpp_temp_x = hh_config[i].Q_kx;
			hhpp_temp_y = hh_config[i].Q_ky;
			hhpp_temp_z = hh_config[i].Q_kz;
			for(int j = 0; j<hh_no; j++)
			{
				if(hh_config[j].Q_kx == hhpp_temp_x && hh_config[j].Q_ky == hhpp_temp_y && hh_config[j].Q_kz == hhpp_temp_z)
				{
					hh_config[j].flag = 0;
				}
			}
			for(int j = 0; j<pp_no; j++)
			{
				if(pp_config[j].Q_kx == hhpp_temp_x && pp_config[j].Q_ky == hhpp_temp_y && pp_config[j].Q_kz == hhpp_temp_z)
				{
					//hh_config[j].flag = 0;
					hhpp_dimension++;
					break;
				}
			}
			//hhpp_dimension++; // 这记录了hh里有多少个不同的总动量
		}	
	}
//	cout<<"hhpp_dimension = "<<hhpp_dimension<<endl;
//	for(int i = 0; i<hh_no; i++)
//	{
//		hh_config[i].flag = 1;
//	}



	CHANNEL hhpp_channel_L[hhpp_dimension]; //先声明！！！！！！！
	CHANNEL hhpp_channel_R[hhpp_dimension]; //先声明！！！！！！！
//	CHANNEL *pp_channel = malloc(sizeof(CHANNEL) * pp_dimension);
//	cout<< "pp_dimension ="<<pp_dimension<<endl;

//	cout<<"pp_dimension = "<<pp_dimension<<endl;



	for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
	{
		hh_config[i].flag = 1;
	}
	for(int i = 0; i<pp_no; i++) //将hh_config[i]的flag全都赋值回1
	{
		pp_config[i].flag = 1;
	}



	int hhpp_temp_dimension = 0;
    int Qx,Qy,Qz,times_L,times_R;


	for(int i = 0; i< hh_no; i++)
	{
		if(hh_config[i].flag == 1)
		{
            Qx = hh_config[i].Q_kx; //这里先不记录hhpp_channel_L[hhpp_temp_dimension].Q，因为这里事实上并未确定在pp里一定能找到匹配的hh
			Qy = hh_config[i].Q_ky;
			Qz = hh_config[i].Q_kz;
			times_L = 0;
			times_R = 0;

			for(int j = i; j<hh_no; j++)
			{
				if(hh_config[j].Q_kx == Qx && hh_config[j].Q_ky == Qy && hh_config[j].Q_kz == Qz)
				{
					hh_config[j].flag = 0;
					times_L ++;
					//hhpp_channel_L[hhpp_temp_dimension].times++;
				}
			}
		    for(int j = 0; j<pp_no; j++)
			{
				if(pp_config[j].Q_kx == Qx && pp_config[j].Q_ky == Qy && pp_config[j].Q_kz == Qz)
				{
					pp_config[j].flag = 0;
					times_R ++;
					//hhpp_channel_R[hhpp_temp_dimension].times++;
				}
			}
            if(times_R != 0)
            {
            	hhpp_channel_L[hhpp_temp_dimension].Qx = Qx;
				hhpp_channel_L[hhpp_temp_dimension].Qy = Qy;
				hhpp_channel_L[hhpp_temp_dimension].Qz = Qz;
				hhpp_channel_L[hhpp_temp_dimension].times = times_L;
            			
            	hhpp_channel_R[hhpp_temp_dimension].Qx = Qx;
				hhpp_channel_R[hhpp_temp_dimension].Qy = Qy;
				hhpp_channel_R[hhpp_temp_dimension].Qz = Qz;
				hhpp_channel_R[hhpp_temp_dimension].times = times_R;
            }
            else
            {
            	continue; //如果times_R =0,即找不到匹配的pp和hh，那么不进行后面的分配，直接进行下一轮循环
            }

		//	cout<<"here = "<<hhpp_channel_L[hhpp_temp_dimension].times<<endl;

			//为每一个hhpp_channel[hh_temp_dimension]动态分配
			tbwf *temp_star_hhpp_L;
			temp_star_hhpp_L = (tbwf *)malloc(hhpp_channel_L[hhpp_temp_dimension].times * sizeof(tbwf));
			hhpp_channel_L[hhpp_temp_dimension].wf = temp_star_hhpp_L;
			
			tbwf *temp_star_hhpp_R;
			temp_star_hhpp_R = (tbwf *)malloc(hhpp_channel_R[hhpp_temp_dimension].times * sizeof(tbwf));
			hhpp_channel_R[hhpp_temp_dimension].wf = temp_star_hhpp_R;		
	
			hhpp_temp_dimension++; 
		}
	}


	for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
	{
		hh_config[i].flag = 1;
	}
	for(int i = 0; i<hh_no; i++) //将hh_config[i]的flag全都赋值回1
	{
		pp_config[i].flag = 1;
	}
//	for(int i = 0; i<hh_dimension; i++)
//	{
	//	cout<<"hh dimension "<<i<<endl;
	//	cout<<"times = "<<hh_channel[i].times<<" "<<hh_channel[i].Qx<<" "<<hh_channel[i].Qy<<" "<<hh_channel[i].Qz<<endl;
//	}

//	cout<<"hh no"<<hh_no<<endl;
	//////////////////////////////////////////
	//**************************************//
	///////////重要判定！！！！！////////////////
	int xx;
	if(hh_dimension != hhpp_dimension)
	{
		cout<<"ERROR!"<<endl;
		cin>>xx;
	}
	//**************************************//
	//////////////////////////////////////////

    for(int i = 0; i<hhpp_dimension; i++) //循环所有hhpp_channel 并且对hhpp_channel[].L_wf R_wf 赋值
    {
    	k = 0;
    	//cout<<"channel = "<<i<<endl;
    	//cout <<" qx="<<hhpp_channel[i].Qx<<" qy="<<hhpp_channel[i].Qy<<" qz="<<hhpp_channel[i].Qz<<endl;

    	for(int j = 0; j< hh_no; j++)
    	{

    		//cout <<" hqx="<<hh_config[j].Q_kx<<" hqy="<<hh_config[j].Q_ky<<" hqz="<<hh_config[j].Q_kz<<endl;
    		
    		if(hh_config[j].Q_kx == hhpp_channel_L[i].Qx && hh_config[j].Q_ky == hhpp_channel_L[i].Qy && hh_config[j].Q_kz == hhpp_channel_L[i].Qz)
			{
				//cout <<" hqx="<<hh_config[j].Q_kx<<" hqy="<<hh_config[j].Q_ky<<" hqz="<<hh_config[j].Q_kz<<endl;
    		   // cout <<" hqkx="<<hh_config[j].q_kx<<" hqky="<<hh_config[j].q_ky<<" hqkz="<<hh_config[j].q_kz<<endl;
				hhpp_channel_L[i].wf[k].q_kx = hh_config[j].q_kx;
				hhpp_channel_L[i].wf[k].q_ky = hh_config[j].q_ky;
			    hhpp_channel_L[i].wf[k].q_kz = hh_config[j].q_kz;
				hhpp_channel_L[i].wf[k].sz1  = hh_config[j].sz1;
				hhpp_channel_L[i].wf[k].sz2  = hh_config[j].sz2;
                k ++ ;  
                //cout <<"erro"<<endl;  
			}
	
			
    	}	   
    	 if (k != hhpp_channel_L[i].times)
			{
				cout<<"wtf erro1#"<<endl;
			}
    	//cout<< "k="<<k<<endl;
    	
    	
	    //cout <<"erro"<<endl;
        k = 0;
    	for(int j = 0; j< pp_no; j++)
    	{
    		if(pp_config[j].Q_kx == hhpp_channel_R[i].Qx && pp_config[j].Q_ky == hhpp_channel_R[i].Qy && pp_config[j].Q_kz == hhpp_channel_R[i].Qz)
			{
				hhpp_channel_R[i].wf[k].q_kx = pp_config[j].q_kx;
				hhpp_channel_R[i].wf[k].q_ky = pp_config[j].q_ky;
				hhpp_channel_R[i].wf[k].q_kz = pp_config[j].q_kz;
				hhpp_channel_R[i].wf[k].sz1  = pp_config[j].sz1;
				hhpp_channel_R[i].wf[k].sz2  = pp_config[j].sz2;
                k ++ ; 
			}
    	}
    	if (k != hhpp_channel_R[i].times)
		{
			cout<< " k = "<<k<<" hhpp_times="<<hhpp_channel_R[i].times<<endl; 
			cout<<"wtf erro2#"<<endl;
		}
    }



	///////////////////////////////////////////////////////////
	//*******************************************************//
	//************** APPLY FOR THE V_ij_ab MATRIX ***********//
	//*******************************************************//
	///////////////////////////////////////////////////////////
	//这里先不填 V,t,H，后面填

	MATRIX_CHANNEL V_ijab[hhpp_dimension];
	MATRIX_CHANNEL t_ijab[hhpp_dimension];
	MATRIX_CHANNEL H_ijab[hhpp_dimension];



	for(int i = 0; i< hhpp_dimension; i++)
	{
		temp_pointer1_V = (double **)malloc(hhpp_channel_L[i].times * sizeof(double*));
		V_ijab[i].matrix = temp_pointer1_V;

		temp_pointer1_V = (double **)malloc(hhpp_channel_L[i].times * sizeof(double*));
		t_ijab[i].matrix = temp_pointer1_V;

		temp_pointer1_V = (double **)malloc(hhpp_channel_L[i].times * sizeof(double*));
		H_ijab[i].matrix = temp_pointer1_V;

		for(int ij = 0; ij < hhpp_channel_L[i].times; ij++)
		{
			temp_pointer2_V = (double *)malloc(hhpp_channel_R[i].times * sizeof(double));
			V_ijab[i].matrix[ij] = temp_pointer2_V;

			temp_pointer2_V = (double *)malloc(hhpp_channel_R[i].times * sizeof(double));
			t_ijab[i].matrix[ij] = temp_pointer2_V;

			temp_pointer2_V = (double *)malloc(hhpp_channel_R[i].times * sizeof(double));
			H_ijab[i].matrix[ij] = temp_pointer2_V;
		}
	}
	


	////////////////////////////////////////////////////////////////////////
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	//*************************  PP-PP- CHANNEL  *************************//
	//********************************************************************//
	//********************************************************************//
	//********************************************************************//
	////////////////////////////////////////////////////////////////////////


	for(int i = 0; i<pp_no; i++)
	{
		pp_config[i].flag = 1;
	}

	CHANNEL pp_channel[hhpp_dimension]; //先声明！！！！！！！

	//int pp_temp_dimension = 0;
	for(int i = 0; i< hhpp_dimension; i++)
	{
		pp_channel[i].Qx = hhpp_channel_L[i].Qx;
		pp_channel[i].Qy = hhpp_channel_L[i].Qy;
		pp_channel[i].Qz = hhpp_channel_L[i].Qz;
		pp_channel[i].times = hhpp_channel_R[i].times;

		//为每一个pp_channel[pp_temp_dimension]动态分配
		tbwf *temp_star;
		temp_star = (tbwf *)malloc(pp_channel[i].times * sizeof(tbwf));
		pp_channel[i].wf = temp_star;
	}



	for(int i = 0; i< hhpp_dimension; i++)
	{
		for(int j = 0; j < pp_channel[i].times; j++)
		{
			
			pp_channel[i].wf[j].q_kx = hhpp_channel_R[i].wf[j].q_kx;
			pp_channel[i].wf[j].q_ky = hhpp_channel_R[i].wf[j].q_ky;
			pp_channel[i].wf[j].q_kz = hhpp_channel_R[i].wf[j].q_kz;
			pp_channel[i].wf[j].sz1  = hhpp_channel_R[i].wf[j].sz1;
			pp_channel[i].wf[j].sz2  = hhpp_channel_R[i].wf[j].sz2;			
		}
//			if(k != pp_channel[loop5].times) 
//			{
//					cout<<"k = "<<k<<" , times = "<<pp_channel[loop5].times<<endl;
//			}
	}

//	cout<<loop5<<endl;
//	cout<<"test count "<<count_1<<endl;

	//// Assignment for the V_cdab ////
	//为每一个pp_channel[pp_temp_dimension]动态分配

	MATRIX_CHANNEL V_cdab[hhpp_dimension];

	for(int i = 0; i< hhpp_dimension; i++)
	{
		temp_pointer1_V = (double **)malloc(pp_channel[i].times * sizeof(double*));
		V_cdab[i].matrix = temp_pointer1_V;

		for(int cd = 0; cd < pp_channel[i].times; cd++)
		{
			temp_pointer2_V = (double *)malloc(pp_channel[i].times * sizeof(double));
			V_cdab[i].matrix[cd] = temp_pointer2_V;
		}
	}
	



////////////// FULFILL THE V_cd_ab MATRIX ////////////
	for(int i = 0; i< hhpp_dimension; i++)
	{
//		cout<<"Q = "<<pp_channel[i].Qx<<" "<<pp_channel[i].Qy<<" "<<pp_channel[i].Qz;
//		cout<<" times = "<<pp_channel[i].times<<endl;
		V_cdab[i].Qx = pp_channel[i].Qx;
		V_cdab[i].Qy = pp_channel[i].Qy;
		V_cdab[i].Qz = pp_channel[i].Qz;

		for(int ab = 0; ab < pp_channel[i].times; ab++)
		{
			for(int cd = 0; cd < pp_channel[i].times; cd++)
			{
				V_cdab[i].matrix[cd][ab]= V_ks_AS(L, pp_channel[i].wf[cd].q_kx, pp_channel[i].wf[cd].q_ky, pp_channel[i].wf[cd].q_kz, pp_channel[i].wf[cd].sz1, pp_channel[i].wf[cd].sz2,\
				pp_channel[i].wf[ab].q_kx,pp_channel[i].wf[ab].q_ky, pp_channel[i].wf[ab].q_kz, pp_channel[i].wf[ab].sz1, pp_channel[i].wf[ab].sz2);
				//cout<<V_cdab[i].matrix[cd][ab]<<" ";
			}
//			cout<<"\n";
		}
//		cout<<endl;
//		cout<<"///////////////"<<endl;
//		cout<<endl;
	}







	/////////////////////////////////////////////////////////////////////////
	//*********************************************************************//
	//******************** FULFILL THE V_ij_ab MATRIX *********************//
	//********************************************************************///
	/////////////////////////////////////////////////////////////////////////

	for(int i = 0; i< hhpp_dimension; i++)
	{
//		cout<<"Q = "<<pp_channel[i].Qx<<" "<<pp_channel[i].Qy<<" "<<pp_channel[i].Qz;
//		cout<<" times = "<<pp_channel[i].times<<endl;
		V_ijab[i].Qx = hhpp_channel_L[i].Qx;
		V_ijab[i].Qy = hhpp_channel_L[i].Qy;
		V_ijab[i].Qz = hhpp_channel_L[i].Qz;

		for(int ij = 0; ij < hhpp_channel_L[i].times; ij++)
		{
			for(int ab = 0; ab < hhpp_channel_R[i].times; ab++)
			{
				V_ijab[i].matrix[ij][ab]= V_ks_AS(L, hhpp_channel_L[i].wf[ij].q_kx, hhpp_channel_L[i].wf[ij].q_ky, hhpp_channel_L[i].wf[ij].q_kz, hhpp_channel_L[i].wf[ij].sz1, hhpp_channel_L[i].wf[ij].sz2,\
													 hhpp_channel_R[i].wf[ab].q_kx, hhpp_channel_R[i].wf[ab].q_ky, hhpp_channel_R[i].wf[ab].q_kz, hhpp_channel_R[i].wf[ab].sz1, hhpp_channel_R[i].wf[ab].sz2);
	//			cout<<V_ijab[i].matrix[ij][ab]<<" ";
			}
	//		cout<<"\n";
		}
	//	cout<<endl;
	//	cout<<"///////////////"<<endl;
	//	cout<<endl;
	}



	/////////////////////////////////////////////////////////////////////////
	//*********************************************************************//
	//******************** FULFILL THE t_ij_ab MATRIX *********************//
	//********************************************************************///
	/////////////////////////////////////////////////////////////////////////

	int kp_x, kp_y, kp_z;
	int sp;
    double f_ii,f_jj,f_aa,f_bb;

	for(int i = 0; i< hhpp_dimension; i++)
	{
//		cout<<"Q = "<<pp_channel[i].Qx<<" "<<pp_channel[i].Qy<<" "<<pp_channel[i].Qz;
//		cout<<" times = "<<pp_channel[i].times<<endl;
		t_ijab[i].Qx = hhpp_channel_L[i].Qx;
		t_ijab[i].Qy = hhpp_channel_L[i].Qy;
		t_ijab[i].Qz = hhpp_channel_L[i].Qz;

		for(int ij = 0; ij < hhpp_channel_L[i].times; ij++)
		{
			kp_x = (int)(((double)hhpp_channel_L[i].Qx + 2 * hhpp_channel_L[i].wf[ij].q_kx) / 2.);
			kp_y = (int)(((double)hhpp_channel_L[i].Qy + 2 * hhpp_channel_L[i].wf[ij].q_ky) / 2.);
			kp_z = (int)(((double)hhpp_channel_L[i].Qz + 2 * hhpp_channel_L[i].wf[ij].q_kz) / 2.);
			sp = hhpp_channel_L[i].wf[ij].sz1;
            f_ii = f_ks(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp);

            kp_x = (int)(((double)hhpp_channel_L[i].Qx - 2 * hhpp_channel_L[i].wf[ij].q_kx) / 2.);
			kp_y = (int)(((double)hhpp_channel_L[i].Qy - 2 * hhpp_channel_L[i].wf[ij].q_ky) / 2.);
			kp_z = (int)(((double)hhpp_channel_L[i].Qz - 2 * hhpp_channel_L[i].wf[ij].q_kz) / 2.);
			sp = hhpp_channel_L[i].wf[ij].sz2;
            f_jj = f_ks(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp);
			
			for(int ab = 0; ab < hhpp_channel_R[i].times; ab++)
			{    				
                kp_x = (int)(((double)hhpp_channel_R[i].Qx + 2 * hhpp_channel_R[i].wf[ab].q_kx) / 2.);
				kp_y = (int)(((double)hhpp_channel_R[i].Qy + 2 * hhpp_channel_R[i].wf[ab].q_ky) / 2.);
				kp_z = (int)(((double)hhpp_channel_R[i].Qz + 2 * hhpp_channel_R[i].wf[ab].q_kz) / 2.);
				sp = hhpp_channel_R[i].wf[ab].sz1;
                f_aa = f_ks(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp);

                kp_x = (int)(((double)hhpp_channel_R[i].Qx - 2 * hhpp_channel_R[i].wf[ab].q_kx) / 2.);
				kp_y = (int)(((double)hhpp_channel_R[i].Qy - 2 * hhpp_channel_R[i].wf[ab].q_ky) / 2.);
				kp_z = (int)(((double)hhpp_channel_R[i].Qz - 2 * hhpp_channel_R[i].wf[ab].q_kz) / 2.);
				sp = hhpp_channel_R[i].wf[ab].sz2;
                f_bb = f_ks(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp);
            

				t_ijab[i].matrix[ij][ab]= V_ks_AS(L, hhpp_channel_L[i].wf[ij].q_kx, hhpp_channel_L[i].wf[ij].q_ky, hhpp_channel_L[i].wf[ij].q_kz, hhpp_channel_L[i].wf[ij].sz1, hhpp_channel_L[i].wf[ij].sz2,\
													 hhpp_channel_R[i].wf[ab].q_kx, hhpp_channel_R[i].wf[ab].q_ky, hhpp_channel_R[i].wf[ab].q_kz, hhpp_channel_R[i].wf[ab].sz1, hhpp_channel_R[i].wf[ab].sz2) \
											/ (f_ii+f_jj-f_aa-f_bb);
		//		cout<<t_ijab[i].matrix[ij][ab]<<" ";
			}
		//	cout<<"\n";
		}
	//	cout<<endl;
	//	cout<<"///////////////"<<endl;
	//	cout<<endl;
	}


	/////////////////////////////////////////////////////////////////////////
	//*********************************************************************//
	//******************** FULFILL THE H_ij_ab MATRIX *********************//
	//********************************************************************///
	/////////////////////////////////////////////////////////////////////////
	int iter = 0;
	double t_no_this, t_last;
	double alpha = 0.9; //*** Here we could try different alpha!!! ****//
	double temp1, temp2, delta;
	double Ec = 0.;
	do
	{	
		temp1 = Ec;
		Ec = 0.;
		double h0, h1, h2, h3, h4;
		for(int i = 0; i< hhpp_dimension; i++)
		{
//		cout<<"Q = "<<pp_channel[i].Qx<<" "<<pp_channel[i].Qy<<" "<<pp_channel[i].Qz;
//		cout<<" times = "<<pp_channel[i].times<<endl;
			H_ijab[i].Qx = hhpp_channel_L[i].Qx;
			H_ijab[i].Qy = hhpp_channel_L[i].Qy;
			H_ijab[i].Qz = hhpp_channel_L[i].Qz;

			for(int ij = 0; ij < hhpp_channel_L[i].times; ij++)
			{
				kp_x = (int)(((double)hhpp_channel_L[i].Qx + 2 * hhpp_channel_L[i].wf[ij].q_kx) / 2.);
				kp_y = (int)(((double)hhpp_channel_L[i].Qy + 2 * hhpp_channel_L[i].wf[ij].q_ky) / 2.);
				kp_z = (int)(((double)hhpp_channel_L[i].Qz + 2 * hhpp_channel_L[i].wf[ij].q_kz) / 2.);
				sp = hhpp_channel_L[i].wf[ij].sz1;
				f_ii = f_ks(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp);
				//cout<<"i "<<f_ii<<endl;//<<" si = "<<sp<<endl;
				double byhand;
				byhand = pow(h_bar,2) * pow((2.*pi/L),2) * (pow(kp_x,2)+pow(kp_y,2)+pow(kp_z,2))/(2*mass);
				byhand += 7. * (200./pow(L,3) * pow(pi/1.487,1.5) + (-91.85)/pow(L,3) * pow(pi/0.465,1.5));

				//cout<<"by hand = "<<byhand<<endl;

				kp_x = (int)(((double)hhpp_channel_L[i].Qx - 2 * hhpp_channel_L[i].wf[ij].q_kx) / 2.);
				kp_y = (int)(((double)hhpp_channel_L[i].Qy - 2 * hhpp_channel_L[i].wf[ij].q_ky) / 2.);
				kp_z = (int)(((double)hhpp_channel_L[i].Qz - 2 * hhpp_channel_L[i].wf[ij].q_kz) / 2.);
				sp = hhpp_channel_L[i].wf[ij].sz2;
				f_jj = f_ks(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp);
	
				for(int ab = 0; ab < hhpp_channel_R[i].times; ab++)
				{

            	    H_ijab[i].matrix[ij][ab] = 0.;

					h0 = V_ijab[i].matrix[ij][ab];
					h1 = 0.;
					h2 = 0.;
					for(int kl = 0; kl< hh_channel[i].times; kl++)
					{
						h1 += 0.5 * t_ijab[i].matrix[kl][ab] * V_ijkl[i].matrix[ij][kl];
					}

					for(int cd = 0; cd < pp_channel[i].times; cd++)
					{
						h2 += 0.5 * V_cdab[i].matrix[cd][ab] * t_ijab[i].matrix[ij][cd];
					}
			//	cout<<"error!"<<endl;


					kp_x = (int)(((double)hhpp_channel_R[i].Qx + 2 * hhpp_channel_R[i].wf[ab].q_kx) / 2.);
					kp_y = (int)(((double)hhpp_channel_R[i].Qy + 2 * hhpp_channel_R[i].wf[ab].q_ky) / 2.);
					kp_z = (int)(((double)hhpp_channel_R[i].Qz + 2 * hhpp_channel_R[i].wf[ab].q_kz) / 2.);
					sp = hhpp_channel_R[i].wf[ab].sz1;
					f_aa = f_ks(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp);

                	kp_x = (int)(((double)hhpp_channel_R[i].Qx - 2 * hhpp_channel_R[i].wf[ab].q_kx) / 2.);
					kp_y = (int)(((double)hhpp_channel_R[i].Qy - 2 * hhpp_channel_R[i].wf[ab].q_ky) / 2.);
					kp_z = (int)(((double)hhpp_channel_R[i].Qz - 2 * hhpp_channel_R[i].wf[ab].q_kz) / 2.);
					sp = hhpp_channel_R[i].wf[ab].sz2;
					f_bb = f_ks(config, spstates_no, L, kp_x, kp_y, kp_z, sp, kp_x, kp_y, kp_z, sp);

					h3 =  (f_aa + f_bb) * t_ijab[i].matrix[ij][ab];
					h4 = -(f_ii + f_jj) * t_ijab[i].matrix[ij][ab];

					H_ijab[i].matrix[ij][ab] = h0 + h1 + h2 + h3 + h4;

					//////////////////////////////
					//**************************//
					//**** ITERATIVE t_ijab ****//
					//**************************//
					//////////////////////////////
					///// RENEW t_ijab MATRIX ELEMENTS /////
				//////**** HERE we should use MIXING !!!****///////
			  //**** t(i) = a * t_no_mixing(i) + (1-a) * t(i-1) ****//
//**** ALTERNATIVE CHOICE, WHICH IS IDENTICAL t(i) = t(i-1) + a * H_bar / (f_aa + f_bb + f_ii + f_jj) *****//

					t_last = t_ijab[i].matrix[ij][ab];
					t_no_this = t_ijab[i].matrix[ij][ab] - H_ijab[i].matrix[ij][ab] / ( f_aa + f_bb - f_ii - f_jj );
					t_ijab[i].matrix[ij][ab] = alpha * t_no_this + (1-alpha) * t_last;
			
				//	cout<<H_ijab[i].matrix[ij][ab]<<" ";
             //   cout<<"f aa bb ii jj "<<f_aa<<" "<<f_bb<<" "<<f_ii<<" "<<f_jj<<endl;
				}
			//	cout<<"\n";
			}
		//	cout<<endl;
		//	cout<<"///////////////"<<endl;
		//	cout<<endl;
		}

		////////////////////////////////
		//****************************//
		//**** correlation energy ****//
		//****************************//
		////////////////////////////////
		for(int i = 0; i< hhpp_dimension; i++)
		{
			for(int ij = 0; ij < hhpp_channel_L[i].times; ij++)
			{
				for(int ab = 0; ab < hhpp_channel_R[i].times; ab++)
				{
					Ec += 0.25 * t_ijab[i].matrix[ij][ab] * V_ijab[i].matrix[ij][ab];

				}
			}
		}
		
		temp2 = Ec; //record the Ec this time
		
		delta = fabs(temp1 - temp2);
		
		//cout<<"Ec/N_"<<iter<<" = "<<temp2/14.<<", d_Ec = "<<delta<<endl;
		
		iter++; //COUNT THE ITERATION TIMES

		if(iter>1000) break; // IF IT CANNOT CONVERGE UNTIL 1000 ITERATIONS, BREAK OUT.

	}while(delta>10e-7);
        cout.setf(ios::right, ios::adjustfield);
        cout<<" , (Ec/N) ="  << setw(7)<< setprecision(5) <<temp2/magic_no;
        cout<<" , (Etot/N) ="<< setw(8)<<(hf_energy+temp2)/magic_no;



	////////////////////////////////
	//****************************//
	//*********** FREE ***********//
	//****************************//
	////////////////////////////////

	free(config);
	free(pp_config);
	free(hh_config);
//	free(ph_config);

	for(int i = 0; i<hhpp_dimension; i++)
	{
		free(pp_channel[i].wf);
	}
	for(int i = 0; i<hh_dimension; i++)
	{
		free(hh_channel[i].wf);
	}
	for(int i = 0; i<hhpp_dimension; i++)
	{
		free(hhpp_channel_L[i].wf);
		free(hhpp_channel_R[i].wf);
	}

	//FREE THE V MATRIX//
	for(int i = 0; i<hhpp_dimension; i++)
	{
		for(int j = 0; j<pp_channel[i].times; j++)
		{
			free(V_cdab[i].matrix[j]);			
		}
		free(V_cdab[i].matrix);
	}
	//FREE THE V MATRIX//
	for(int i = 0; i<hh_dimension; i++)
	{
		for(int j = 0; j<hh_channel[i].times; j++)
		{
			free(V_ijkl[i].matrix[j]);			
		}
		free(V_ijkl[i].matrix);

	}
	for(int i = 0; i<hhpp_dimension; i++)
	{
		for(int j = 0; j<hhpp_channel_L[i].times; j++)
		{
			free(V_ijab[i].matrix[j]);	
			free(t_ijab[i].matrix[j]);
			free(H_ijab[i].matrix[j]);		
		}
		free(V_ijab[i].matrix);
		free(t_ijab[i].matrix);
		free(H_ijab[i].matrix);
	}


}

int main()
{
	double rho;
	clock_t start, end;
	for (int i = 0; i < 9; i++ )
	{
		rho = 0.04 + double(i)/50.;
		start = clock();
		CCD_neutronmatter(rho);
		end = clock();
		cout<<" , TIME = "<<(double)(end-start)/CLOCKS_PER_SEC<<" s"<<endl;
	}

}
