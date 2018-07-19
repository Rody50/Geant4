#include <iostream>
#include <math.h>
using namespace std;

int count1Bits(int N) //Calculate the "1" counts of an integer written in form of binary form. 求一个整数写成二进制后共有多少个1，即全部核子数
{
	int n=0;
	do{
		N &= N-1;
		n++;
	}while(N!=0);
	return n;
}

int factorial(int n) // calculate the factorial
{
	int result = n;
	do
	{
		result *= (n-1);
		n--;
	} while(n != 1);

	return result;
}

int combinator(int n, int m) //calculate the combinator
{
	int result;
	result = factorial(n) / (factorial(m) * factorial(n-m));
	return result;
}

void print_bin(int N, int n_max) //the bin-representation, type an integer in form of binary form. 以二进制形式打印一个整数，即显示成bit-representation形式的组态
{
	double f = 0;
	int digits = 0;
	f = frexp(n_max, &digits);

    for(int i = digits-1; i >= 0; i--)
    {
        if( ( N & (1<<i) ) != 0) cout<<"1";
        else cout<<"0";
    }
    cout<<endl;
}

/* Creat the configuration space. The reference state correspond to the N_max. All the other states are beyond */
/* this state, and at the same time the N values are smaller than N_max. So we use N-- here.*/

void creat_config(int particle_no, int degen_lev)//creat the full configurations.生成全部组态空间
{
	int n = combinator(2*degen_lev, particle_no);
	int *config = (int *) malloc (n*sizeof(int));
	int N_max;//N_max is corresponding value of the reference state (phi_0)
	N_max = pow(2, 2*degen_lev) - pow(2, 2*degen_lev - particle_no);//equals to 2^(2*P)-2^(2P-A). according to the geometric series sum formula.
	int N = N_max;
	int i=0;
	int number_1;

	do
	{
		number_1 = count1Bits(N);
		if(number_1==particle_no) // find the configurations have the same "1" occupied numbers who equal A.
		{
			config[i] = N;
			N--;
			i++;
		}
		else
		{
			N--;
		}
	}while(N != 0);

	free(config);

	for(int i = 0; i < n ; i++)
	{
		print_bin(config[i], config[0]); //print the configuration in form of binary form.
	}
}

int main()
{
	creat_config(4,4);
	return 0;
}
