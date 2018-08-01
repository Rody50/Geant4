
#include <iostream>

void Permute(int * i1, int * i2, int ind0, int ind1)
{
	for(int i = 0; i < 4; i++)
	{
		if(i1[i] == ind0) i1[i] = ind1;
		else if(i1[i] == ind1) i1[i] = ind0;
		
		if(i2[i] == ind0) i2[i] = ind1;
		else if(i2[i] == ind1) i2[i] = ind0;
	}
}

 class A{
 public:
 	A(){
 		fC[0] = 'a';
 		fC[1] = 'b';
 		fC[2] = 'c';
 		fC[3] = 'd';
 	}

 	char fC[5];
	void operator()(int a, int b, int c, int d){
		std::cout << fC[a] << " " << fC[b];
		std::cout << fC[c] << " " << fC[d] << std::endl;
	}

 };

int main()
{
	int i1[4], i2[4], i3[4], i4[4], i5[4];

	i1[0] = 2; i1[1] = 0; i1[2] = 3; i1[3] = 4;
	i2[0] = 3; i2[1] = 4; i2[2] = 1; i2[3] = 2;
	i3[0] = 0; i3[1] = 1; i3[2] = 4; i3[3] = 3;
	i5[0] = 4; i5[1] = 2; i5[2] = 9; i5[3] = 9;

	Permute(i2, i5, 2, 3);

//	for (int i = 0; i < 4; i++)
//		std::cout << i2[i] << " ";
//	std::cout << std::endl;
//	for (int i = 0; i < 4; i++)
//		std::cout << i5[i] << " ";
//		std::cout << std::endl;

	A a;
	int i = 0, j = 1, k = 2, l = 3;
	a(i, j, k, l);

}