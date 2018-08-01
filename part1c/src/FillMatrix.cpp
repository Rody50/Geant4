
#include "FillMatrix.h"
#include "Tensor2.h"
#include "Tensor4.h"
#include "CCMInf.h"


FillMatrix::FillMatrix()
{}

void FillMatrix::FillV(Tensor4 & v, double g)
{
	int iMax[4];

	for (int i = 0; i < 4; i++)
		iMax[i] = v.GetL(i);

	for (int a = 0; a < iMax[0]; a++)
		for (int b = 0; b < iMax[1]; b++)
			for (int i = 0; i < iMax[2]; i++)
				for (int j = 0; j < iMax[3]; j++)
				{
					if (IsPair(i, j) && IsPair(a, b) && i < j && a < b )
						{
							v(a, b, i, j) = - g / 2.;
							v(a, b, j, i) =   g / 2.;
							v(b, a, i, j) =   g / 2.;
							v(b, a, j, i) = - g / 2.;
						}	
				}
}

void FillMatrix::FillV(mat & vpppp, mat & vpphh,
	mat & vhhhh, P_group_t * p, InfiniteMatterSP *InfSP_p)
{
	vector<qState> &s = InfSP_p->GetfStates();
	int size = p->np;
	int hhhhA = 0, pphhA = 0, ppppA = 0;
	std::cout << "error? #0" << std::endl;

	for (int a = 0; a < size; a++) // loop over row
	{
		int hhhhI = 0, pphhI = 0, ppppI = 0;
		for (int i = 0; i < size; i++) //  loop over column
		{
				std::cout << "error? #1" << std::endl;

			qState * sa = & s[p->pr[a].i];
			qState * sb = & s[p->pr[a].j];
			qState * si = & s[p->pr[i].i];
			qState * sj = & s[p->pr[i].j];

			bool isHoleA = sa->isHole;
			bool isHoleB = sb->isHole;
			bool isHoleI = si->isHole;
			bool isHoleJ = sj->isHole;

			if (isHoleA && isHoleB && isHoleI && isHoleJ)
			{
				vhhhh(hhhhA + hhhhI * size) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
					std::cout << "error? #2" << std::endl;

				hhhhI++;
			}
			if ((!isHoleA && !isHoleB && isHoleI && isHoleJ) || 
				(isHoleA && isHoleB && !isHoleI && !isHoleJ))
			{	
				getchar();
				vpphh(pphhA + pphhI * size) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
					std::cout << "error? #3" << std::endl;

				pphhI++;
			}
			if (!isHoleA && !isHoleB && !isHoleI && !isHoleJ)
			{
				vpppp(ppppA + ppppI * size) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
					std::cout << "error? #4" << std::endl;

				ppppI++;
			}
		}
		if (hhhhI) hhhhA++; 
		if (pphhI) pphhA++; 
		if (ppppI) ppppA++;
	}
}

void FillMatrix::FillF(Tensor2 & f, int A, double g, double d)
{
	bool isHole = false;
	if(A == 0) isHole = true;

	for (int i = 0; i < f.GetL(1); i++)
	{
		if (!isHole) f(i, i) = (i + A) / 2 * d;

		if (isHole) f(i, i) = i / 2 * d - g / 2.;
	}
}

void FillMatrix::FillT(Tensor4 & t, Tensor4 & v, Tensor2 & fp,
	Tensor2 & fh)
{
	int iMax[4];

	for (int i = 0; i < 4; i++)
		iMax[i] = t.GetL(i);

	for (int a = 0; a < iMax[0]; a++)
		for (int b = 0; b < iMax[1]; b++)
			for (int i = 0; i < iMax[2]; i++)
				for (int j = 0; j < iMax[3]; j++)
				{
					if (IsPair(i, j) && IsPair(a, b) && i < j && a < b)
						{
							t(a, b, i, j) = v(a, b, i, j) / 
								(fh(i, i) + fh(j, j) - fp(a, a) - fp(b, b));
							t(a, b, j, i) = - t(a, b, i, j);
							t(b, a, i, j) = - t(a, b, i, j);
							t(b, a, j, i) =   t(a, b, i, j);
						}

				}
}

void FillMatrix::FillT(mat & tpphh,
	P_group_t * p, CCMInf *CCMInf_p)
{
	InfiniteMatterSP *InfSP_p = CCMInf_p->GetfInfSP();
	vector<qState> &s = InfSP_p->GetfStates();
	int size = p->np;
	int pphhA = 0;

	for (int a = 0; a < size; a++) // loop over row
	{
		int pphhI = 0;
		for (int i = 0; i < size; i++) //  loop over column
		{
			qState * sa = & s[p->pr[a].i];
			qState * sb = & s[p->pr[a].j];
			qState * si = & s[p->pr[i].i];
			qState * sj = & s[p->pr[i].j];

			bool isHoleA = sa->isHole;
			bool isHoleB = sb->isHole;
			bool isHoleI = si->isHole;
			bool isHoleJ = sj->isHole;

			if (!isHoleA && !isHoleB && isHoleI && isHoleJ)
			{	
				tpphh(pphhA + pphhI * size) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i])
					/ CCMInf_p->FSum(&p->pr[a], &p->pr[i]);
				pphhI++;
			}

		}
		if (pphhI) pphhA++; 	
	}	
}
// void FillMatrix::Print()
// {
// 	fConfigObj.Print();
// 	v->print();
// }

bool FillMatrix::IsPair(int a, int b)
{
	return floor(a / 2.) == floor(b / 2.);
}