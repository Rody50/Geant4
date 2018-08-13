
#include "FillMatrix.h"
#include "Tensor2.h"
#include "Tensor4.h"
#include "CCMInf.h"
#include "IMSRGInf.h"


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

	for (int a = 0; a < size; a++) // loop over row
	{
		int hhhhI = 0, pphhI = 0, ppppI = 0;
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

			if (isHoleA && isHoleB && isHoleI && isHoleJ)
			{
				vhhhh(hhhhA, hhhhI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				hhhhI++;
			}
			if ((!isHoleA && !isHoleB && isHoleI && isHoleJ))
			{
				vpphh(pphhA, pphhI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				pphhI++;
			}
			if (!isHoleA && !isHoleB && !isHoleI && !isHoleJ)
			{
				vpppp(ppppA, ppppI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				ppppI++;
			}
		}
		if (hhhhI) hhhhA++; 
		if (pphhI) pphhA++; 
		if (ppppI) ppppA++;
	}
}

void FillMatrix::FillV(mat & v,
	P_group_t *p, InfiniteMatterSP *InfSP_p)
{
	int npp = p->npp, nph = p->nph, nhh = p->nhh;

	mat vpppp(p->npp, p->npp, fill::zeros);
	mat vpphh(p->npp, p->nhh, fill::zeros);
	mat vhhhh(p->nhh, p->nhh, fill::zeros);
	mat vphph(p->nph, p->nph, fill::zeros);
	mat vhhph(p->nhh, p->nph, fill::zeros);
	mat vppph(p->npp, p->nph, fill::zeros);

	vector<qState> &s = InfSP_p->GetfStates();
	int size = p->np;
	int ppppA = 0, pphhA = 0, hhhhA = 0, phphA = 0, hhphA = 0, ppphA = 0;

	for (int a = 0; a < size; a++) // loop over row
	{
		int ppppI = 0, pphhI = 0, hhhhI = 0, phphI = 0, hhphI = 0, ppphI = 0;
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

			// pppp
			if (!isHoleA && !isHoleB && !isHoleI && !isHoleJ)
			{
				vpppp(ppppA, ppppI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				ppppI++;
			}
			// pphh
			if ((!isHoleA && !isHoleB && isHoleI && isHoleJ))
			{
				vpphh(pphhA, pphhI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				pphhI++;
			}
			// hhhh
			if (isHoleA && isHoleB && isHoleI && isHoleJ)
			{
				vhhhh(hhhhA, hhhhI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				hhhhI++;
			}
			// phph
			if (!isHoleA && isHoleB && !isHoleI && isHoleJ)
			{
				vphph(phphA, phphI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				phphI++;
			}
			// hhph
			if ((isHoleA && isHoleB && !isHoleI && isHoleJ))
			{
				vhhph(hhphA, hhphI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				hhphI++;
			}
			// ppph
			if (!isHoleA && !isHoleB && !isHoleI && isHoleJ)
			{
				vppph(ppphA, ppphI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i]);
				ppphI++;
			}
		}
		if (ppppI) ppppA++;
		if (pphhI) pphhA++; 
		if (hhhhI) hhhhA++; 
		if (phphI) phphA++; 
		if (hhphI) hhphA++; 
		if (ppphI) ppphA++;
	}
	v(span(0, nhh-1), span(0, nhh+nph+npp-1)) = 
		join_rows(join_rows(vhhhh, vhhph), vpphh.t());
	
	if (nph != 0)
		v(span(nhh, nhh+nph-1), span(0, nhh+nph+npp-1))  = 
			join_rows(join_rows(vhhph.t(), vphph),vppph.t());
	
	v(span(nhh+nph, nhh+nph+npp-1), span(0, nhh+nph+npp-1)) = 
		join_rows(join_rows(vpphh, vppph), vpppp);
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

void FillMatrix::FillF(mat & f,
			P_group_t *p, IMSRGInf *imsrg_p)
{
	InfiniteMatterSP *InfSP_p = imsrg_p->GetfInfSP();
	vector<qState> &s = InfSP_p->GetfStates();

	int npp = p->npp, nph = p->nph, nhh = p->nhh;
	int size = p->np;

	mat fpppp(p->npp, p->npp, fill::zeros);
	mat fhhhh(p->nhh, p->nhh, fill::zeros);
	mat fphph(p->nph, p->nph, fill::zeros);

	int ppppA = 0, hhhhA = 0, phphA = 0;

	for (int a = 0; a < size; a++) // loop over row
	{
		int ppppI = 0, hhhhI = 0, phphI = 0;
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

			// pppp
			if (!isHoleA && !isHoleB && !isHoleI && !isHoleJ)
			{
				fpppp(ppppA, ppppI) = imsrg_p->F(sa) + imsrg_p->F(sb);
				ppppI++;
			}
			// hhhh
			if (isHoleA && isHoleB && isHoleI && isHoleJ)
			{
				fhhhh(hhhhA, hhhhI) = imsrg_p->F(sa) + imsrg_p->F(sb);
				hhhhI++;
			}
			// phph
			if (!isHoleA && isHoleB && !isHoleI && isHoleJ)
			{
				fphph(phphA, phphI) = imsrg_p->F(sa) + imsrg_p->F(sb);
				phphI++;
			}
		}
		if (ppppI) ppppA++;
		if (hhhhI) hhhhA++; 
		if (phphI) phphA++; 
	}
	cout << "#1" << endl;
	getchar();
	f(span(0, nhh-1), span(0, nhh-1)) = fhhhh;
	f.print();
	cout << "#2" << endl;
	getchar();
	if (nph != 0)
		f(span(nhh, nhh+nph-1), span(nhh, nhh+nph-1))  = fphph;
	cout << "#3" << endl;
	getchar();
	f(span(nhh+nph, nhh+nph+npp-1), span(nhh+nph, nhh+nph+npp-1)) = fpppp;
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
		qState * sa = & s[p->pr[a].i];
		qState * sb = & s[p->pr[a].j];
		bool isHoleA = sa->isHole;
		bool isHoleB = sb->isHole;
		if(isHoleA || isHoleB) continue;
		int pphhI = 0;

		for (int i = 0; i < size; i++) //  loop over column
		{
			qState * si = & s[p->pr[i].i];
			qState * sj = & s[p->pr[i].j];
			bool isHoleI = si->isHole;
			bool isHoleJ = sj->isHole;

			if (isHoleI && isHoleJ)
			{	
				tpphh(pphhA, pphhI) = InfSP_p->Minnesota(&p->pr[a], &p->pr[i])
					/ (-CCMInf_p->FSum(&p->pr[a], &p->pr[i]));
				pphhI++;
			}

		}
		if (pphhI) pphhA++; 	
	}	
}

void FillMatrix::FillEta(mat & eta, mat & matUpdate, 
	P_group_t * p, IMSRGInf * imsrg_p)
{
	InfiniteMatterSP *InfSP_p = imsrg_p->GetfInfSP();
	vector<qState> &s = InfSP_p->GetfStates();
	mat eta_temp(p->npp, p->nhh, fill::zeros);

	int size = p->np, nhh = p->nhh, npp = p->npp, nph = p->nph;
	int pphhA = 0;

	for (int a = 0; a < size; a++) // loop over row
	{
		qState * sa = & s[p->pr[a].i];
		qState * sb = & s[p->pr[a].j];
		bool isHoleA = sa->isHole;
		bool isHoleB = sb->isHole;
		
		if(isHoleA || isHoleB) continue;
		int pphhI = 0;

		for (int i = 0; i < size; i++) //  loop over column
		{
			qState * si = & s[p->pr[i].i];
			qState * sj = & s[p->pr[i].j];
			bool isHoleI = si->isHole;
			bool isHoleJ = sj->isHole;

			if (isHoleI && isHoleJ)
			{
				eta_temp(pphhA, pphhI) = matUpdate(pphhA, pphhI) 
					/ (imsrg_p->FSum(&p->pr[a], &p->pr[i]));
				pphhI++;
			}
		}
		if (pphhI) pphhA++; 
	}

	eta(span(0, nhh-1), span(nhh+nph, nhh+nph+npp-1)) = eta_temp.t();
	
	eta(span(nhh+nph, nhh+nph+npp-1), span(0, nhh-1)) = eta_temp;
}

bool FillMatrix::IsPair(int a, int b)
{
	return floor(a / 2.) == floor(b / 2.);
}
