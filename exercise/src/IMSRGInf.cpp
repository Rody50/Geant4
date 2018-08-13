#include "IMSRGInf.h"
#include <fstream>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::ofstream;
using std::to_string;
using std::flush;

IMSRGInf::IMSRGInf(P_group_t *p, InfiniteMatterSP *infSP, double deltaS) : 
	fA(p->nhh), fNu(p->npp), fP(p), fDeltaS(deltaS),
	fHamil0(p->np, p->np, fill::zeros),
	fHamilS(p->np, p->np, fill::zeros),
	fF(p->np, p->np, fill::zeros),
	fV(p->np, p->np, fill::zeros),
	fEta(p->np, p->np, fill::zeros),
	fOmega(p->np, p->np, fill::zeros),		
	fInfSP(infSP)
{
	cout << "npp: " << p->npp << "\tnhh: " << p->nhh << "\tnph: " << p->nph
		 << "\tnp: " << p->np << endl;
	FillMatrix::FillV(fV, p, fInfSP);
	FillMatrix::FillEta(fEta, fV, p, this);
	FillMatrix::FillF(fF, p, this);
	//cout << "The eta matrix: " << endl;
	//fEta.print();
	fHamil0 = fF + fV;
	fHamilS = fHamil0;
	//cout << "The omega matrix:" << endl;
	//fOmega.print();
}

inline int r2(int *p){
	return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
}

inline double fact(int n)
{
	int s = 1;
	for (int i = 1; i <= n; i++)
	{
		s *= i;
	}
	return s;
}

double IMSRGInf::F(qState *s)
{	
	vector<qState> &states = fInfSP->GetfStates();

	int n[3] = {s->nx, s->ny, s->nz};

	const double u = fInfSP->GetfEnUn();

	double v = 0.;

	for(auto &t : states){
		if(t.isHole){
			v += fInfSP->Minnesota(*s, t, *s, t);
		}
	}
	return r2(n) * u + v;
}

double IMSRGInf::FSum(pair_t *tpp, pair_t *shh)
{
	vector<qState> &states = fInfSP->GetfStates();
	qState *a = &states[tpp->i];
	qState *b = &states[tpp->j];
	qState *i = &states[shh->i];
	qState *j = &states[shh->j];

	int na[3] = {a->nx, a->ny, a->nz};
	int nb[3] = {b->nx, b->ny, b->nz};
	int ni[3] = {i->nx, i->ny, i->nz};
	int nj[3] = {j->nx, j->ny, j->nz};
	const double u = fInfSP->GetfEnUn();

	double va = 0., vb = 0., vi = 0., vj = 0.;

	for(auto &t : states){
		if(t.isHole){
			va += fInfSP->Minnesota(*a, t, *a, t);
			vb += fInfSP->Minnesota(*b, t, *b, t);
			vi += fInfSP->Minnesota(*i, t, *i, t);
			vj += fInfSP->Minnesota(*j, t, *j, t);
		}
	}

	double fa = r2(na) * u + va;
	double fb = r2(nb) * u + vb;
	double fi = r2(ni) * u + vi;
	double fj = r2(nj) * u + vj;

	return fa + fb - fi - fj;
}

mat IMSRGInf::Commute(const mat &a, const mat &b)
{
	return a * b - b * a;
}

void IMSRGInf::UpdateOm()
{
	fOmega += fEta + fEta * fDeltaS 
		+ 1 / 2 * Commute(fEta * fDeltaS, fOmega)
		+ 1 / 12 * (Commute(fEta * fDeltaS, Commute(fEta * fDeltaS, fOmega))
		+ Commute( fOmega, Commute(fOmega, fEta * fDeltaS)))
		- 1 / 24 * Commute(fOmega, Commute( fEta * fDeltaS,
		Commute(fEta * fDeltaS, fOmega)));
}

void IMSRGInf::UpdateEta()
{
	int nhh = fP->nhh, npp = fP->npp, nph = fP->nph;

	fEta(span(0, nhh-1), span(nhh+nph, nhh+nph+npp-1)) = 
		fHamilS(span(0, nhh-1), span(nhh+nph, nhh+nph+npp-1));
	fEta(span(nhh+nph, nhh+nph+npp-1), span(0, nhh-1)) =
		fHamilS(span(nhh+nph, nhh+nph+npp-1), span(0, nhh-1));
}

double IMSRGInf::Iterate()
{
	double corr_en;

	while (norm(fEta, "fro") > 1E-5)
	{	
		mat temp = fHamilS;
		int n = 0;
		UpdateOm();
		while (norm(temp, "fro") > 0.5)
		{
			temp = Commute(fOmega, temp);
			fHamilS += 1 / fact(n) * temp;
			n++;
		}
		UpdateEta();
	}

	fHamilS.print();

	// int n0 = 0;
	// mat temp0 = fHamil0;
	// while (norm(temp0, "fro") > 0.5)
	// {
	// 	temp = Commute(fOmega, temp0);
	// 	fHamil0 += 1 / fact(n0) * temp0;
	// 	n0++;
	// }	
	// Sum(fHamil0);

	return 0.;
}
