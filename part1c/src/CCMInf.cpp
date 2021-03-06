#include "CCMInf.h"
#include <fstream>
#include <iostream>
#include <string>

using std::ofstream;
using std::to_string;
using std::flush;


CCMInf::CCMInf(P_group_t *p, InfiniteMatterSP *infSP) : fN(p->npp),
	fA(p->nhh), fNu(p->npp), fP(p),
	fHamil(p->npp, p->nhh, fill::zeros),
	fVpppp(p->npp, p->npp, fill::zeros),	
	fVpphh(p->npp, p->nhh, fill::zeros), 
	fVhhpp(p->nhh, p->npp, fill::zeros), 
	fVhhhh(p->nhh, p->nhh, fill::zeros),
	fTpphh(p->npp, p->nhh, fill::zeros),
	fInterm(p->npp, p->nhh, fill::zeros),
	fInfSP(infSP)
{
	FillMatrix::FillV(fVpppp, fVpphh, fVhhhh,
		p, fInfSP);
	FillMatrix::FillT(fTpphh, p, this);
}

inline int r2(int *p){
	return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
}

double CCMInf::FSum(pair_t *tpp, pair_t *shh)
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
	

#ifdef FSUM
	//if(fabs(fabs(-fa - fb + fi + fj) - 480.) < 50){
		cout << "fa: " << fa << "\tfb: " << fb << "\tfi: " << fi << "\tfj: " << fj << endl;
		cout << "va: " << va << "\tvb: " << vb << "\tvi: " << vi << "\tvj: " << vj << endl;
		cout << "-fa - fb + fi + fj: " << -fa - fb + fi + fj << endl;
		getchar();
	//}
#endif

	return fa + fb - fi - fj;
}

void CCMInf::Interm()
{
	vector<qState> &s = fInfSP->GetfStates();
	int size = fP->np;
	int pphhA = 0;

	for (int a = 0; a < size; a++) // loop over row
	{
		qState * sa = & s[fP->pr[a].i];
		qState * sb = & s[fP->pr[a].j];
		
		bool isHoleA = sa->isHole;
		bool isHoleB = sb->isHole;
		if(isHoleA || isHoleB) continue;
		int pphhI = 0;
		for (int i = 0; i < size; i++) //  loop over column
		{

			qState * si = & s[fP->pr[i].i];
			qState * sj = & s[fP->pr[i].j];
			bool isHoleI = si->isHole;
			bool isHoleJ = sj->isHole;

			if (isHoleI && isHoleJ)
			{	
				fInterm(pphhA, pphhI) = fTpphh(pphhA, pphhI)
					* FSum(&(fP->pr[a]), &(fP->pr[i]));
				pphhI++;
			}

		}
		if (pphhI) pphhA++;
	}
}

void CCMInf::ComputeH()
{
	Interm();

	fHamil = fVpphh + fInterm + 0.5 * (fVpppp * fTpphh + fTpphh * fVhhhh);
}

double CCMInf::SolveT()
{
	double factor = 1.;

	double corr_en_pre = 1.;
	double corr_en = 0.;
	vector<double> correlationEn;

	int NStepsMax = 1E5;
	int NSteps = 0;

//	ofstream corr_en_file(("results/corr_en-steps" + to_string(fN) 
//		+ "_" + to_string(fA) + ".txt").c_str());

	vector<qState> &s = fInfSP->GetfStates();

	while (fabs(corr_en - corr_en_pre) > 1E-3 || NSteps > NStepsMax)
	{	
		ComputeH();
		
		corr_en_pre = corr_en;
		corr_en = 0.;

		int size = fP->np;
		int pphhA = 0;

		for (int a = 0; a < size; a++) // loop over row
		{
			qState * sa = & s[fP->pr[a].i];
			qState * sb = & s[fP->pr[a].j];
			bool isHoleA = sa->isHole;
			bool isHoleB = sb->isHole;
			if(isHoleA || isHoleB) continue;
			int pphhI = 0;
			for (int i = 0; i < size; i++) //  loop over column
			{

				qState * si = & s[fP->pr[i].i];
				qState * sj = & s[fP->pr[i].j];


				bool isHoleI = si->isHole;
				bool isHoleJ = sj->isHole;

				if (isHoleI && isHoleJ)
				{
					corr_en += 0.25 * fVpphh(pphhA, pphhI)
						* fTpphh(pphhA, pphhI);

					fTpphh(pphhA, pphhI) += factor
						* fHamil(pphhA, pphhI)
						/ (-FSum(&(fP->pr[a]), &(fP->pr[i])));
					pphhI++;
				}

			}
			if (pphhI) pphhA++;
		}

		cout << "\tcorr_en: " << corr_en;
		cout << " \tNSteps: " << NSteps;
		cout << "\r" << flush;

		correlationEn.push_back(corr_en);
		NSteps++;
//		corr_en_file << corr_en << endl;
	}
//	corr_en_file.close();
	
	return corr_en;
}

void CCMInf::Print()
{
	fHamil.print();
}
