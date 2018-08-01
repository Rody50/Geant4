#include "CCMInf.h"
#include <fstream>
#include <string>

using std::ofstream;
using std::to_string;
using std::flush;

CCMInf::CCMInf(P_group_t *p, InfiniteMatterSP *infSP) : fN(p->N()),
	fA(p->A()), fNu(p->N() - p->A()), fP(p),
	fHamil(p->N() - p->A(), p->A()),
	fVpppp(p->N() - p->A(), p->N() - p->A()),	
	fVpphh(p->N() - p->A(), p->A()), 
	fVhhpp(p->A(), p->N() - p->A()), 
	fVhhhh(p->A(), p->A()),
	fTpphh(p->N() - p->A(), p->A()),
	fInterm(p->N() - p->A(), p->A()),
	fInfSP(infSP)
{
	FillMatrix::FillV(fVpppp, fVpphh, fVhhhh,
		p, fInfSP);
	FillMatrix::FillT(fTpphh, p, this);
	//fVpppp.print();
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

	int na[3] = {a->nx, a->nx, a->nz};
	int nb[3] = {b->nx, b->nx, b->nz};
	int ni[3] = {i->nx, i->nx, i->nz};
	int nj[3] = {j->nx, j->nx, j->nz};
	const double u = fInfSP->GetfEnUn();

	const int np = fP->np;
	double va = fInfSP->Minnesota(tpp, tpp), vb = va;
	double vi = fInfSP->Minnesota(shh, shh), vj = vi;

	double fa = r2(na) * u + va;
	double fb = r2(nb) * u + vb;
	double fi = r2(ni) * u + vi;
	double fj = r2(nj) * u + vj;

	return fa + fb - fi - fj;
}

void CCMInf::Interm()
{
	vector<qState> &s = fInfSP->GetfStates();
	int size = fP->np;
	int pphhA = 0;

	for (int a = 0; a < size; a++) // loop over row
	{
		int pphhI = 0;
		for (int i = 0; i < size; i++) //  loop over column
		{
			qState * sa = & s[fP->pr[a].i];
			qState * sb = & s[fP->pr[a].j];
			qState * si = & s[fP->pr[i].i];
			qState * sj = & s[fP->pr[i].j];

			bool isHoleA = sa->isHole;
			bool isHoleB = sb->isHole;
			bool isHoleI = si->isHole;
			bool isHoleJ = sj->isHole;

			if (!isHoleA && !isHoleB && isHoleI && isHoleJ)
			{	
				fInterm(pphhA + pphhI * size) = fTpphh(pphhA + pphhI * size)
					* this->FSum(&(fP->pr[a]), &(fP->pr[i]));
				pphhI++;
			}

		}
		if (pphhI) pphhA++;
	}
}

void CCMInf::ComputeH()
{
	Interm();

	fHamil = fVpphh + fInterm + fVpppp * fTpphh + fTpphh * fVhhhh;
}

double CCMInf::SolveT()
{
	double factor = 1.;

	double corr_en_pre = 1.;
	double corr_en = 0.;
	vector<double> correlationEn;

	int NStepsMax = 1E5;
	int NSteps = 0;

	ofstream corr_en_file(("results/corr_en-steps" + to_string(fN) 
		+ "_" + to_string(fA) + ".txt").c_str());

	vector<qState> &s = fInfSP->GetfStates();

	while (fabs(corr_en - corr_en_pre) > 1E-5 || NSteps > NStepsMax)
	{	
		ComputeH();
		
		corr_en_pre = corr_en;
		corr_en = 0.;

		int size = fP->np;
		int pphhA = 0;

		for (int a = 0; a < size; a++) // loop over row
		{
			int pphhI = 0;
			for (int i = 0; i < size; i++) //  loop over column
			{
				qState * sa = & s[fP->pr[a].i];
				qState * sb = & s[fP->pr[a].j];
				qState * si = & s[fP->pr[i].i];
				qState * sj = & s[fP->pr[i].j];

				bool isHoleA = sa->isHole;
				bool isHoleB = sb->isHole;
				bool isHoleI = si->isHole;
				bool isHoleJ = sj->isHole;

				if (!isHoleA && !isHoleB && isHoleI && isHoleJ)
				{	
					corr_en += 0.25 * fVpphh(pphhA + pphhI * size)
						* fTpphh(pphhA + pphhI * size);

					fTpphh(pphhA + pphhI * size) += factor
						* fHamil(pphhA + pphhI * size)
						/ this->FSum(&(fP->pr[a]), &(fP->pr[i]));
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
		corr_en_file << corr_en << endl;
	}
	corr_en_file.close();

	cout << "The converged correlation energy is: " << corr_en << endl;
	
	return corr_en;
}

void CCMInf::Print()
{
	fHamil.print();
}
