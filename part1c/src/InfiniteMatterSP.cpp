#include "InfiniteMatterSP.h"
#include "CCMInf.h"

InfiniteMatterSP::InfiniteMatterSP(int Nmax, int nShell, double rho) : 
		fNmax(Nmax), fNShell(nShell), fA(0), fNSP(0), fRho(rho), fL(0),
		fEnUn(0)
{
	if(nShell > 20){
		cout << "nShell too large (has to be < 20): " << nShell << endl;
		exit(1);
	}
	const int nP = 4*Nmax + 1;
	fP = new P_group_t[nP*nP*nP]; // -2fNmax ~ 2fNmax
}

InfiniteMatterSP::~InfiniteMatterSP()
{
	delete [] fP;
}

inline int r2(int x, int y, int z)
{
	return x*x+y*y+z*z;
}

inline double r(int x, int y, int z)
{
	return sqrt(r2(x,y,z));
}

inline int delta(int i, int j)
{
	return i == j;
}

static const double RSHELL[20] = {
	r(0,0,0),r(1,0,0),r(1,1,0),r(1,1,1),r(2,0,0),
	r(2,1,0),r(2,1,1),r(2,2,0),r(2,2,1),r(3,1,0),
	r(3,1,1),r(2,2,2),r(3,2,0),r(3,2,1),r(4,0,0),
	r(3,2,2),r(3,3,0),r(3,3,1),r(4,2,0),r(4,2,1)
};

void InfiniteMatterSP::GenerateSP()
{
	qState tmp;
	for (int nx = -fNmax; nx <= fNmax; nx++){
		tmp.nx = nx;
		for (int ny = -fNmax; ny <= fNmax; ny++){
			tmp.ny = ny;
			for (int nz = -fNmax; nz <= fNmax; nz++){
				tmp.nz = nz;
				// tell shell id
				for(int i = 0; i < 20; i++){
					if(tmp.r() == RSHELL[i]){
						tmp.nshell = i;
						tmp.isHole = i <= fNShell;
						break;
					} // end if
				} // end for over i
				tmp.isUp = 1;
				fStates.push_back(tmp);
				tmp.isUp = 0;
				fStates.push_back(tmp);
				if(tmp.isHole) fA += 2;
			} // end for over nz
		} // end for over ny
	} // end for over nx
	fNSP = fStates.size();

	fL = pow(fA / fRho, 1 / 3.);
	fEnUn = 2 * PI * PI * hbarc * hbarc / (nmass * fL * fL); 
//	cout << "fL: " << fL << endl; getchar();

}

P_group_t *InfiniteMatterSP::P(int nX, int nY, int nZ)
{
	const int PM = 2*fNmax; // maximum P
	const int nP = 4*fNmax + 1; // maximum P
	const int n[3] = {nX+PM, nY+PM, nZ+PM};
	//cout << "n[0]: " << nX << "\tn[1]: " << nY; // DEBUG
	//cout << "\tn[2]: " << nZ << endl; // DEBUG
	//cout << "PM: " << PM << endl; // DEBUG
	//if(nX == -4) getchar(); // DEBUG
	const int index = (n[0]*nP + n[1])*nP + n[2];
	if(index >= nP*nP*nP){ cout << "BANG!1" << endl;
		cout << "n0: " << n[0] << "\tn1: " << n[1] << "\tn2: " << n[2] << endl;
	}
	return &fP[index];
}

void InfiniteMatterSP::ConstructPairs()
{
	vector<qState> &s = fStates;
	cout << fNSP << endl;
	for(int i = 0; i < fNSP; i++){
		for(int j = 0; j < fNSP; j++){
			if(i == j) continue;
			qState *a = &s[i], *b = &s[j];
			const int n[3] = {a->nx+b->nx, a->ny+b->ny, a->nz+b->nz};
			P_group_t *p = P(n[0], n[1], n[2]);
			p->P[0] = n[0]; p->P[1] = n[1]; p->P[2] = n[2];
			pair_t pr;
			pr.i = i; pr.j = j;

			if(a->isHole && b->isHole) p->nhh++;
			if(!a->isHole && b->isHole) p->nph++;
			if(a->isHole && !b->isHole) p->nph++;
			if(!a->isHole && !b->isHole) p->npp++;
			if(p->np+1 > 10000) cout << "BANG!2" << endl;
			p->pr[p->np++] = pr;
//			a->print(); b->print();
//			p->print(s); getchar(); // DEBUG
		} // end for over j
	} // end for over i
	const int n = 4*fNmax + 1; // DEBUG
	int nn = 0;
	if(0)
	for(int i = 0; i < n; i++){ // DEBUG
		for(int j = 0; j < n; j++){
			for(int k = 0; k < n; k++){
				cout << "i: " << i; // DEBUG
				cout << "\tj: " << j; // DEBUG
				cout << "\tk: " << k << endl; // DEBUG
				P_group_t *p = P(i - 2*fNmax, j - 2*fNmax, k - 2*fNmax);
				p->print(fStates);
				getchar();
				if(p->nhh != 0)
					cout << nn++ << endl;
			} // end for over k
		} // end for over j
	} // DEBUG
} // end of member function ConstructPair

double P_group_t::CorrelationEnergy(InfiniteMatterSP *InfSP){
	return CCMInf(this, InfSP).SolveT();
}

double InfiniteMatterSP::CorrelationEnergy()
{
	const int n = 4*fNmax + 1;
	double corr_en = 0.;
	
	//P(-2, 0, 0)->print(fStates);
	//return P(-2, 0, 0)->CorrelationEnergy(this);
	
	
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < n; k++){
//				cout << "i: " << i;
//				cout << "\tj: " << j;
//				cout << "\tk: " << k << endl;
//				cout << "\tTotal momentum group: ";
				P_group_t *p = P(i - 2*fNmax, j - 2*fNmax, k - 2*fNmax);
				//p->print(fStates);
				if(p->nhh == 0) continue;
					corr_en += P(i - 2*fNmax, j - 2*fNmax, k - 2*fNmax)->CorrelationEnergy(this);
//				cout << "corr_en: " << corr_en << endl;
				//getchar();
			} // end for over k
		} // end for over j
	} // end for over i
	
	return corr_en / fA;
}

void InfiniteMatterSP::Print() const
{
//	cout.precision(3);
	cout << "The quantum numbers are: " << endl;
	cout << "nx \t" << "ny \t" << "nz \t" <<  "isUp \t" << "r \t" << "nshell \t" << "isHole" << endl;
	
	for (auto t : fStates)
		cout << t.nx << " \t" << t.ny << " \t" << t.nz << " \t" << t.isUp << " \t" << t.r() << " \t" << t.nshell << " \t" << t.isHole << endl;
	
	cout << "The number of single particle states is: " << fStates.size() << endl;
	
	for(auto &t : RSHELL) cout << t*t << endl;
	cout << "Number of particles: " << fA << endl;
	cout << "Done!" << endl;
}

double InfiniteMatterSP::Minnesota(const qState &ti, const qState &tj, 
const qState &si, const qState &sj)
{
	short spin[2][2] = { // [pair][qState]
		{ti.spin(), tj.spin()},
		{si.spin(), sj.spin()}
	};
	if(spin[0][0]*spin[0][1] > 0) return 0.; // spin triplet
	if(spin[1][0]*spin[1][1] > 0) return 0.; // spin triplet

	// anti-symmetrization
	short sign = 1;
	if(spin[0][0] == +1 && spin[0][1] == -1 && 
		 spin[1][0] == +1 && spin[1][1] == -1) sign = +1.; // <+-|V|+->

	if(spin[0][0] == +1 && spin[0][1] == -1 && 
		 spin[1][0] == -1 && spin[1][1] == +1) sign = -1.; // <+-|V|-+>

	if(spin[0][0] == -1 && spin[0][1] == +1 && 
		 spin[1][0] == +1 && spin[1][1] == -1) sign = -1.; // <-+|V|+->

	if(spin[0][0] == -1 && spin[0][1] == +1 && 
		 spin[1][0] == -1 && spin[1][1] == +1) sign = +1.; // <-+|V|-+>

	int qij[3], qji[3];


	if(0) if(ti.nx+tj.nx != si.nx+sj.nx)
	{
		cout << "0" << endl;
		cout << "ti.nx: " << ti.nx << "\ttj.nx: " << tj.nx << endl;
		cout << "si.nx: " << si.nx << "\tsj.nx: " << sj.nx << endl;
		cout << "ti.ny: " << ti.ny << "\ttj.ny: " << tj.ny << endl;
		cout << "si.ny: " << si.ny << "\tsj.ny: " << sj.ny << endl;
		cout << "ti.nz: " << ti.nz << "\ttj.nz: " << tj.nz << endl;
		cout << "si.nz: " << si.nz << "\tsj.nz: " << sj.nz << endl;
		getchar();
		return 0;
	}
	if(0) if(ti.ny+tj.ny != si.ny+sj.ny)
	{
		cout << "1" << endl;
		cout << "ti.nx: " << ti.nx << "\ttj.nx: " << tj.nx << endl;
		cout << "si.nx: " << si.nx << "\tsj.nx: " << sj.nx << endl;
		cout << "ti.ny: " << ti.ny << "\ttj.ny: " << tj.ny << endl;
		cout << "si.ny: " << si.ny << "\tsj.ny: " << sj.ny << endl;
		cout << "ti.nz: " << ti.nz << "\ttj.nz: " << tj.nz << endl;
		cout << "si.nz: " << si.nz << "\tsj.nz: " << sj.nz << endl;
		getchar();
		return 0;
	}
	
	if(0) if(ti.nz+tj.nz != si.nz+sj.nz)
	{
		cout << "2" << endl;
		cout << "ti.nx: " << ti.nx << "\ttj.nx: " << tj.nx << endl;
		cout << "si.nx: " << si.nx << "\tsj.nx: " << sj.nx << endl;
		cout << "ti.ny: " << ti.ny << "\ttj.ny: " << tj.ny << endl;
		cout << "si.ny: " << si.ny << "\tsj.ny: " << sj.ny << endl;
		cout << "ti.nz: " << ti.nz << "\ttj.nz: " << tj.nz << endl;
		cout << "si.nz: " << si.nz << "\tsj.nz: " << sj.nz << endl;
		getchar();
		return 0;
	}


	qij[0] = ti.nx-tj.nx-si.nx+sj.nx;
	qij[1] = ti.ny-tj.ny-si.ny+sj.ny;
	qij[2] = ti.nz-tj.nz-si.nz+sj.nz;

	qji[0] = tj.nx-ti.nx-si.nx+sj.nx;
	qji[1] = tj.ny-ti.ny-si.ny+sj.ny;
	qji[2] = tj.nz-ti.nz-si.nz+sj.nz;

	int q2ij = r2(qij[0], qij[1], qij[2]);
	int q2ji = r2(qji[0], qji[1], qji[2]);

	static double fact = PI / fL * PI / fL;
	const double v = 0.5 *
		( // ij term
		(VR / pow(fL, 3.) * pow(PI / KAPPAR, 3. / 2)
		* exp(-q2ij * fact / (4 * KAPPAR))
		+ VS / pow(fL, 3.) * pow(PI / KAPPAS, 3. / 2)
		* exp(-q2ij * fact / (4 * KAPPAS)))
		
		 + 
		 
		// ji term
		(VR / pow(fL, 3.) * pow(PI / KAPPAR, 3. / 2)
		* exp(-q2ji * fact / (4 * KAPPAR))
		+ VS / pow(fL, 3.) * pow(PI / KAPPAS, 3. / 2)
		* exp(-q2ji * fact / (4 * KAPPAS)))
		)
		*(delta(spin[0][0], spin[1][0]) * delta(spin[0][1], spin[1][1]) 
		- delta(spin[0][0], spin[1][1]) * delta(spin[0][1], spin[1][0]));
		
		return v; // * sign;
}
double InfiniteMatterSP::Minnesota(const pair_t * t, const pair_t * s)
{
//	t->print(fStates); s->print(fStates);
	return Minnesota(fStates[t->i], fStates[t->j], fStates[s->i], fStates[s->j]);
}

