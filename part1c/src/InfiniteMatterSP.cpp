#include "InfiniteMatterSP.h"
#include "CCM.h"

InfiniteMatterSP::InfiniteMatterSP(int Nmax, int nShell, double rho) : 
		fNmax(Nmax), fNShell(nShell), fA(0), fNSP(0), fRho(rho)
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

inline double r(int x, int y, int z)
{
	return sqrt(x*x+y*y+z*z);
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
}

P_group_t *InfiniteMatterSP::P(int nX, int nY, int nZ)
{
	const int PM = 2*fNmax; // maximum P
	const int nP = 4*fNmax + 1; // maximum P
	const int n[3] = {nX+PM, nY+PM, nZ+PM};
//	cout << "n[0]: " << nX << "\tn[1]: " << nY; // DEBUG
//	cout << "\tn[2]: " << nZ << endl; // DEBUG
//	cout << "PM: " << PM << endl; // DEBUG
//	if(nX == -4) getchar(); // DEBUG
	const int index = (n[0]*nP + n[1])*nP + n[2];
	if(index >= nP*nP*nP) cout << "BANG!" << endl;
	return &fP[index];
}

void InfiniteMatterSP::ConstructPairs()
{
	vector<qState> &s = fStates;
	cout << fNSP << endl;
	for(int i = 0; i < fNSP; i++){
		for(int j = i + 1; j < fNSP; j++){
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
			if(p->np+1 > 10000) cout << "BANG!" << endl;
			p->pr[p->np++] = pr;
//			a->print(); b->print();
//			p->print(s); getchar(); // DEBUG
		} // end for over j
	} // end for over i
	const int n = pow(4*fNmax + 1, 3); // DEBUG
	for(int i = 0; i < n; i++){ // DEBUG
		for(int j = 0; j < n; j++){
			for(int k = 0; k < n; k++){
				cout << "i: " << i; // DEBUG
				cout << "\tj: " << j; // DEBUG
				cout << "\tk: " << k << endl; // DEBUG
				P(i - 2*fNmax, j - 2*fNmax, k - 2*fNmax)->print(fStates);
				getchar();
			} // end for over k
		} // end for over j
	} // DEBUG
} // end of member function ConstructPair

double P_group_t::CorrelationEnergy(){
	return CCM(this).SolveT();
}

double InfiniteMatterSP::CorrelationEnergy()
{
	const int n = pow(4*fNmax + 1, 3);
	
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			for(int k = 0; k < n; k++){
				cout << "i: " << i;
				cout << "\tj: " << j;
				cout << "\tk: " << k << endl;
				corr_en += P(i - 2*fNmax, j - 2*fNmax, k - 2*fNmax)->CorrelationEnergy();
				cout << "corr_en: " << corr_en << endl;
				getchar();
			} // end for over k
		} // end for over j
	} // end for over i
	
	return corr_en;
}

void InfiniteMatterSP::Print() const
{
	cout.precision(3);
	cout << "The quantum numbers are: " << endl;
	cout << "nx \t" << "ny \t" << "nz \t" <<  "isUp \t" << "r \t" << "nshell \t" << "isHole" << endl;
	
	for (auto t : fStates)
		cout << t.nx << " \t" << t.ny << " \t" << t.nz << " \t" << t.isUp << " \t" << t.r() << " \t" << t.nshell << " \t" << t.isHole << endl;
	
	cout << "The number of single particle states is: " << fStates.size() << endl;
	
	for(auto &t : RSHELL) cout << t*t << endl;
	cout << "Number of particles: " << fA << endl;
	cout << "Done!" << endl;
}

inline int delta(int i, int j)
{
	return i == j;
}

double InfiniteMatter::Minnesota(const pair_t * t, const pair_t * s)
{
	double L = pow(fNSP * fRho, 1 / 3.); 
	double fact = 2 * PI / L;

	double npRelx = rp(t.i, fStates);
	double npRely = rp(t.i, fStates);
	double npRelz = rp(t.i, fStates);
	double nqRelx = rp(t.j, fStates);
	double nqRely = rp(t.j, fStates);
	double nqRelz = rp(t.j, fStates);

	double nrRelx = rp(s.i, fStates);
	double nrRely = rp(s.i, fStates);
	double nrRelz = rp(s.i, fStates);
	double nsRelx = rp(s.j, fStates);
	double nsRely = rp(s.j, fStates);
	double nsRelz = rp(s.j, fStates); 

	double qx = fact * (npRelx - nqRelx + nrRelx + nsRelx);
	double qy = fact * (npRely - nqRely + nrRely + nsRely);
	double qz = fact * (npRelz - nqRelz + nrRelz + nsRelz);

	return (VR / pow(L, 3.) * pow(PI / KAPPAR, 2 / 3.) 
		* exp((qx*qx + qy*qy + qz*qz) / (4 * KAPPAR))
		//
		+ VS / pow(L, 3.) * pow(PI / KAPPAS, 2 / 3.) 
		* exp((qx*qx + qy*qy + qz*qz) / (4 * KAPPAS)))
		//
		* delta(knRelx + nqRelx, nrRelx + nsRelx)
		* delta(knRely + nqRely, nrRely + nsRely)
		* delta(knRelz + nqRelz, nrRelz + nsRelz)
		* (fStates[t.i].isUp && fStates[s.i].isUp 
			* fStates[t.j].isUp && fStates[s.j].isUp
			- fStates[t.i].isUp && fStates[s.j].isUp
			* fStates[t.j].isUp && fStates[s.i].isUp);
}
