#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#define M_PI       3.14159265358979323846
using namespace std;

int Divide_potential(double* ueff, double* u_att, double* u_rep, int M, int N_sig, int* n_0) {
	int i, nmin, nzero;
	double u_min = 10.0;
	for(i=0;i<M;i++){
		u_min>ueff[i] ? (u_min=ueff[i], nmin=i) : i;
	}
	for(i=0;i<nmin;i++){
		u_rep[i] = ueff[i] - u_min;
		u_att[i] = u_min;
	}
	for(i=nmin;i<M; i++) {
		u_rep[i] = 0;
		u_att[i] = ueff[i];
	}
	//Find sigma, where U_eff(\sigma) \approx 0
	nzero = 2*N_sig;
	for(i=2;i<2*N_sig;i++){
		if ((ueff[i-1]>=0)&&(ueff[i]<0))
			nzero = i;
	}
	*n_0 = nzero;
	return 0;
}

int calculate_sk(double* hr, double *sk, double rho, double dr, int M)
{
	double dk = (M_PI/M)/dr;
	double* k = new double [M];
	double* r = new double [M];
	int i, j;
	for (i = 0; i<M; i++) sk[i] = 0.0;
	for (i=0; i<M; i++) r[i] = dr * (i+0.5);
	for (i=0; i<M; i++) {
		k[i] = dk * (i+1.0);
		for (j=0; j<M; j++) {
			sk[i] = sk[i] + r[j] * sin(r[j]*k[i]) * hr[j] *dr / k[i]; 
		}
		sk[i] = sk[i] * 4.0 * M_PI * rho + 1.0; 
	}
	delete [] k;
	delete [] r;

	return 0;
}

int Baxtor(double* ueff, double* sk, double* prs, double* r, int M, int N_sig, double beta, double rho, double dr, double mix, double Pe, bool Lable_HI) 
{
	double* u_att = new double [M];
	double* u_rep = new double [M];
	int n_0 = 2*N_sig;
	int* p_n0 = &n_0;
	int i, j;
	double temp, err = 1.0;
	Divide_potential (ueff, u_att, u_rep, M, N_sig, p_n0);

	double* cr = new double [M];
	double* hr = new double [M];
	double* jr = new double [M];
	double* ir = new double [M];
	double* h1 = new double [M];
	double* dqr = new double [M];
	double* tp = new double [M];
	double* qr = new double [M];
	double* k = new double [M];
	ifstream dqr_ini;
	dqr_ini.open("Dqr.dat");
	for(i=0;i<M;i++) dqr[i]=0.0;
	for(i=0; i<((M<1000)?M:1000);i++){
		dqr_ini >> temp >> dqr[i];
	}
	for(i=0;i<M;i++){
		cr[i] = -beta * ueff[i];
		hr[i] = -1.0;
		jr[i] = 0.0;
		ir[i] = 0.0;
		h1[i] = 0.0;
	}
	if (n_0 < 1.0 * N_sig)	//appropriate cut off
		N_sig = n_0;
	do 
	{
		// h(r>sigma)
		for (i=N_sig;i<M;i++){
			tp[i] = (-dqr[i]+jr[i])/r[i];
			hr[i] = mix * hr[i] + (1.0-mix) * tp[i];
		}
		//c(r<sigma)
		for (i=0;i<N_sig;i++){
			tp[i] = (-dqr[i]+ir[i])/r[i];
			cr[i] = mix * cr[i] + (1.0-mix) * tp[i];
		}
		//q'(r<sigma)
		for (i=0;i<N_sig;i++) {
			tp[i] = -r[i]*hr[i]+jr[i];
			dqr[i] = mix * dqr[i] + (1.0-mix) * tp[i];
		}
		//q'(r>sigma)
		for (i=N_sig;i<M;i++) {
			tp[i] = -r[i]*cr[i]+ir[i];
			dqr[i] = mix * dqr[i] + (1.0-mix) * tp[i];
		}
		//q(r)
		qr[M-1] = 0.0;
		for (i=M-2;i>=0;i--) {
			qr[i] = qr[i+1] - dr * dqr[i];
		}
		//Ir
		for (i=0; i<M; i++) ir[i] = 0.0;
		for (i=0; i < M-1; i++) {
			for (j = i+1; j<M; j++) {
				ir[i] = ir[i] + dr * qr[j] *dqr[j-i];
			}
			ir[i] = 2.0 * M_PI * rho * ir[i]; 
		}
		//Jr
		for (i=0;i<M;i++) jr[i] = 0.0;
		for (i=0;i<M;i++) {
			for (j=0;j<M; j++) {
				jr[i] = jr[i] + dr*qr[j] * (r[i]-r[j]) * hr[abs(i-j)];
			}
			jr[i] = jr[i] * 2.0 * M_PI * rho;
		}
	//closure
		//c(r>sigma)
		for (i=N_sig; i<M; i++ ) {
			cr[i] = (exp(-beta*u_rep[i])-1.0) * (hr[i]+1.0) / exp(-beta * u_rep[i]) - beta * u_att[i];
		}
		//h(r<sig)
		for (i=0; i<N_sig; i++) {
			hr[i] = -1.0 + (beta*u_att[i]+cr[i]) * exp(-beta*u_rep[i]) / (exp(-beta*u_rep[i])-1.0); //SMSA
		}

		//err
		err = 0.0;
		for (i=0; i<M; i++) {
			err += (hr[i] - h1[i])*(hr[i] - h1[i]);
		}
		err *= dr;
		cout << "err: " << err <<endl;
		for(i=0; i<M; i++) h1[i] = hr[i]; 
	} while (err > 1.0e-8);
	calculate_sk(hr, sk, rho, dr, M);

	ofstream fout_gr;
	ofstream fout_sk;
	double dk= M_PI/(M*dr);
	for(i=0;i<M;i++) k[i] = dk*(i+1);
	char fn[30], fs[30] ;
	if (Lable_HI == 0){
		sprintf(fn, "data/gr_%5.3f_%5.1f.dat", rho, Pe) ;
		sprintf(fs, "data/sk_%5.3f_%5.1f.dat", rho, Pe) ;
		fout_gr.open(fn);
		fout_sk.open(fs);
	}
	else {
		sprintf(fn, "data/gr_%5.3f_%5.1f_HI.dat", rho, Pe) ;
		sprintf(fs, "data/sk_%5.3f_%5.1f_HI.dat", rho, Pe) ;
		fout_gr.open(fn);
		fout_sk.open(fs);
	}
	for (i=0;i<M;i++){
		fout_gr << r[i] <<"    "<< hr[i]+1.0 <<endl;
		fout_sk << k[i] <<"    "<< sk[i] <<endl;
	}
	fout_gr.close();
	fout_sk.close();

	//Pressure
	double pres=0 ;
	for (i=0; i<M-1; i++) {
		//pres += (ueff[i+1] - ueff[i]) * (hr[i] + 1.0) * pow(r[i], 3);
		pres += r[i] * r[i] * (exp(-beta * ueff[i]) - 1.0 );
	}
	//*prs = 1.0 - 2.0*M_PI*beta*rho * pres /3.0;
	*prs = 1.0 - 2.0*M_PI * rho * pres * dr;

	delete [] u_rep;
	delete [] u_att;
	delete [] cr;
	delete [] hr;
	delete [] jr;
	delete [] ir;
	delete [] h1;
	delete [] dqr;
	delete [] tp;
	delete [] qr;
	delete [] k;
	return 0;
}
