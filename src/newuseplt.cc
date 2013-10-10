//newuseplt.cc
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "newplt.h"
#include <R.h>
#include <Rmath.h>
#include <iomanip>

using namespace std;

extern "C" {void plt_main(double *p, double *q, int *n, double *tt, double *cause, double *group, 
	double *tjp, int *njp, double *ny1, double *f11, double *sdf11, double *ny2, double *f21, double *sdf21,
	double *dif, double *sddif, double *pvdif, double *rr, double *sdrr, double *pvrr, double *Or, double *sdor, double *pvor, 
	int *size, double *ave, double *avese, double *aveu95, double *avel95, double *avepval, int *nbound, double *cbcut, double *wt,
	int *conf_bd, int *n_sim, int *causeS);}

void plt_main(double *p, double *q, int *n, double *tt, double *cause, double *group, 
	double *tjp, int *njp, double *ny1, double *f11, double *sdf11, double *ny2, double *f21, double *sdf21,
	double *dif, double *sddif, double *pvdif, double *rr, double *sdrr, double *pvrr, double *Or, double *sdor, double *pvor,
	int *size, double *ave, double *avese, double *aveu95, double *avel95, double *avepval, int *nbound, double *cbcut, double *wt,
	int *conf_bd, int *n_sim, int *causeS){
	
	PLT plt(p, q, n, tt, cause, group, size, njp, tjp);
	
	//cout<<"Can you see me?"<<endl;
	
	int N=*n;   //define global and private
	int NOGP=2; //define global and private
	int NJP;    //define global and private
	NJP = *njp;
	int ntau1=0, ntau2=0;

	double f1diff[N+1], sgf1diff[N+1], pvf1diff[N+1], rr11[N+1], sgrr11[N+1], pvrr11[N+1], rr12[N+1], sgrr12[N+1], pvrr12[N+1];
	// //rr_ij, i=cause, j=1, rr; j=2, or.
	double **ny = new double *[NOGP];
	for(int i=0;i<NOGP;i++) ny[i] = new double[N+1];
	double **f1 = new double *[NOGP];
	for(int i=0;i<NOGP;i++) f1[i] = new double[N+1];
	double **sgf1 = new double *[NOGP];
	for(int i=0;i<NOGP;i++) sgf1[i] = new double[N+1];

	double ***xf1=new double**[NOGP];
	for(int i=0;i<NOGP;i++){
		xf1[i] = new double*[N+1];
	}
	for(int i=0;i<NOGP;i++){
		for(int j=0; j<N; j++){
			xf1[i][j] = new double[N+1];
		}
	}	

	void bubble_sort(double arrin[], int n);
	plt.CIF_est(ny, f1, sgf1, xf1);
	
	for(int j=0; j<=NJP; j++){                //quantities sent back to R
		ny1[j] = ny[0][j];
		f11[j] = f1[0][j];
		sdf11[j] = sqrt(sgf1[0][j]);
		ny2[j] = ny[1][j];
		f21[j] = f1[1][j];
		sdf21[j] = sqrt(sgf1[1][j]);
	}

	for(int j=1;j<=NJP;j++){
		if((f1[0][j]>0)&&(f1[1][j]>0)){
			ntau1=j;
			for(int k=NJP; k>=j ; k--){
				if((f1[0][k]<1)&&(f1[1][k]<1)){
					ntau2=k;
					break;
				}
			}
			if(ntau2>0) break;
		}
	}
	nbound[0] = ntau1+1; nbound[1] =ntau2+1;  //+1 consistent with starting t from 0, quantities return to R

	plt.G_est(ntau1, ntau2, f1, sgf1,
		      f1diff, sgf1diff, pvf1diff, rr11, sgrr11, pvrr11, rr12, sgrr12, pvrr12);
	

	for(int j=0;j<=NJP;j++){                    //quantities sent back to R
		dif[j] = f1diff[j];
		sddif[j] = sqrt(sgf1diff[j]);
		pvdif[j] = pvf1diff[j];
		rr[j] = rr11[j];
		sdrr[j] = sqrt(sgrr11[j]);
		pvrr[j] = pvrr11[j];
		Or[j] = rr12[j];
		sdor[j] = sqrt(sgrr12[j]);
		pvor[j] = pvrr12[j];
	}
	
	plt.WT(wt, ntau1, ntau2);
	plt.Ave_est(ntau1, ntau2, f1, xf1, rr11, rr12, wt, ave, avese, aveu95, avel95, avepval);
	if(*conf_bd==1){
		plt.CB(ntau1, ntau2, f1, sgf1, xf1, sgrr11, sgrr12, cbcut, n_sim);
	}

	//release
	for(int i=0;i<NOGP;i++) delete[] ny[i];  
	delete[] ny;
	for(int i=0;i<NOGP;i++) delete[] f1[i];
	delete[] f1;
	for(int i=0;i<NOGP;i++) delete[] sgf1[i];
	delete[] sgf1;
	
	for(int i=0; i<NOGP; i++){
		for(int j=0; j<N; j++){
			delete[] xf1[i][j];
		}
	}
	for(int i=0; i<NOGP; i++){
		delete[] xf1[i];
	}
	delete[] xf1; 
}
