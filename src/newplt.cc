//newplt.cc
#include <iostream>
#include <fstream>
#include <cmath>
#include "newplt.h"
#include <R.h>
#include <rmath.h>
#include <iomanip>
#include <new>

using namespace std;

void bubble_sort(double arrin[], int n){
	double temp;
	for(int i = 0; i < n-1; i++){
		for(int j = 0; j < n-i-1; j++){
			if(arrin[j]>arrin[j+1]){
				temp= arrin[j];
				arrin[j] = arrin[j+1];
				arrin[j+1] = temp;
			}
		}
	}
	return;
}

PLT::PLT(double *p, double *q, int *n, double *tt, double *cause, double *group, int *size, int *njp, double *tjp){
	N = *n;
	P = *p;
	Q = *q;
	TT = Array1D<double> (N,0.);     
	Cause = Array1D<double>(N,0.);  
	IGP = Array1D<double> (N,0.);
	TJP = Array1D<double> (N+1,0.);
	NSIZE = Array1D<int> (2,0);
	T = Array2D<double> (2,N,0.);  
	IC = Array3D<double>(2,N,3,0.);
	IDX = Array1D<int>(N+1,0);
	TTJP = Array1D<double> (N+1,0.);

	NOGP=2;
	Z95=1.96;
	//NRUN=500;
	
	for(int i=0; i<N; i++){
		TT[i] = tt[i];
		Cause[i] = cause[i];
	}

	int ct1=0, ct2=0;
	
	for(int i=0; i<N; i++){
		if(group[i]==1){
			IGP[i]=1;
			ct1 += 1;
		}
		else{
			IGP[i]=2;
			ct2 += 1;
		}

		int ig=0;
		int ns=0;
		ig = IGP[i]-1;   
		if(ig==0) ns=ct1-1;
		else ns=ct2-1;

		T[ig][ns]=tt[i];

		IC[ig][ns][1]=0;
		if(cause[i]==1) IC[ig][ns][1]=1;
		IC[ig][ns][2]=0;
		if(cause[i]==2) IC[ig][ns][2]=1;
		IC[ig][ns][0] = IC[ig][ns][1]+IC[ig][ns][2];
	}
	
	int jp1=0;
	for(int i=0; i<N; i++){                     //TJP 0:NJP
		if((tt[i]>TJP[jp1])&&(cause[i]==1)){     //should be cause[i]>=0 if censoring matters, like FG model; cause[i]>=1 if calculate 2 causes together
			jp1 += 1;
			TJP[jp1] = tt[i];                    //private member
			tjp[jp1] = tt[i];                    //make a copy return to R
		}
	}
	NJP = jp1;                                   //private member
	*njp = jp1;                                   //make a copy return to R
	NSIZE[0] = ct1; NSIZE[1] = ct2;             //private member
	size[0] =ct1; size[1] =ct2;                 //make a copy return to R
	
	int jp2=0;

	for(int i=0; i<N; i++){                     
		if((tt[i]>TTJP[jp2])&&(cause[i]>=1)){     //should be cause[i]>=0 if censoring matters, like FG model; cause[i]>=1 if calculate 2 causes together
			jp2 += 1;
			TTJP[jp2] = tt[i];                    //private member
		}
	}
	TNJP = jp2;
	
	int ct=0;
	for(int i=0; i<=TNJP; i++){
		for(int j=ct; j<=NJP; j++){
			if(TJP[j]==TTJP[i]){
				IDX[ct]=i;
				ct++;
				break;
			}
		}
	}
}
void PLT::CIF_est(double **ny, double **f1, double **sgf1, double ***xf1){
	Array2D<double> surv(NOGP,N+1,0.);
	Array3D<double> nd (NOGP,2,N+1,0.);
	Array3D<double> dlambda(NOGP,2,N+1,0.);
	
	double ***dmhat0=new double**[NOGP];
	for(int i=0;i<NOGP;i++){
		dmhat0[i] = new double*[N+1];
	}
	for(int i=0;i<NOGP;i++){
		for(int j=0; j<N; j++){
			dmhat0[i][j] = new double[N+1];
		}
	}

	double ***dmhat1=new double**[NOGP];
	for(int i=0;i<NOGP;i++){
		dmhat1[i] = new double*[N+1];
	}
	for(int i=0;i<NOGP;i++){
		for(int j=0; j<N; j++){
			dmhat1[i][j] = new double[N+1];
		}
	}	
	double ***xs=new double**[NOGP];
	for(int i=0;i<NOGP;i++){
		xs[i] = new double*[N+1];
	}
	for(int i=0;i<NOGP;i++){
		for(int j=0; j<N; j++){
			xs[i][j] = new double[N+1];
		}
	}

	for(int kg=0; kg<NOGP; kg++){          //initialization
		for(int i=0;i<NSIZE[kg];i++){
			for(int j=0;j<=TNJP;j++){
				xs[kg][i][j] = 0.;
				xf1[kg][i][j] = 0.;
				dmhat0[kg][i][j] = 0.;
				dmhat1[kg][i][j] = 0.;
			}
		}
		for(int j=0;j<N;j++){
			f1[kg][j]=0.;
			sgf1[kg][j]=0.;
			ny[kg][j]=0.;
		}
	}
	
	for(int kg=0; kg<NOGP; kg++){  
		surv[kg][0]=1.0;
		ny[kg][0]=NSIZE[kg];
		for(int j=1; j<=TNJP;j++){
			for(int i=0;i<NSIZE[kg];i++){ 		
				if(T[kg][i]==TTJP[j]){ 				
					if(IC[kg][i][0] == 1.) nd[kg][0][j]=nd[kg][0][j]+1.;
					if(IC[kg][i][1] == 1.) nd[kg][1][j]=nd[kg][1][j]+1.;
				}
				if(T[kg][i] >= TTJP[j]) ny[kg][j]=ny[kg][j]+1.0;
			}
			if(ny[kg][j]>0.0){
				dlambda[kg][0][j] = nd[kg][0][j]/ny[kg][j];
				dlambda[kg][1][j] = nd[kg][1][j]/ny[kg][j];
			}
			else{
				dlambda[kg][0][j]=0.;
				dlambda[kg][1][j]=0.;
			}
			surv[kg][j]=surv[kg][j-1]*(1.0-dlambda[kg][0][j]);            //j starts from 1, diff from others
			f1[kg][j]=f1[kg][j-1]+surv[kg][j-1]*dlambda[kg][1][j]; 

			for(int i=0;i<NSIZE[kg];i++){
				if(T[kg][i] >= TTJP[j]){
					if((T[kg][i]==TTJP[j])&(IC[kg][i][0]==1)) dmhat0[kg][i][j]=1.0-dlambda[kg][0][j];
					else dmhat0[kg][i][j]=0.0-dlambda[kg][0][j];
					if((T[kg][i]==TTJP[j])&(IC[kg][i][1]==1)) dmhat1[kg][i][j]=1.0-dlambda[kg][1][j];
					else dmhat1[kg][i][j]=0.0-dlambda[kg][1][j];
				}
				if(ny[kg][j] > 0.0){
					xs[kg][i][j]=xs[kg][i][j-1] + dmhat0[kg][i][j]/ny[kg][j];
					xf1[kg][i][j]=xf1[kg][i][j-1] + surv[kg][j-1]*((dmhat1[kg][i][j]/ny[kg][j]) - (xs[kg][i][j-1]*dlambda[kg][1][j]));
				}
				else{
					xs[kg][i][j]=xs[kg][i][j-1];
					xf1[kg][i][j]=xf1[kg][i][j-1];
				}
				sgf1[kg][j]=sgf1[kg][j]+pow(xf1[kg][i][j],2.0);
			}	
		}
	}
		
	for(int kg=0; kg<NOGP; kg++){ 
		for(int j=1; j<=TNJP;j++){
			if(j<=NJP){
				ny[kg][j] = ny[kg][IDX[j]];
				f1[kg][j] = f1[kg][IDX[j]];
				sgf1[kg][j] = sgf1[kg][IDX[j]];
				for(int i=0;i<NSIZE[kg];i++){
					xf1[kg][i][j]=xf1[kg][i][IDX[j]];
				}
			}
			else{
				ny[kg][j] = 0.;
				f1[kg][j] = 0.;
				sgf1[kg][j] = 0.;
				for(int i=0;i<NSIZE[kg];i++){
					xf1[kg][i][j] = 0.;
				}
			}
		}
	}
			

	for(int i=0; i<NOGP; i++){                        //release 
		for(int j=0; j<N; j++){
			delete[] dmhat0[i][j];
		}
	}
	for(int i=0; i<NOGP; i++){
		delete[] dmhat0[i];
	}
	delete[] dmhat0; 
		
	for(int i=0; i<NOGP; i++){
		for(int j=0; j<N; j++){
			delete[] dmhat1[i][j];
		}
	}
	for(int i=0; i<NOGP; i++){
		delete[] dmhat1[i];
	}
	delete[] dmhat1; 	
		
	for(int i=0; i<NOGP; i++){
		for(int j=0; j<N; j++){
			delete[] xs[i][j];
		}
	}
	for(int i=0; i<NOGP; i++){
		delete[] xs[i];
	}
	delete[] xs;
	
	return;
}

void PLT::G_est(int ntau1, int ntau2, double **f1, double **sgf1, double f1diff[], double sgf1diff[], double pvf1diff[], 
			double rr11[], double sgrr11[], double pvrr11[], double rr12[], double sgrr12[], double pvrr12[]){
	double stat;
	double tmp1,tmp2;
	double u1,v1;
	double sglog11=0.;
	double sglog12=0.;
	double cc=10e-6;
	
	for(int j=0; j<=NJP; j++){
		u1 = f1[0][j];   
		v1 = f1[1][j];
		f1diff[j] = f1[0][j]-f1[1][j];
		sgf1diff[j] = sgf1[0][j]+sgf1[1][j];
		if(sgf1diff[j] > cc){
			stat = pow(f1diff[j],2)/sgf1diff[j];
			pvf1diff[j] = pchisq(stat,1.,0,0); 
		}
		
		if((f1[0][j])>cc && (f1[1][j]>cc)){
			rr11[j]=f1[0][j]/f1[1][j];
			tmp1=1.0/v1;
			tmp2=0.0-u1/pow(v1,2.0);
			sgrr11[j]=pow(tmp1,2.0)*sgf1[0][j]+pow(tmp2,2.0)*sgf1[1][j];
			sglog11=pow((1.0/rr11[j]),2)*sgrr11[j];
		}
		if(sglog11>cc){
			stat=pow(log(rr11[j]),2)/sglog11;
			pvrr11[j]=pchisq(stat,1.,0,0);
		}

		if((f1[0][j]>cc)&&(f1[1][j]>cc)&&(f1[0][j]<1-cc)&&(f1[1][j]<1-cc)){
			tmp1=u1/(1.0-u1);
			tmp2=v1/(1.0-v1);
			rr12[j]=tmp1/tmp2;
			tmp1=(1.0-v1)/(v1*pow(1.0-u1,2.0));
			tmp2=0.0-(u1/(pow(v1,2.0)*(1.0-u1)));
			sgrr12[j]=pow(tmp1,2.0)*sgf1[0][j] + pow(tmp2,2.0)*sgf1[1][j];
			sglog12=pow(1.0/rr12[j],2.0)*sgrr12[j];
		}
		if(sglog12>cc){
			stat=pow(log(rr12[j]),2.0)/sglog12;
			pvrr12[j]=pchisq(stat,1.,0,0);
		}
	}
}

void PLT::WT(double *wt, int ntau1, int ntau2){
	double wab;
	Array1D<double> surv(N+1,0.);
	Array1D<double> w(N+1,0.);
	Array1D<double> ny(N+1,0.);
	Array1D<double> f(N+1,0.);
	Array2D<double> nd(2,N+1,0.);
	Array2D<double> dlambda(2,N+1,0.);

	surv[0]=1.0;
	ny[0]=N;
	wab = 0.0;

	for(int j=1;j<=NJP;j++){
		for(int i=0;i<N;i++){
			if(TT[i]==TJP[j]){
				if(Cause[i] > 0) nd[0][j]=nd[0][j]+1;
				if(Cause[i] == 1) nd[1][j]=nd[1][j]+1;
			}
			if(TT[i] >= TJP[j]) ny[j]=ny[j]+1;
		}
		if(ny[j]>0.001){
			dlambda[0][j] = nd[0][j]/ny[j];
			dlambda[1][j] = nd[1][j]/ny[j];
		}
		else{
			dlambda[0][j]=0.;
			dlambda[1][j]=0.;
		}
		surv[j]=surv[j-1]*(1.0-dlambda[0][j]);                
		f[j]=f[j-1]+surv[j-1]*dlambda[1][j];         
	}

	for(int j=ntau1;j<ntau2;j++){
		w[j]=pow(1-f[j]/f[ntau2],P)*pow(f[j]/f[ntau2],Q);
		wab += w[j]*(TJP[j+1]-TJP[j]);
	}
	for(int j=ntau1; j<ntau2; j++){
		wt[j] = w[j]/wab;
	}
}

void PLT::Ave_est(int ntau1, int ntau2, double **f1, double ***xf1, double rr11[], double rr12[], double *wt, 
				 double *ave, double *avese, double *aveu95, double *avel95, double *avepval){
	double avef1d, sgavef1d, pv1d;
	double ave1dcil, ave1dciu;
	Array1D<double> avef1(NOGP,0.);
	Array2D<double> avexf1(NOGP,N,0.);
	Array2D<double> avexrr11(NOGP,N,0.);
	Array2D<double> avexrr12(NOGP,N,0.);
	Array1D<double> sgavef1(NOGP,0.);

	double stat;
	double tmp1,tmp2;
	double u1,v1;
	double sglog11=0.;
	double sglog12=0.;
	double averr11=0., averr12=0.;
	
	double sgaverr11=0.; 
	double sgaverr12=0.;

	for(int j=0; j<=NJP; j++){
		u1=f1[0][j];
		v1=f1[1][j];
		if((ntau1<=j)&&(j<ntau2)){               
			for(int kg=0;kg<NOGP;kg++){
				avef1[kg]=avef1[kg] + f1[kg][j]*(TJP[j+1]-TJP[j])*wt[j];
				for(int i=0;i<NSIZE[kg];i++){
					avexf1[kg][i]=avexf1[kg][i] + xf1[kg][i][j]*(TJP[j+1]-TJP[j])*wt[j];
				}
			}
			if((rr11[j]>0)&&(rr12[j]>0)){
				averr11=averr11 + rr11[j]*(TJP[j]-TJP[j-1])*wt[j];
				averr12=averr12 + rr12[j]*(TJP[j]-TJP[j-1])*wt[j];
				for(int kg=0; kg<NOGP; kg++){
					for(int i=0;i<NSIZE[kg];i++){
						if(kg==0){
							tmp1=1.0/v1;
							tmp2=(1.0-v1)/(v1*pow((1-u1),2.0));
							avexrr11[kg][i]=avexrr11[kg][i] + tmp1*xf1[kg][i][j]*(TJP[j+1]-TJP[j])*wt[j];
							avexrr12[kg][i]=avexrr12[kg][i] + tmp1*xf1[kg][i][j]*(TJP[j+1]-TJP[j])*wt[j];
						}
						if(kg==1){
							tmp1=0.0-u1/pow(v1,2);
							tmp2=0.0-u1/(pow(v1,2)*(1-u1));
							avexrr11[kg][i]=avexrr11[kg][i] - tmp1*xf1[kg][i][j]*(TJP[j+1]-TJP[j])*wt[j];
							avexrr12[kg][i]=avexrr12[kg][i] - tmp2*xf1[kg][i][j]*(TJP[j+1]-TJP[j])*wt[j];
						}
					}
				}
			}
		}
	}

	for(int kg=0; kg<NOGP; kg++){
		for(int i=0;i<NSIZE[kg];i++){
			sgavef1[kg] = sgavef1[kg]+pow(avexf1[kg][i],2.0);
			sgaverr11=sgaverr11 + pow(avexrr11[kg][i],2.0);
			sgaverr12=sgaverr12 + pow(avexrr12[kg][i],2.0);
		}
	}
	
	avef1d=avef1[0] - avef1[1];
    sgavef1d=sgavef1[0] + sgavef1[1];
    ave1dcil=avef1d - Z95*sqrt(sgavef1d);
    ave1dciu=avef1d + Z95*sqrt(sgavef1d);
    stat=pow(avef1d,2.0)/sgavef1d;
    pv1d=pchisq(stat,1.,0,0);
	ave[0] = avef1d; avese[0]=sqrt(sgavef1d); aveu95[0]=ave1dciu; avel95[0]=ave1dcil; avepval[0]=pv1d; 

	double aver11cil=0., aver11ciu=0.;
	double pvaverr11=0.;
	sglog11=pow((1.0/averr11),2.0)*sgaverr11;
    tmp1=log(averr11) - Z95*sqrt(sglog11);
    tmp2=log(averr11) + Z95*sqrt(sglog11);
    aver11cil=exp(tmp1);
    aver11ciu=exp(tmp2);
    stat=pow(log(averr11),2.0)/sglog11;
    pvaverr11=pchisq(stat,1.,0,0);
	ave[1] = averr11; avese[1]=sqrt(sgaverr11); aveu95[1]=aver11ciu; avel95[1]=aver11cil; avepval[1]=pvaverr11; 

	double aver12cil=0., aver12ciu=0.;
	double pvaverr12=0.;
	sglog12=pow((1.0/averr12),2.0)*sgaverr12;
    tmp1=log(averr12) - Z95*sqrt(sglog12);
    tmp2=log(averr12) + Z95*sqrt(sglog12);
    aver12cil=exp(tmp1);
    aver12ciu=exp(tmp2);
    stat=pow(log(averr12),2.0)/sglog12;
    pvaverr12=pchisq(stat,1.,0,0);
	ave[2] = averr12; avese[2]=sqrt(sgaverr12); aveu95[2]=aver12ciu; avel95[2]=aver12cil; avepval[2]=pvaverr12; 
}


void PLT::CB(int ntau1, int ntau2, double **f1, double **sgf1, double ***xf1, double sgrr11[], double sgrr12[], double *cbcut, int *n_sim){
	int NRUN = *n_sim;
	Array2D<double> gn(NOGP,N,0.);
	double wb1diff[NRUN];	
	double wbrr11[NRUN];
	double wbrr12[NRUN];
	double wmax1d, wmaxrr11, wmaxrr12;
	double xd12, sdd12;
	double tmp1, tmp2, wtmp1d, u1, v1; 
	double xrr11, xrr12, wtmprr11, wtmprr12;
	
	for(int i=0; i<NRUN; i++){
		wb1diff[i]=0;
		wbrr11[i]=0;
		wbrr12[i]=0;
	}

    GetRNGstate();
	for(int nr=0; nr<NRUN; nr++){
		for(int kg=0; kg<NOGP; kg++){
			for(int i=0; i<NSIZE[kg];i++){
				gn[kg][i] = rnorm(0.0,1.0);
			}
		}
		
		wmax1d=0.,  wmaxrr11=0., wmaxrr12=0.; 
		for(int j=ntau1; j<=ntau2; j++){
			sdd12=sqrt(sgf1[0][j]+sgf1[1][j]);
			xd12=0.;
			for(int i=0;i<NSIZE[0]; i++){
				xd12  = xd12 + xf1[0][i][j]*gn[0][i];
			}
			for(int i=0;i<NSIZE[1];i++){
				xd12  = xd12 - xf1[1][i][j]*gn[1][i];
			}
			wtmp1d=0.;
			if(sdd12>0.) wtmp1d=xd12/sdd12;   //standardized confidence band
			if(abs(wtmp1d)>wmax1d) wmax1d=abs(wtmp1d);

			u1=f1[0][j];
			v1=f1[1][j];
			xrr11=0.;
			xrr12=0.; 
			if((sgrr11[j]>0.)&&(sgrr12[j]>0.)){
				for(int i=0; i<NSIZE[0]; i++){
					tmp1=1.0/v1;
					tmp2=(1.0-v1)/(v1*pow((1.0-u1),2.0));
					xrr11=xrr11+tmp1*xf1[0][i][j]*gn[0][i];
					xrr12=xrr12+tmp2*xf1[0][i][j]*gn[0][i];
				}
				for(int i=0; i<NSIZE[1]; i++){
					tmp1=0.0-(u1/pow(v1,2));
					tmp2=0.0-u1/(pow(v1,2)*(1.0-u1));
					xrr11=xrr11+tmp1*xf1[1][i][j]*gn[1][i];
					xrr12=xrr12+tmp2*xf1[1][i][j]*gn[1][i];
				}
				wtmprr11=xrr11/sqrt(sgrr11[j]);
				wtmprr12=xrr12/sqrt(sgrr12[j]);
			}
			else{
				wtmprr11=0.;
				wtmprr12=0.;
			}
			if(abs(wtmprr11)>wmaxrr11) wmaxrr11=abs(wtmprr11);
			if(abs(wtmprr12)>wmaxrr12) wmaxrr12=abs(wtmprr12);		
		}

		wb1diff[nr]=wmax1d;
		wbrr11[nr]=wmaxrr11;
		wbrr12[nr]=wmaxrr12;
	}
	PutRNGstate();

	bubble_sort(wb1diff,NRUN);
	bubble_sort(wbrr11,NRUN);
	bubble_sort(wbrr12,NRUN);

	int icut=0.95*NRUN;	
	cbcut[0] = wb1diff[icut]; 
	cbcut[1] = wbrr11[icut];
	cbcut[2] = wbrr12[icut];
}



