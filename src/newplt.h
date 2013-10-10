//newplt.h
#ifndef NEWPLT_H
#define NEWPLT_H
#include "jamatnt/tnt.h"
#include "jamatnt/tnt_array1d.h"
#include "jamatnt/tnt_array2d.h"
#include "jamatnt/tnt_array3d.h"
#include "jamatnt/tnt_array1d_utils.h"
#include "jamatnt/tnt_array2d_utils.h"
#include "jamatnt/tnt_array3d_utils.h"

using namespace TNT;

class PLT{
	int N;					//Total number of obs
	int NOGP;
	int NJP;
	int TNJP;
	//int NRUN;				//# of simulated normal
	double P, Q;            //weight parameter
	double Z95;				//95% quantile normal
	Array1D<double> TT;      //All time
	Array1D<double> Cause;
	Array1D<double> IGP;	     //Group indictor
	Array1D<int> NSIZE;   //size for each group
	Array1D<double> TJP;  
	Array1D<double> TTJP;
	Array2D<double> T;      //time for each group
	Array3D<double> IC;     //indictor for group i, individual j, cause k
	Array1D<int> IDX;
	
public:
	PLT(double *p, double *q, int *n, double *tt, double *group, double *cause, int *size, int *njp, double *tjp);
	void CIF_est(double **ny, double **f1, double **sgf1, double ***xf1);
	void G_est(int ntau1, int ntau2, double **f1, double **sgf1, double f1diff[], double sgf1diff[], double pvf1diff[], 
			double rr11[], double sgrr11[], double pvrr11[], double rr12[], double sgrr12[], double pvrr12[]);
	void WT(double *wt, int ntau1, int ntau2);		
	void Ave_est(int ntau1, int ntau2, double **f1, double ***xf1, double rr11[], double rr12[], 
				double *wt, double *ave, double *avese, double *aveu95, double *avel95, double *avepval);
	void CB(int ntau1, int ntau2, double **f1, double **sgf1, double ***xf1, double sgrr11[], double sgrr12[], double *cbcut, int *n_sim);
	
};

#endif
