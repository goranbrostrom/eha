#ifndef WEIBREG_H
#define WEIBREG_H

void F77_NAME(wfunc)(int *, int *, double *, int *, int *, double *,
		     int *, double *, double *, double *, int *, double *,
		     double *, double *, double *, int *);

void F77_NAME(weibnr)(int *, double *, int *, int *, int *, int *,
		      double *, double *, int *, double *, double *,
		      double *, double *, double *,
		      double *, int *, int *,
		      int *, int *);

/*
static double we_fun(int n, double *beta, void *vex);
static void gwe_fun(int n, double *beta, double *dloglik, void *vex);
*/
void sw_fun(int *order,
	    int *bdim, int *mb, double *b, 
	    int *nn, double *z, double *time0, double *time, int *ind, 
	    double *offset, int *ns, int *nstra,
      	    double *f, double *fp, double *fpp, int *iok);

void weibsup(int *iter, double *eps, int *printlevel,
	     int *ns, int *nstra, int *nn, int *ncov, int *bdim,
	     double *time0, double *time, int * ind,
	     double *covar, double *offset,
	     double *init, double *beta, double *lambda, double *lambda_sd,
	     double *shape, double *shape_sd,
	     double *loglik, double *dloglik, double *variance, double *sctest,
	     int *conver, int *fail);

typedef struct{
    int *ns;
    int *nstra;
    double *pfix;
    int *mb;
    int *nn;
    double *z;
    double *time0;
    double *time;
    int *ind;
    double *offset;
    double *f;
    double *fp;
    double *fpp;
    int *iok;
} Exts;

#endif
