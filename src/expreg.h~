#ifndef EXPREG_H
#define EXPREG_H

void F77_NAME(wfunc)(int *, int *, double *, int *, int *, double *,
		     int *, double *, double *, double *, int *, double *,
		     double *, double *, double *, int *);

void F77_NAME(expnr)(int *, double *, int *, int *, int *, int *,
		    double *, double *, int *, double *, double *, double *,
		    double *, double *, double *, double *, double *,
		    double *,
		    int *, int *);

static double e_fun(int n, double *beta, void *vex);
static void ge_fun(int n, double *beta, double *dloglik, void *vex);

void expsup(int *iter, double *eps, int *printlevel,
	    int *nn, int *ncov, int *bdim,
	    double *time0, double *time, int * ind,
	    double *covar, double *offset, double *shape,
	    double *init, double *beta, double *lambda, double *lambda_sd,
	    double *loglik, double *dloglik, double *variance, double *sctest,
	    int *conver, int *fail);

typedef double optimfn(int n, double *par, void *ex);
typedef void optimgr(int n, double *par, double *gr, void *ex);

void vmmin(int n, double *x, double *Fmin,
	   optimfn fn, optimgr gr, int maxit, int trace,
	   int *mask, double abstol, double reltol, int nREPORT,
	   void *ex, int *fncount, int *grcount, int *fail);

typedef struct{
    double *pfix;
    int *mb;
    int *nn;
    double *z;
    double *time0;
    double *time;
    int *ind;
    double *offset;
    int *iok;
} Exts;

#endif
