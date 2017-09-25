/* 'Exponential'regression. Göran Broström (1982-2003) */

#include <R.h>
/* #include <R_ext/Applic.h> */
#include "expreg.h"

static double e_fun(int n, double *beta, void *vex){

    int ord, ipfixed;
    Exts *ex;
    double f;
    double *dummy = NULL;


    ex = vex;

    ord = 0;
    ipfixed = 1;

    F77_CALL(wfunc)(&ord, &ipfixed, ex->pfix, &n, ex->mb, beta,
		    ex->nn, ex->z, ex->time0, ex->time, ex->ind, ex->offset,
		    &f, dummy, dummy, ex->iok);

    return(f);
    
}
    
	
static void ge_fun(int n, double *beta, double *dloglik, void *vex){

    int ord, ipfixed, i;
    Exts *ex;
    double f;
    double *dummy = NULL;
    double *fp;
/*    Rprintf("Into [ge_fun]\n"); */

    fp = Calloc(n, double);

    ex = vex;
    ord = 1;
    ipfixed = 1;

    F77_CALL(wfunc)(&ord, &ipfixed, ex->pfix, &n, ex->mb, beta,
		    ex->nn, ex->z, ex->time0, ex->time, ex->ind, ex->offset,
		    &f, fp, dummy, ex->iok);

    for (i = 0; i < n; i++) dloglik[i] = fp[i];
    Free(fp);

    return;
}
    
	
void expsup(int *iter, double *eps, int *printlevel,
	    int *nn, int *ncov, int *bdim,
	    double *time0, double *time, int * ind,
	    double *covar, double *offset, double *shape,
	    double *init, double *beta, double *lambda, double *lambda_sd,
	    double *loglik, double *dloglik, double *variance, double *sctest,
	    int *conver, int *fail){

    Exts *ex;
    int ord, i, j;
    int iok;
    int maxiter;
    int trace;
    int *mask;
    int events;
    int nREPORT = 1;
    int fncount, grcount;
    int ipfixed = 1;
    double Fmin;
    double zb, s, d, ap, bdz;

    ex = (Exts *)R_alloc(1, sizeof(Exts));
    mask = (int *)R_alloc(*bdim, sizeof(int));

    for (i = 0; i < *bdim; i++){
	mask[i] = 1;
    }

    maxiter = 1000;
    trace = *printlevel;

    iok = 0;

/* Fill in 'ex': */
    ex->pfix = shape;
    ex->mb = ncov;
    ex->nn = nn;
    ex->z = covar;
    ex->time0 = time0;
    ex->time = time;
    ex->ind = ind;
    ex->offset = offset;
    ex->iok = &iok;

    for (i = 0; i < *ncov; i++) beta[i] = init[i];
    *lambda = 0.0;
    events = 0;
    for (i = 0; i < *nn; i++){
	zb = offset[i];
	for (j = 0; j < *ncov; j++){
	    zb += beta[j] * covar[j + i * (*ncov)];
	}
	*lambda += exp(zb) * (time[i] - time0[i]);
	events += ind[i];
    }
    if (events <= 0) error("No events\n");
    if (*lambda <= 0.0) error("No (or negative) exposure time!\n");
    *lambda = (double)events / *lambda;
    beta[*ncov] = log(*lambda);

    Fmin = 0.0;
    s = 0.0;
    d = 0.0;
    bdz = 0.0;
    ap = log(*lambda);
    for (i = 0; i < *nn; i++){
	zb = offset[i];
	for (j = 0; j < *ncov; j++){
	    if (ind[i])
		bdz += beta[j] * covar[j + i * (*ncov)];
	    zb += beta[j] * covar[j + i * (*ncov)];
	}
	s += *lambda * exp(zb) * (time[i] - time0[i]);
	d += ind[i];
	Fmin += ind[i] * (ap + zb);
	Fmin -= *lambda * exp(zb) * (time[i] - time0[i]);
    }

    ord = 0; /* get initial loglik value in 'loglik[0]' */

    F77_CALL(wfunc)(&ord, &ipfixed, ex->pfix, bdim, ex->mb, beta,
		    ex->nn, ex->z, ex->time0, ex->time, ex->ind, ex->offset,
		    &Fmin, dloglik, variance, ex->iok);
    loglik[0] = -Fmin; /* NOTE! */

    vmmin(*bdim, beta, &Fmin,   
	  e_fun, ge_fun, maxiter, trace,  
	  mask, *eps, *eps, nREPORT,
	  ex, &fncount, &grcount, fail);

    loglik[1] = -Fmin;
    ord = 2;
    F77_CALL(wfunc)(&ord, &ipfixed, ex->pfix, bdim, ex->mb, beta,
		    ex->nn, ex->z, ex->time0, ex->time, ex->ind, ex->offset,
		    &Fmin, dloglik, variance, ex->iok);
 
    F77_CALL(expnr)(iter, eps, printlevel, nn, ncov, bdim,
		    time0, time, ind, covar, offset, shape,
		    beta, lambda, lambda_sd, &Fmin, dloglik,
		    variance,
		    conver, fail);
    loglik[1] = Fmin;
}
