#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>

#include "phreg.h"
#include "phfun.h"
#include "loglik_phexp.h"

/*
int dist;      
ph0_fun *S0;   
ph0_fun *f0;   
ph0_fun *h0;   
ph0_fun *f0_t; 
ph0_fun *h0_t;
ph0_fun *h0_tt;
*/

static double phexp_fun(int n, double *beta, void *vex){
/* Here shape is FIXED, i.e. NOT estimated, and NOT part of beta! */
    /* the fixed shape is found in ex->pfix. */

    int i;
    Exts *ex;
    double f, alpha, gamma, res;
    int nn, start, mb;
    
    ex = vex;
    mb = *(ex->mb);

    f = 0.0;
    
    for (i = 0; i < *(ex->ns); i++){
	start = ex->nstra[i];
	nn = ex->nstra[i+1] - start;
	alpha = beta[mb + i];     /* log(scale) */
	gamma = log(ex->pfix[i]); /* log(shape) */

	loglik_phexp(&dist,
		     &mb, beta, &alpha, &gamma,
		     &nn, ex->z + mb * start, 
		     ex->time0 + start, ex->time + start, 
		     ex->ind + start, 
		     ex->offset + start,
		     &res);
	f += res;
    }
    return(f);
}
	
static void gphexp_fun(int n, double *beta, double *dloglik, void *vex){
/* Here shape is FIXED, i.e. NOT estimated, and NOT part of beta! */
    int i, j;
    double alpha, gamma;
    double *fp;
    Exts *ex;
    int start, nn, mb;

    ex = vex;

    mb = *(ex->mb);

    fp = Calloc(mb + 1, double);

    for (j = 0; j < n; j++) dloglik[j] = 0.0;
    for (i = 0; i < *(ex->ns); i++){
	start = ex->nstra[i];
	nn = ex->nstra[i+1] - ex->nstra[i];
	alpha = beta[mb + i];
	gamma = log(ex->pfix[i]);
	d_loglik_phexp(&dist, &mb, beta, &alpha, &gamma,
		       &nn, ex->z + *(ex->mb) * start,
		       ex->time0 + start, ex->time + start,
		       ex->ind + start, ex->offset + start, fp);
	for (j = 0; j < mb; j++) dloglik[j] += fp[j];
	dloglik[mb + i] += fp[mb];
	/* dloglik[mb + 2*i + 1] += fp[mb + 1]; inte med */
    }
    Free(fp);
}
	
static void g2phexp_fun(int n, double *beta, double *d2loglik, void *vex){
/* Here shape is FIXED, i.e. NOT estimated, and NOT part of beta! */
/* n is the dimension of beta, i.e., mb + ns. */

    int i, j, m, bdim;
    double alpha, gamma;
    double *fpp;
    Exts *ex;
    int start, nn, mb;

    bdim = n; /* = mb + ns */
    ex = vex;

    mb = *(ex->mb);

    fpp = Calloc((mb + 1) * (mb + 1), double);

    for (j = 0; j < bdim * bdim; j++) d2loglik[j] = 0.0;
    if (*(ex->ns) != 1) error("Stratification not allowed here\n");
    for (i = 0; i < *(ex->ns); i++){
	start = ex->nstra[i];
	nn = ex->nstra[i+1] - ex->nstra[i];
	alpha = beta[mb + i];
	gamma = log(ex->pfix[i]);
	d2_loglik_phexp(&dist, &mb, beta, &alpha, &gamma,
			&nn, ex->z + *(ex->mb) * start,
			ex->time0 + start, ex->time + start,
			ex->ind + start, ex->offset + start, fpp);
	for (j = 0; j < mb; j++){
/* ToDo: Only correct for "No strata": */                        /* NOTE!!! */
	    d2loglik[j + mb * bdim] = fpp[j + mb * bdim];
	    d2loglik[mb + j * bdim] = fpp[mb + j * bdim];
/*
	    d2loglik[j + (mb + 1) * bdim] = fpp[j + (mb + 1) * (mb + 2)];
	    d2loglik[mb + 1 + j * bdim] = fpp[mb + 1 + j * (mb + 2)];
*/
	    for (m = 0; m < mb; m++){
		d2loglik[j + m * bdim] += fpp[j + m * (mb + 1)];
	    }
	}
	/* Måste kollas !!! */
	d2loglik[mb + 2*i + (mb + 2*i) * bdim] += fpp[mb + mb * (mb + 1)];
/*	
	d2loglik[mb + 2*i + 1 + (mb + 2 * i + 1) * bdim] += 
	    fpp[mb + 1 + (mb + 1) * (mb + 2)];
*/
	/* scale x shape: */
/*
	d2loglik[mb + 2*i + (mb + 2 * i + 1) * bdim] += 
	    fpp[mb + (mb + 1) * (mb + 2)];
	d2loglik[mb + 2*i + 1 + (mb + 2 * i) * bdim] += 
	    fpp[mb + 1 + mb * (mb + 2)];
*/
    }
    Free(fpp);
}

static void phexp_nr(int iter, double eps, int printlevel,
		  int bdim, double *beta, 
		  double *loglik, double *dloglik, double *variance,
		  int *conver, int *fail, Exts *ex){
    
    int ione = 1;
    int job = 1; /* Gives inverse from dpodi */ 
    int	itmax, info;
    double one = 1.0;
    double L2;
    double *det, *db;
    int i, j;
    void *vex;

    vex = ex;

    det = Calloc(2, double);
    db = Calloc(bdim, double);

    itmax = iter;
    iter = 0;

    *conver = 0;

    *loglik = phexp_fun(bdim, beta, vex);
    gphexp_fun(bdim, beta, dloglik, vex);
    g2phexp_fun(bdim, beta, variance, vex);

    for (i = 0; i < bdim; i++){
	dloglik[i] = -dloglik[i];
    }

    *loglik = -*loglik;
    
    while ((iter < itmax) & !(*conver)){
	F77_CALL(dcopy)(&bdim, dloglik, &ione, db, &ione); 
	F77_CALL(dpofa)(variance, &bdim, &bdim, &info);
	if (!info){
            F77_CALL(dposl)(variance, &bdim, &bdim, db);
	}else{
	    Rprintf("fail in [dpofa]; info = %d\n", info); 
	    *fail = info;
            return;
	}
	L2 = F77_CALL(dnrm2)(&bdim, db, &ione);
	if (L2 < eps) *conver = 1;
	if (printlevel){
            Rprintf("*** Iteration %d\n", iter);
            Rprintf("L2 = %f\n", L2);
            Rprintf("loglik = %f\n", *loglik);
	}
/* Update beta: */
	F77_CALL(daxpy)(&bdim, &one, db, &ione, beta, &ione);
/* Calculate f, fp, fpp at new beta: */
	*loglik = phexp_fun(bdim, beta, vex);
	gphexp_fun(bdim, beta, dloglik, vex);
	g2phexp_fun(bdim, beta, variance, vex);

	for (i = 0; i < bdim; i++) {
	    dloglik[i] = -dloglik[i];
	}
	*loglik = -*loglik;
	
/* Next iteration: */
	iter++;
    }
/* Done! The afterwork: */
    if (printlevel){
      Rprintf("Hessian (in [phnr]):\n");
      for (i = 0; i < bdim; i++){
	for (j = 0; j < bdim; j++){
	  Rprintf("%f ", variance[i + j * (bdim)]);
	}
	Rprintf("\n");
      }
    }


    F77_CALL(dpofa)(variance, &bdim, &bdim, &info);
    if (!info){
	F77_CALL(dpodi)(variance, &bdim, &bdim, det, &job);
	for (i = 1; i < bdim; i++){
	    for (j = 0; j < i; j++){
		variance[i + j * bdim] = variance[j + i * bdim];
	    }
	}
    }else{
	*fail = info;
	Rprintf("Inverse failed [dpodi 2]; info = %d\n", info);
	return;
    } 

    if (printlevel){
	Rprintf("*** Iteration %d\n", iter);
	if (*conver){
	    Rprintf("[nr_ph] Convergence!\n");
	}else{
	    Rprintf("[nr_ph] NO Convergence!\n");
	}
	Rprintf("loglik = %f\n", *loglik);
    }

    Free(db);
    Free(det);
}

void phexpsup(int *iter, double *eps, int *printlevel,
	      int *ns, int *nstra, int *nn, int *ncov, int *bdim,
	      double *time0, double *time, int * ind,
	      double *covar, double *offset, double *shape, int *dis, 
	      double *init, double *beta, double *lambda, double *lambda_sd,
	      double *loglik, double *dloglik, double *variance, double *sctest,
	      int *conver, int *fail){
/* Here shape is FIXED, i.e. NOT estimated, and NOT part of beta! */
    Exts *ex;
    int i, j;
    int iok;
    int maxiter;
    int trace;
    int *mask;
    int events;
    int nREPORT = 1;
    int fncount, grcount;
    double *xin;
    double Fmin;

    void *vex;

    int info;
    int job = 1;
    double * det;

    det = Calloc(2, double);
    dist = *dis;

    if (dist == 0){
	S0 = &S0_weibull;
	f0 = &f0_weibull;      
	h0 = &h0_weibull;      
	f0_t = &f0_t_weibull;     
	h0_t = &h0_t_weibull;     
	h0_tt = &h0_tt_weibull;    
    }else if (dist == 1){
	S0 = &S0_loglogistic;
	f0 = &f0_loglogistic;      
	h0 = &h0_loglogistic;      
	f0_t = &f0_t_loglogistic;     
	h0_t = &h0_t_loglogistic;     
	h0_tt = &h0_tt_loglogistic;
    }else if (dist == 2){    
	S0 = &S0_lognormal;
	f0 = &f0_lognormal;      
	h0 = &h0_lognormal;      
	f0_t = &f0_t_lognormal;     
	h0_t = &h0_t_lognormal;     
	h0_tt = &h0_tt_lognormal;   
    }else if ((dist == 3) || (dist == 4)){    
	S0 = &S0_ev;
	f0 = &f0_ev;      
	h0 = &h0_ev;      
	f0_t = &f0_t_ev;     
	h0_t = &h0_t_ev;     
	h0_tt = &h0_tt_ev;   
    }else{
	error("Unknown distribution");
    }

    ex = (Exts *)R_alloc(1, sizeof(Exts));
    vex = ex; /* NOTE!!! */
    mask = (int *)R_alloc(*bdim, sizeof(int));
    xin = (double *)R_alloc(*bdim, sizeof(double));

    maxiter = 1000;
    trace = *printlevel;

    iok = 0;

/* Fill in 'ex': */
    ex->ns = ns;
    ex->nstra = nstra;
    ex->pfix = shape;
    ex->mb = ncov;
    ex->nn = nn;
    ex->z = covar;
    ex->time0 = time0;
    ex->time = time;
    ex->ind = ind;
    ex->f = loglik + 1; /* Note: 'loglik' is a vector of length 2 */ 
    ex->fp = dloglik;
    ex->fpp = variance;
    ex->offset = offset;
    ex->iok = &iok;

/* Calculate simple start values for first call to 'vmmin': */
    for (i = 0; i < *ncov; i++) {
	xin[i] = init[i];
	beta[i] = init[i];
    }

    for (i = *ncov; i < *bdim; i++){
	xin[i] = 0.0;
	beta[i] = 0.0;
    }

    *lambda = 0.0;
    events = 0;
    for (i = 0; i < *nn; i++){
	*lambda += time[i] - time0[i];
	events += ind[i];
    }

    if (events <= 0) error("No events\n");
    if (*lambda <= 0.0) error("No (or negative) exposure time!\n");
    *lambda = *lambda / (double)events;
    for (i = 0; i < *ns; i++){
	j = *ncov + i; 
	xin[j] = log(*lambda);
	beta[j] = log(*lambda);
    }

/* Done with initial values; the same p and lambda in all strata. */

    /* Temporary hack for check: */
/*
    for (i = 0; i < *bdim; i++) xin[i] = 0.1;
    xin[*bdim - 2] = 0.0;
    Fmin = phexp_fun(*bdim, xin, vex);
    gphexp_fun(*bdim, xin, dloglik, vex);
    g2phexp_fun(*bdim, xin, variance, vex);
    Rprintf("\n!!b = 0.1:\n");
    Rprintf("f = %f\n", Fmin);
    Rprintf("fp: ");
    for (i = 0; i < *bdim; i++){
	Rprintf("%f ", dloglik[i]);
    }

    Rprintf("\n\n fpp = \n");
    for (i = 0; i < *bdim; i++){
	for (j = 0; j < *bdim; j++){
	    Rprintf("%f ", variance[j + *bdim * i]);
	}
	Rprintf("\n");
    }
*/
/*
    trace = 0; 
    nmmin(*bdim, xin, beta, &Fmin, phexp_fun, 
	  fail, *eps, *eps, vex, 
	  1.0, 0.5, 2.0, trace,
	  &fncount, maxiter);
    trace = *printlevel;

    Fmin = phexp_fun(*bdim, beta, vex);
    gphexp_fun(*bdim, beta, dloglik, vex);
    loglik[1] = -Fmin;
    if (trace){
	Rprintf("\nEfter 'nmmin' [phexpreg]: loglik = %f\n", -Fmin);
	Rprintf(" beta och dloglik:\n");
	for (i = 0; i < *bdim; i++){
	    Rprintf("%f, %f\n", beta[i], dloglik[i]);
	}
    }
*/

/* 'Mask out' the regression coefficients: */
    for (i = 0; i < *ncov; i++){
	mask[i] = 0;
    }
    if ((dist == 1) || (dist == 2)) mask[0] = 1; /*Intercept */ 
    for (i = *ncov; i < *bdim; i++){
	mask[i] = 1;
    }

/* Estimate the 'null' model (only scale and shape): */
    vmmin(*bdim, beta, &Fmin,
	  phexp_fun, gphexp_fun, maxiter, trace,  
	  mask, *eps, *eps, nREPORT,
	  vex, &fncount, &grcount, fail);

    if (trace)
	Rprintf("\nOnly scale and shape: loglik = %f\n", -Fmin);

    loglik[0] = -Fmin;

/* Now, 'mask in' the regression coefficients: */
    for (i = 0; i < *bdim; i++){
	mask[i] = 1;
    }
/* Estimate the full model: */
    vmmin(*bdim, beta, &Fmin,
	  phexp_fun, gphexp_fun, maxiter, trace,  
	  mask, *eps, *eps, nREPORT,
	  vex, &fncount, &grcount, fail);
    if (trace)
	Rprintf("\nAfter 'vmmin': loglik = %f\n", -Fmin);
    loglik[1] = -Fmin;
    if (trace){
	Rprintf("\n[phexpreg] After vmmin; beta is\n");
	for (i = 0; i < *bdim; i++){
	    Rprintf("%f ", beta[i]);
	}
	Rprintf("\n\n");


	Rprintf("Fmin efter vmmin: %f\n", Fmin);
    }
/*
    nmmin(*bdim, xin, beta, &Fmin, phexp_fun, 
	  fail, *eps, *eps, vex, 
	  1.0, 0.5, 2.0, trace,
	  &fncount, maxiter);
    Rprintf("Fmin efter nmmin: %f\n", -Fmin);
    Rprintf("\n[phexpreg] After Nmmin; beta is\n");
    for (i = 0; i < *bdim; i++){
	Rprintf("%f ", beta[i]);
    }
    Rprintf("\n\n");
*/	
    gphexp_fun(*bdim, beta, dloglik, vex);
    if (trace){
	Rprintf("\n[phexpreg] After vmmin; score is\n");
	for (i = 0; i < *bdim; i++){
	    Rprintf("%f ", dloglik[i]);
	}
	Rprintf("\n\n");
	
    }
    g2phexp_fun(*bdim, beta, variance, vex);
    if (trace){
	Rprintf("Hessian (in [phexpreg]):\n");
	for (i = 0; i < *bdim; i++){
	    for (j = 0; j < *bdim; j++){
		Rprintf("%f ", variance[i + j * (*bdim)]);
	    }
	    Rprintf("\n");
	}
	Rprintf("\n");
    }
/* Testing... */
    job = 1;
    F77_CALL(dpofa)(variance, bdim, bdim, &info);
    if (!info){
	F77_CALL(dpodi)(variance, bdim, bdim, det, &job);
	for (i = 1; i < *bdim; i++){
	    for (j = 0; j < i; j++){
		variance[i + j * (*bdim)] = variance[j + i * (*bdim)]; 
	    }
	}
    }else{
	*fail = info;
	Rprintf("Inverse failed [phexpreg]; info = %d\n", info);
	return;
    } 
/* End testing */

/* Take some Newton-Raphson steps: */


    phexp_nr(*iter, *eps, *printlevel, 
	  *bdim, beta, 
	  (loglik + 1), dloglik, variance, 
	  conver, fail, ex);
    if (trace){
	Rprintf("Fmin efter Newton-Raphson: %f\n", -loglik[1]);
	Rprintf("\n[phexpreg] After Newton-Raphson; beta is\n");
	for (i = 0; i < *bdim; i++){
	    Rprintf("%f ", beta[i]);
	}
	Rprintf("\n\n");
    }

    if (trace){
	Rprintf("Variance (in [phexpreg]) after N-R:\n");
	for (i = 0; i < *bdim; i++){
	    for (j = 0; j < *bdim; j++){
		Rprintf("%f ", variance[i + j * (*bdim)]);
	    }
	    Rprintf("\n");
	}
	
	Rprintf("\nScore: ");
	for (i = 0; i < *bdim; i++) Rprintf("%f ", dloglik[i]);
	Rprintf("\n");
	Rprintf("\nbeta: ");
	for (i = 0; i < *bdim; i++) Rprintf("%f ", beta[i]);
	Rprintf("\n");
	if (trace){
	    Rprintf("\nAfter Newton-Raphson: loglik = %f\n", loglik[1]);
	    Rprintf("fail = %d\n", *fail);
	}
    }

    Free(det);
}
