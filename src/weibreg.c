#include <R.h>
#include <R_ext/Applic.h>
#include "weibreg.h"


static double we_fun(int n, double *beta, void *vex){

    int ord, ipfixed, i;
    Exts *ex;
    double f, ftot, *w_beta;
    int nn, start, mb, w_n;
    double *dummy = NULL;
    
    ex = vex;
    mb = *(ex->mb);
    w_n = mb + 2; 
    w_beta = Calloc(mb + 2, double); 

    ord = 0;
    ipfixed = 0;
    ftot = 0.0;
    
    for (i = 0; i < mb; i++) w_beta[i] = beta[i];
    for (i = 0; i < *(ex->ns); i++){
	start = ex->nstra[i];
	nn = ex->nstra[i+1] - start;
	w_beta[mb] = beta[mb + 2 * i];
	w_beta[mb + 1] = beta[mb + 2 * i + 1];
    
	F77_CALL(wfunc)(&ord, &ipfixed, ex->pfix, &w_n, ex->mb, w_beta,
			&nn, ex->z + mb * start, 
			ex->time0 + start, ex->time + start, 
			ex->ind + start, 
			ex->offset + start,
			&f, dummy, dummy, ex->iok);
	ftot += f;
    }
    Free(w_beta);
    return(ftot);
}
	
static void gwe_fun(int n, double *beta, double *dloglik, void *vex){

    int ord, ipfixed, i, j;
    double *fp, *w_beta;
    Exts *ex;
    int start, nn, mb, w_n;
    double f;
    double *dummy = NULL;

    ex = vex;

    mb = *(ex->mb);
    w_n = mb + 2;

    fp = Calloc(mb + 2, double);
    w_beta = Calloc(mb + 2, double);

    ord = 1;
    ipfixed = 0;


    for (j = 0; j < n; j++) dloglik[j] = 0.0;
    for (i = 0; i < mb; i++) w_beta[i] = beta[i];
    for (i = 0; i < *(ex->ns); i++){
	start = ex->nstra[i];
	nn = ex->nstra[i+1] - ex->nstra[i];
	w_beta[mb] = beta[mb + 2 * i];
	w_beta[mb + 1] = beta[mb + 2 * i + 1];
	F77_CALL(wfunc)(&ord, &ipfixed, ex->pfix, &w_n, ex->mb, w_beta,
			&nn, ex->z + *(ex->mb) * start,
			ex->time0 + start, ex->time + start,
			ex->ind + start, ex->offset + start,
			&f, fp, dummy, ex->iok);
	for (j = 0; j < mb; j++) dloglik[j] += fp[j];
	dloglik[mb + 2*i] += fp[mb];
	dloglik[mb + 2*i + 1] += fp[mb + 1];
    }
    Free(fp);
    Free(w_beta);
}
    
void sw_fun(int *order,
	    int *bdim, int *mb, double *beta, 
	    int *nn, double *z, double *time0, double *time, int *ind, 
	    double *offset, int *ns, int *nstra,
      	    double *f, double *fp, double *fpp, int *iok){

/* This function deals _only_ with true Weibull,   */
/* i.e., 'shape' is estimated.                     */
/* Data _must_ be sorted according to stratum.     */
/* (if ns == 1, no sorting is of course necessary) */
/* 'nstra' contains the 'cumsum' strata sizes,     */ 
/* starting with 0, length (ns + 1).               */
 
    int  i, j, m;
    double w_f, *w_fp, *w_fpp, *w_beta;
    int start, w_nn, w_bdim;
    int sc_row, sh_row;

    int ipfixed = 0;   /* 'shape', AKA 'p' is NOT fixed */ 
    double pfix = 0.0;
 
/* If no strata, do it quickly and return: */
    if (*ns == 1){ /* No strata */
	F77_CALL(wfunc)(order, &ipfixed, &pfix, bdim, mb, beta,
			nn, z,
			time0, time,
			ind, offset,
			f, fp, fpp, iok);
	return;
    }

/* If strata, a little elaboration is necessary.       */
/* We must call 'wfunc' for each stratum and add up.   */
/* Remember, each stratum has its own scale and shape. */

    w_bdim = *mb + 2;

    w_fp = Calloc(w_bdim, double);
    w_fpp = Calloc(w_bdim * w_bdim, double);
    w_beta = Calloc(w_bdim, double);

    *f = 0.0;
    for (j = 0; j < *bdim; j++) fp[j] = 0.0;
    for (j = 0; j < (*bdim) * (*bdim); j++) fpp[j] = 0.0;

    for (j = 0; j < *mb; j++) w_beta[j] = beta[j];
    for (i = 0; i < *ns; i++){
	start = nstra[i];
	w_nn = nstra[i+1] - nstra[i];
	sc_row = *mb + 2 * i;
	sh_row = sc_row + 1;
	w_beta[*mb] = beta[sc_row];
	w_beta[*mb + 1] = beta[sh_row];
	F77_CALL(wfunc)(order, &ipfixed, &pfix, &w_bdim, mb, w_beta,
			&w_nn, z + (*mb) * start,
			time0 + start, time + start,
			ind + start, offset + start,
			&w_f, w_fp, w_fpp, iok);
	*f += w_f;
	for (j = 0; j < *mb; j++) fp[j] += w_fp[j];
	fp[*mb + 2*i] += w_fp[*mb];
	fp[*mb + 2*i + 1] += w_fp[*mb + 1];
/* And then the hessian, a little bit complicated: */
	for (j = 0; j < *mb; j++){
	    for (m = 0; m <= j; m++){
		fpp[m + j * (*bdim)] += w_fpp[m + j * w_bdim];
	    }
	}

	for (m = 0; m < *mb; m++){
	    fpp[m + sc_row * (*bdim)] += w_fpp[m + (*mb) * w_bdim];
	    fpp[m + sh_row * (*bdim)] += w_fpp[m + ((*mb) + 1) * w_bdim];
	}

	fpp[sc_row + sc_row * (*bdim)] += w_fpp[*mb + (*mb) * w_bdim];
	fpp[sc_row + sh_row * (*bdim)] += w_fpp[*mb + ((*mb) + 1) * w_bdim];
	fpp[sh_row + sh_row * (*bdim)] += 
	    w_fpp[*mb + 1 + ((*mb) + 1) * w_bdim];
	
    }

    /* Fill in "the other" half: */
    for (j = 0; j < *bdim; j++){
	for (m = (j + 1); m < *bdim; m++){
	    fpp[m + j * (*bdim)] = fpp[j + m * (*bdim)];
	}
    }

    Free(w_fp);
    Free(w_fpp);
    Free(w_beta);
}


	
void weibsup(int *iter, double *eps, int *printlevel,
	     int *ns, int *nstra, int *nn, int *ncov, int *bdim,
	     double *time0, double *time, int * ind,
	     double *covar, double *offset,
	     double *init, double *beta, double *lambda, double *lambda_sd,
	     double *shape, double *shape_sd,
	     double *loglik, double *dloglik, double *variance, double *sctest,
	     int *conver, int *fail){

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
    /*    int zero = 0;
    int two = 2;
    double pfix = 1.0;
    */
    ex = (Exts *)R_alloc(1, sizeof(Exts));
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
    *lambda = (double)events / *lambda;
    for (i = 0; i < *ns; i++){
	j = *ncov + 2 * i; 
	xin[j] = log(*lambda);
	beta[j] = log(*lambda);
    }

/* Done with initial values; the same p and lambda in all strata. */

    /* Temporary hack for check: */
/*
    for (i = 0; i < *bdim; i++) xin[i] = 0.1;
    xin[*bdim - 2] = 0.0;
    Fmin = we_fun(*bdim, xin, ex);
    gwe_fun(*bdim, xin, dloglik, ex);
    F77_CALL(wfunc)(&two, &zero, &pfix, bdim, ncov, xin,
		    nn, covar, time0, time, ind, offset,
		    &Fmin, dloglik, variance, &iok);
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
    nmmin(*bdim, xin, beta, &Fmin, we_fun, 
	  fail, *eps, *eps, ex, 
	  1.0, 0.5, 2.0, trace,
	  &fncount, maxiter);
    trace = *printlevel;

    loglik[1] = -Fmin;
    gwe_fun(*bdim, beta, dloglik, ex);
    Rprintf("\nEfter 'nmmin': loglik = %f\n", -Fmin);
    Rprintf(" beta och dloglik:\n");
    for (i = 0; i < *bdim; i++){
	Rprintf("%f, %f\n", beta[i], dloglik[i]);
    }
*/

/* 'Mask out' the regression coefficients: */
    for (i = 0; i < *ncov; i++){
	mask[i] = 0;
    }

    for (i = *ncov; i < *bdim; i++){
	mask[i] = 1;
    }

/* Estimate the 'null' model (only scale and shape): */
    vmmin(*bdim, beta, &Fmin,
	  we_fun, gwe_fun, maxiter, trace,  
	  mask, *eps, *eps, nREPORT,
	  ex, &fncount, &grcount, fail);

    if (trace)
	Rprintf("\nOnly scale and shape: loglik = %f\n", -Fmin);

    loglik[0] = -Fmin;

/* Now, 'mask in' the regression coefficients: */
    for (i = 0; i < *bdim; i++){
	mask[i] = 1;
    }
/* Estimate the full model: */
    vmmin(*bdim, beta, &Fmin,
	  we_fun, gwe_fun, maxiter, trace,  
	  mask, *eps, *eps, nREPORT,
	  ex, &fncount, &grcount, fail);

    if (trace)
      Rprintf("\nAfter 'vmmin': loglik = %f\n", -Fmin);
    loglik[1] = -Fmin;

    gwe_fun(*bdim, beta, dloglik, ex);
    
    if (trace){
      Rprintf("\nEfter 'vmmin': loglik = %f\n", -Fmin);
      Rprintf(" beta och dloglik:\n");
      for (i = 0; i < *bdim; i++){
	Rprintf("%f, %f\n", beta[i], dloglik[i]);
      }
    }

/*
    ord = 1;

    F77_CALL(wfunc)(&ord, &ipfixed, ex->pfix, bdim, ex->mb, beta,
		    ex->nn, ex->z, ex->time0, ex->time, ex->ind, ex->offset,
		    ex->f, ex->fp, ex->fpp, ex->iok);
*/


/*
    ord = 2;

    sw_fun(&ord, 
	   bdim, ncov, beta, 
	   nn, covar, time0, time, ind, 
	   offset, ns, nstra,
	   &Fmin, dloglik, variance, fail);

    loglik[1] = -Fmin;
*/

/* Take some Newton-Raphson steps: */
    F77_CALL(weibnr)(iter, eps, printlevel, nn, ncov, bdim,
		     time0, time, ind, covar, offset,
		     beta, loglik + 1, dloglik,
		     variance, ns, nstra,
		     conver, fail);
    if (trace){
      Rprintf("Variance (in [weibreg]) after N-R:\n");
      for (i = 0; i < *bdim; i++){
	for (j = 0; j < *bdim; j++){
	  Rprintf("%f ", variance[i + j * (*bdim)]);
	}
	Rprintf("\n");
      }

    Rprintf("Score: ");
    for (i = 0; i < *bdim; i++) Rprintf("%f ", dloglik[i]);
    Rprintf("\n");
    }
    if (trace){
	Rprintf("\nAfter Newton-Raphson: loglik = %f\n", loglik[1]);
	Rprintf("fail = %d\n", *fail);
    }
 }
