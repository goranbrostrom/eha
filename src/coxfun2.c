#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif
#include <R.h>
#include <Rmath.h>
/* #include "sup.h" */
#include "eha_zeroin.h"
#include "coxfun2.h"

/* extern RS_fun *eha_rs; */

static void cox_obs_rs2(int what, int antevents, double *weights, double *lin,
                        double *x, int p, double *b,
                        /* Return: */
                        double *loglik, double *dloglik){
    
    int i;
    /* double one = 1.0; */
    double wght;
    int ione = 1;
    for (i = 0; i < antevents; i++){
	wght = weights[i];
	*loglik += wght * lin[i];
	if (what >= 1) F77_CALL(daxpy)(&p, &wght, (x + p * i), &ione,
				       dloglik, &ione);
    }
}
    
void breslow_rs2(int *what, /* RiskSet *risk, */
		 int *antevents,
		 int *size,
		 double *weights,
		 double *x,
		 double *lin,
		 int *p,
		 double *b, double *e_frac,
		 /* Return: */
		 double *loglik, double *dloglik, 
		 double *d2loglik){

/* Calculates the "breslow" variation of the partial likelihood, */
/* and eventually its first and second partial derivatives.      */

    int i;
    char up = 'U';

    double zero = 0.0;
    int izero = 0;
    int ione = 1;
    double sumscore;
    double *sumdscore;
    double *sumd2score;
    
    int p2 = (*p) * (*p);
    double alpha;
    double *wsc;
    
    /* if (risk->out) return; */

    /* First the "observed" part (common to the 'efron' method): */
    
    cox_obs_rs2(*what, /* risk, */
		*antevents,
		/* int *eventset, */
		weights,
		lin,
		x,
		*p,
		b,
		loglik, dloglik);
    
    /* Then the "expected" part: */

    /* Initialize: */
    sumdscore = Calloc(*p, double);
    sumd2score = Calloc(p2, double);
    wsc = Calloc(*size, double);
    sumscore = 0.0;
    if (*what >= 1){
	F77_CALL(dcopy)(p, &zero, &izero, sumdscore, &ione);
	if (*what >= 2){
	    F77_CALL(dcopy)(&p2, &zero, &izero, sumd2score, &ione);
	}
    }

    /* Go thru riskset: */

    for (i = 0; i < *size; i++){
	wsc[i] = weights[i] * exp(lin[i]);
	/* sumscore += score[who]; */
	sumscore += wsc[i];
	if (*what >= 1){ /* First derivatives: */
	    F77_CALL(daxpy)(p, (wsc + i), (x + (*p) * i), &ione, 
				sumdscore, &ione);

	    if (*what >= 2){ /* Second derivatives: */
		F77_CALL(dsyr)(&up, p, (wsc + i), 
			       (x + (*p) * i), &ione,
			       sumd2score, p FCONE);
	    }
	}
    }

    /* Add in: */

    *loglik -= *antevents * log(sumscore);/*KOLLA!*/
    if (*what >= 1){
	alpha = -(double)(*antevents) / sumscore;
	F77_CALL(daxpy)(p, &alpha, sumdscore, &ione, 
			dloglik, &ione);
	if (*what >= 2){
	    alpha = -alpha;
	    F77_CALL(daxpy)(&p2, &alpha, sumd2score, &ione,
			    d2loglik, &ione);
	    alpha = -alpha / sumscore;
	    F77_CALL(dsyr)(&up, p, &alpha, sumdscore, &ione,
			   d2loglik, p FCONE);
	}
    }
    Free(wsc);
    Free(sumd2score);
    Free(sumdscore);
}

void efron_rs2(int *what, /* RiskSet *risk, */
	       int *antevents,
	       int *size,
	       double *weights,
	       double *x, /* size * p */
	       double *lin,
	       int *p,
	       double *b, double *e_frac,
	       /* Return: */
	       double *loglik, double *dloglik, 
	       double *d2loglik){
  
    /*************************************************************
C     This subroutine calculates the 'efron' variation of the
C     log likelihood function, and its first and second order
C     partial derivatives. In ONE risk set!
    *************************************************************/

    int i, r, who;
    double sumscore = 0.0;

/*************************************************************
C     +++
C     Local (note the deviation from strict standard here!):
**************************************************************/

    double escore;
    double *edscore, *ed2score;
    double w, ws;

    char up = 'U';

    double *temp; /*(antcov)*/
    int p2 = (*p) * (*p);

    int izero = 0;
    int ione = 1;
    
    double *sumdscore;
    double *sumd2score;
    
    double zero = 0.0;
    double one = 1.0;

    double alpha;

    double *wsc;

    /* if (risk->out) return; */

    /* First the "observed" part (common to the 'breslow' method): */

    cox_obs_rs2(*what, /* risk, */
		*antevents,
		/* eventset, */
		weights,
		lin,
		x,
		*p,
		b,
		/* Return: */
		loglik, dloglik);

    /* Then the "expected" part: */

    sumdscore = Calloc(*p, double);
    sumd2score = Calloc(p2, double);
    wsc = Calloc(*size, double);
    edscore = Calloc(*p, double);
    ed2score = Calloc(p2, double);

    temp = Calloc(*p, double);
		     

/*     Reset to zero: */
    sumscore = zero;
    escore = zero;
    if (*what >= 1){
	F77_CALL(dcopy)(p, &zero, &izero, sumdscore, &ione);
        F77_CALL(dcopy)(p, &zero, &izero, edscore, &ione);
	if (*what >= 2){
	    F77_CALL(dcopy)(&p2, &zero, &izero, sumd2score, &ione);
	    F77_CALL(dcopy)(&p2, &zero, &izero, ed2score, &ione);
	}
    }

/*     Go thru riskset(rs, j): */
    for (i = 0; i < *size; i++){
	who = i;
	wsc[i] = weights[i] * exp(lin[i]);
	/* sumscore += score[who]; */
	sumscore += wsc[i];
	if (*what >= 1){
	    F77_CALL(daxpy)(p, (wsc + i), (x + (*p) * who),
			    &ione, sumdscore, &ione);
	    if (*what >= 2){
		F77_CALL(dsyr)(&up, p, (wsc + i),
                               (x + (*p) * who), &ione, 
                               sumd2score, p FCONE);
/*		F77_CALL(dger)(&p, &p, (score + who), (x + p * who), &ione,
		(x + p * who), &ione, sumd2score, &p); */
		
	    }
	}
    }
               
    if (*antevents == 1){ /* No ties */
/*     Add into loglik: */
	*loglik -= log(sumscore);
	if (*what >= 1){
/*     Add into dloglik: */
	    /*  alpha = -one / sumscore; */
	    alpha = - 1.0 / sumscore; /* 090405 */
	    F77_CALL(daxpy)(p, &alpha, 
			    sumdscore, &ione, dloglik, &ione);
	    if (*what >= 2){
/*     Add into d2loglik: */
		alpha = 1.0 / sumscore; 
		F77_CALL(daxpy)(&p2, &alpha, sumd2score, &ione, 
				d2loglik, &ione);
		alpha = -1.0 /  (sumscore * sumscore);
		F77_CALL(dsyr)(&up, p, &alpha, sumdscore, &ione, 
		d2loglik, p FCONE); 
/*		F77_CALL(dger)(&p, &p, &alpha, sumdscore, &ione, 
			       sumdscore, &ione, d2loglik, &p); */
	    }
	}
    }else{
/*     +++ IF TIES: */
	
/* Go thru events and create escore, edscore and ed2score: */
	for (i = 0; i < *antevents; i++){
	    who = i;
	    escore = escore + wsc[i];
	    if (*what >= 1){ /* first derivatives */
		F77_CALL(daxpy)(p, (wsc + i), (x + (*p) * who), &ione,
				edscore, &ione);
		if (*what >= 2){ /* second derivatives */
		    F77_CALL(dsyr)(&up, p, (wsc + i), 
				   (x + (*p) * who), &ione,
				   ed2score, p FCONE); 
/*	       F77_CALL(dger)(&p, &p, (score + who), (x + p * who), &ione,
	       (x + p * who), &ione, ed2score, &p); */
		}
	    }
	}


	for (r = 0; r < *antevents; r++){
/* WRONG: (Fortran)   w = (double)(r - ione) / (double)(risk->antevents);*/
	    w = (double)r / (double)(*antevents);
	    ws = w * escore;
/*     Add into loglik: */
	    *loglik -= log(sumscore - ws);
	    if (*what >= 1){
/*     Add into dloglik: */
		F77_CALL(dcopy)(p, sumdscore, &ione, temp, &ione);
		alpha = -w;
		F77_CALL(daxpy)(p, &alpha, edscore, &ione, temp, &ione);
		alpha = one / (sumscore - ws);
		F77_CALL(dscal)(p, &alpha, temp, &ione);
		alpha = -one;
		F77_CALL(daxpy)(p, &alpha, temp, &ione, dloglik, &ione);
		if (*what >= 2){
/*     Add into d2loglik: */
		    alpha = one / (sumscore - ws);
		    F77_CALL(daxpy)(&p2, &alpha, sumd2score, &ione,
				    d2loglik, &ione);
		    alpha = -w / (sumscore - ws);
		    F77_CALL(daxpy)(&p2, &alpha, ed2score, &ione,
				    d2loglik, &ione);
		    alpha = -one;
		    F77_CALL(dsyr)(&up, p, &alpha, temp, &ione,
				   d2loglik, p FCONE);
		}
	    }
	}
    }
    Free(temp);
    Free(ed2score);
    Free(ed2score);
    Free(wsc);
}

/* Maybe later ...

void mppl_rs2(int what, RiskSet *risk,
	     double *b, double e_frac,
	     double *loglik, double *dloglik, double *d2loglik){

    if (risk->antevents == risk->size) return;
    if (risk->out) return;
    if (risk->antevents == 1){
	breslow_rs(what, risk,
		   b, e_frac,
		       loglik, dloglik, d2loglik);
    }else if (risk->antevents <= e_frac * risk->size){
	efron_rs(what, risk, 
		 b, e_frac,
		 loglik, dloglik, d2loglik);
    }else{
	ml_rs(what, risk,
	      b, e_frac,
	      loglik, dloglik, d2loglik);
    }
}
*/
