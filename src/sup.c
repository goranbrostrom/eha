#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "sup.h"
#include "coxfun.h"

/* Global variables */

/* Fixed: */
int p = 0;
int nn = 0;
double *x = 0; /* Covariates: (nn x p) */
/* double *offset = 0; */
/* Variable: */ 
double *lin = 0;
/* double *score = 0; */
double *sumdscore = 0;
double *sumd2score = 0;

RS_fun *eha_rs;

/* End global variables */

static void get_prob(int nn, 
                     int size, 
                     int *riskset, 
                     double *score, 
                     double *prob){

    double totsize;
    int i, who;

    totsize = 0.0;
    for (i = 0; i < size; +i++){
        who = riskset[i];
        prob[i] = score[who];
        totsize += prob[i];
    }
    for (i = 0; i < size; i++){
        prob[i] /= totsize;
    }

/* Sort prob in decreasing order, apply same permutation to riskset: */

    revsort(prob, riskset, size);
}

static void sample_events(int antevents, 
                          int size, 
                          int *riskset, 
                          int *eventset,
                          double *prob){
/***************************************************************
C     This is an adaption  of the C function 
C     ProbSampleNoReplace found in random.c 
C     in the  R  sources. This is not exactly 
C     the correct procedure! Will be fixed (some day ...)
****************************************************************/

    int k, n, i, j;
    double z, rt, mass, totalmass;
    int *perm;
    double *p;

    perm = Calloc(size, int);
    p = Calloc(size, double);
    
    if (antevents >= size){
        for (i = 0; i < size; i++){
            eventset[i] = riskset[i];
        }
/* This should _never_ happen: */
        if (antevents > size) 
        error("Error in [sample_events]. Mail a bug report!");
    }

    n = size - 1;

/* Take copies to protect prob and riskset: */
    for (i = 0; i < size; i++){
        perm[i] = riskset[i];
        p[i] = prob[i];
    }

    totalmass = 1.0;

    for(i = 0; i < antevents; i++){
        z = unif_rand();
        rt = z * totalmass;
        mass = 0.0;
        j = 0;
        while ( (mass < rt) & (j < n) ){
            j++;
            mass += p[j];
        }
        eventset[i] = perm[j];
        totalmass = totalmass - p[j];
        for (k = j; k < n; k++){
            p[k] = p[k + 1];
            perm[k] = perm[k + 1];
        }
        n--;
    }
    Free(p);
    Free(perm);
}

/*     Start values:
static void fill_haz(int totrs, RiskSet *risks){
    
    int rs, antevents, size;

    for (rs = 0; rs < totrs; rs++){
	antevents = (risks + rs)->antevents;
	size = (risks + rs)->size;
	if (antevents < size)
	    (risks + rs)->hazard = 
		1.0 - exp(-exp((risks + rs)->gamma));
	else
	    (risks + rs)->hazard = 1.0;
    }
}
*/

static void inv_hess(double *h22, int *fail){

    int i, j;
    char up = 'U';;

/* Invert J, AKA h22: */ 
    
    F77_CALL(dpotrf)(&up, &p, h22, &p, fail);
    if (!(*fail)){
	F77_CALL(dpotri)(&up, &p, h22, &p, fail);
	if (!(*fail)){
	    for ( i = 1; i < p; i++){
		for (j = 0; j < i; j++){
		    h22[i + j*p] = h22[j + i*p];
		}
	    }
	}else{
	    Rprintf("[dpotri] info = %d\n", *fail);
	    error("No inverse in [inv_hess]");
	}
    }else{
	Rprintf("[dpotrf] info = %d\n", *fail);
	error("No inverse in [inv_hess]");
    }

}


static void n_r(int prl, int itmax, int *iter, double eps, 
		int totrs, RiskSet *risks, double e_frac,
		double *b, double *db, 
		double *ll, double *dll, double *d2ll, 
		double *sctest, 
		int *f_conver, int *conver, int *fail){	

    int what;
    int ione = 1;
    double one = 1.0;
    char uplo = 'U';
    int info;
    double L2;
    double ll_prev;

    int i;

    what = 2;

    *iter = 0;
    *conver = 0; 
    *f_conver = 0;
    *fail = 0;

    ll_prev = *ll;
    while ( (*iter < itmax) & (!(*conver)) ){

/* C +++ Do it with Lapack routine 'dposv': */ 

	F77_CALL(dcopy)(&p, dll, &ione, db, &ione);
        F77_CALL(dposv)(&uplo, &p, &ione, d2ll, &p, db, &p, &info);
	if (info){
	    *fail = info;
	    return;
	}
 
/* The score test statistic: */
	if (*iter == 0) *sctest = F77_CALL(ddot)(&p, db, &ione, 
						 dll, &ione);
	

	L2 = F77_CALL(dnrm2)(&p, db, &ione);
	if (L2 < eps) *conver = 1;
	if (prl == 1){
            Rprintf(" \n");
	    Rprintf("*** Iteration %d\n", *iter);
	    Rprintf("L2 = %f\n", L2);
	    Rprintf("loglik = %f\n", *ll);
	    for (i = 0; i < p; i++){
		Rprintf("beta[%d] = %f; dll[%d] = %f\n", 
			i, b[i], i, dll[i]);
	    } 
	}

	F77_CALL(daxpy)(&p, &one, db, &ione, b, &ione);

	coxfun(what, totrs, risks, e_frac, 
	       b, ll, dll, d2ll);

	if (fabs(*ll / ll_prev - one) < eps) *f_conver = 1; 
	ll_prev = *ll;
	(*iter)++;
    }
/*  Done iterating!    */
}

static void fill_in(RiskSet *risks, 
		    int ns,
		    int *antrs,
		    int *antevents,
		    int *eventset,
		    int *size,
		    int *riskset,
		    double *offset,
		    double *weights,
		    double *hazard){

    int str, rs, j, eindx, rindx, i;
    double tmp;

    rs = -1;
    eindx = 0;
    rindx = 0;
    for (str = 0; str < ns; str++){
	for (j = 0; j < antrs[str]; j++){
	    rs++;
	    risks[rs].antevents = antevents[rs];
	    risks[rs].eventset = eventset + eindx;
	    risks[rs].size = size[rs];
	    risks[rs].riskset = riskset + rindx;
	    risks[rs].offset = offset + rindx;
	    risks[rs].weights = weights + rindx;
	    tmp = 0.0;
	    for (i = 0; i < antevents[rs]; i++) 
		tmp += risks[rs].weights[i];
	    risks[rs].rs_weight = tmp / antevents[rs];

	    /* if (risks[rs].antevents == risks[rs].size){ */
	    if (1 == risks[rs].size){
		risks[rs].out = 1;
		risks[rs].hazard = 1.0;
		hazard[rs] = 1.0;
		/* warning("[fill_in] Risk set of size 1."); */ 
	    }else{
		risks[rs].out = 0;
		hazard[rs] = (double)antevents[rs] /
		    (double)size[rs];
		risks[rs].gamma = log(-log1p(-hazard[rs]));
	    }
	    eindx += risks[rs].antevents;
	    rindx += risks[rs].size;
	}
    }
}

/* *** The 'MAIN' subroutine: */

void sup(int *meth, 
	 int *iter, 
	 double *eps, 
	 int *prl, 
	 int *totevent, 
	 int *totrs, 
	 int *ns, 
	 int *antrs, 
	 int *antevents, 
	 int *size,
	 double *weights,
	 int *totsize, 
	 int *eventset, 
	 int *riskset, 
	 int *nn_in, 
	 int *p_in, 
	 double *covar, 
	 double *offset_in,
	 double *startbeta,
	 int *boot,
	 double *efrac,
	 double *beta,
	 double *sd_beta,
	 double *loglik, /* Note: length 2! */ 
	 double *variance, 
	 double *sctest,
	 double *hazard,
	 int *conver, 
	 int *f_conver, 
	 int *fail){


/***************************************************************************
C     meth      : 0 = "efron",
C                 1 = "breslow",
C                 2 = "mppl",
C                 3 = "ml"
C     iter      : On input = maxiter; on output = actual No. of iterations.
C     eps       : Convergence criterion; L2 < eps ===> convergence
C     prl       : Print level; 0 = nothing, 1 = more.
C     totevent  : Total number of events.
C     totrs     : Total number of risk sets.
C     ns        : Number of strata.
C
C     antrs     : antrs(i) = No. of risk sets in stratum i, i = 1, ns.
C     antevents : number of events in each riskset.
C     size      : Size of each risk set.
C
C     totsize   : Sum of the risk set sizes.
C     eventset  : pointers to events in risk sets (length totevent).
C     riskset   : pointers to members of risk sets (length totsize).
C
c     nn        : No. of spells.
C     p         : No. of covariates.
C     covar     : matrix of covariates (nn x p).
C     offset    : Vector of offsets (nn).

C     startbeta : Start values for b (p).
C     boot      : No. of bootstrap samples. 
C     beta      : Vector of coefficients (p x (1 + boot)); return value.
C
C     loglik         : return value.
C     h2 (dll)   : return value.
C     h22 (d2ll) : return value.
C     sctest         : score test statistic, return value.
C
C     hazard : estimated hazard atoms (totrs).
C
C     conver : 1 if convergence, 
C              0 otherwise.
C     fail   : 1 if failure (i.e., linear dependency among covariates
C                                    or singular hessian).
C              0 if success
C +++ 
**************************************************************************/
  
/*  Variables needed for ML: */
    double ll;
    double *db = 0;
    double *b = 0;

    double e_frac = *efrac; /* Choosing between 'efron' & 'ml' */

    int what;

    int itmax, i, j, indx;

    double zero = 0.0;
    double one = 1.0;

    int izero = 0;
    int ione = 1;
    char trans = 'N';
    

    RiskSet *risks = 0;

    double *dll;  /* "dll" */
    double *d2ll; /* "d2ll" */
    int p2; /* = p * p */

/* Variables for boot: */
    int eindx, sindx, rsindx, rs, rep;
    double *prob;
    int *save_eventset = 0;

    double *score = 0; /* for bootstrapping */

/* Correct for C-zero! */

    for (j = 0; j < *totevent; j++) eventset[j] -= 1;
    for (j = 0; j < *totsize; j++) riskset[j] -= 1;

    
/* Fill in globals: */
    p = *p_in;
    nn = *nn_in;
    x = covar;
/*    offset = offset_in; */

    if (p > 0){
	lin = (double *)R_alloc((long)nn, sizeof(double));
	F77_CALL(dcopy)(&nn, &zero, &izero, lin, &ione);

	score = (double *)R_alloc((long)nn, sizeof(double));
	F77_CALL(dcopy)(&nn, &zero, &izero, score, &ione); 

	sumdscore = (double *)R_alloc((long)p, sizeof(double));
	F77_CALL(dcopy)(&p, &zero, &izero, sumdscore, &ione); 
	
	p2 = p * p;
	sumd2score = (double *)R_alloc((long)p2, sizeof(double));
	F77_CALL(dcopy)(&p2, &zero, &izero, sumd2score, &ione); 
    }

/* Fill in risk sets: */    
    risks = (RiskSet *)R_alloc((long)(*totrs), sizeof(RiskSet));

    fill_in(risks, *ns, antrs, antevents, eventset, 
	    size, riskset, offset_in, weights, hazard);

    /* Enough if no covariates! */
    if (p <= 0) {
	*conver = 1;
	*f_conver = 1;
	*fail = 0;
	return;
    }

/***************************/
/*    if (*meth <= 1){
	ml = 0;
	method = *meth;
    }else{
	ml = 1;
	method = *meth - 2;
    }
*/
    switch (*meth){
    case 0:
	eha_rs = &efron_rs;
	break;
    case 1:
	eha_rs = &breslow_rs;
	break;
    case 2: 
	eha_rs = &mppl_rs;
	break;
    case 3:
	eha_rs = &ml_rs;
	break;
    default:
	error ("Wrong method! (Should never happen; file a bug report)");
	break;
    }

    b = (double *)R_alloc((long)(p), sizeof(double));
    db = (double *)R_alloc((long)(p), sizeof(double));
    dll = (double *)R_alloc((long)(p), sizeof(double));
    d2ll = (double *)R_alloc((long)p2, sizeof(double));


    F77_CALL(dcopy)(&p, startbeta, &ione, b, &ione);

    itmax = *iter;

    what = 2;
    coxfun(what, *totrs, risks, e_frac, 
	   b, &ll, dll, d2ll);
      
    loglik[0] = ll;
    loglik[1] = ll;

    n_r(*prl, itmax, iter, *eps, 
	*totrs, risks, e_frac, 
	b, db,
	&ll, dll, d2ll, 
	sctest,
	f_conver, conver, fail);	

    if (*fail) return;
    loglik[1] = ll;
/*                     */
/*  The 'afterwork':   */
         
    F77_CALL(dcopy)(&p, b, &ione, beta, &ione);
/* "First column" of beta is the solution! */

    if (*prl == 1){
	if (*conver){
	    Rprintf("Convergence: Hessian is\n");
	    indx = -1;
	    for (i = 0; i < p; i++){
		Rprintf("\n ");
		for (j = 0; j < p; j++){
		    indx = indx + 1;
		    Rprintf("%f ", d2ll[indx]); 
		}
	    }
	    Rprintf("\n");
	}else
            Rprintf("NOTE: No Convergence!\n");
	Rprintf("loglik = %f\n", ll);
}
      
    loglik[1] = ll;

/* Get the variance(p, p) matrix in h22: */
    
    inv_hess(d2ll, fail);

    if (*prl == 1){
	Rprintf("Variance is:\n");
	indx = -1;
	for (i = 0; i < p; i++){
	    Rprintf("\n ");
	    for (j = 0; j < p; j++){
		indx = indx + 1;
		Rprintf("%f ", d2ll[indx]); 
	    }
	}
	Rprintf("\n");
    }
    if (*fail){
	warning("No variance (this should not happen!)");
	return;
    }else{
	F77_CALL(dcopy)(&p2, d2ll, &ione, variance, &ione);
	/* Fill in 'hazard': */
	indx = -1;
	for (i = 0; i < *ns; i++){
	    for (j = 0; j < antrs[i]; j++){
		indx++;
		get1_gam(risks + indx);
		hazard[indx] = risks[indx].hazard;
	    }
	}
	/* fill_haz(*totrs, risks); */
	/* Get the se(beta): */
	for (j = 0; j < p; j++){
	    sd_beta[j] = sqrt(d2ll[(p + 1) * j]);
	}
}

/* Bootstrapping? */
/* Needs fixing in case of time-varying offset! */
    if (*boot){
	
	prob = (double *)R_alloc((long)(*totsize), sizeof(double));
/*
C +++ Calculate score(i), i = 1, nn:
C     Only needed for calculation of selection probabilities. */

/*	F77_CALL(dcopy)(&nn, offset, &ione, score, &ione); */
	for (j = 0; j < nn; j++) score[j] = 0.0;
	F77_CALL(dgemv)(&trans, &nn, &p, &one, covar, &nn, beta, &ione, &one,  
			score, &ione);
	
	for (j = 0; j < nn; j++)
	    score[j] = exp(score[j]);
	
	save_eventset = (int *)R_alloc((long)(*totevent), sizeof(int));
	/* Save 'eventset' */
	for (j = 0; j < *totevent; j++) save_eventset[j] = eventset[j];

/*
C +++ get the (sorted) sample probabilities.
C     sort 'riskset' accordingly. 
C     Note the consequences in calling function!!
*/
	sindx = 0;
	rsindx = -1;
	for (rs = 0; rs < *ns; rs++){
	    for (j = 0; j < antrs[rs]; j++){
		rsindx++;
		get_prob(nn, size[rsindx], (riskset + sindx), 
			 score, (prob + sindx));
		sindx = sindx + size[rsindx];
	    }
	}
	GetRNGstate();
	
	for (rep = 0; rep < *boot; rep++){
	    if (*prl == 1)
		Rprintf("rep = %d\n", rep);
	    F77_CALL(dcopy)(&p, beta, &ione, b, &ione);
	    sindx = 0;
	    eindx = 0;
	    rsindx = -1;
/* Sample riskset, put in eventset:  */
	    for (rs = 0; rs < *ns; rs++){
		for (j = 0; j < antrs[rs]; j++){
		    rsindx++;
		    sample_events(antevents[rsindx], size[rsindx],
				  (riskset + sindx), (eventset + eindx), 
				  (prob + sindx));
		    sindx += size[rsindx];
		    eindx += antevents[rsindx];
		}
	    }
	    /* Get bootstrap sample: */

	    n_r(*prl, itmax, iter, *eps,
		*totrs, risks, e_frac,
		b, db, 
		&ll, dll, d2ll,
		sctest, f_conver, conver, fail);	
	
	    
	    inv_hess(d2ll, fail);

	    F77_CALL(dcopy)(&p, b, &ione, (beta + (rep + 1) * p), &ione);
	    for (j = 0; j < p; j++){
		sd_beta[(rep + 1) * p + j] = sqrt(d2ll[(p + 1) * j]);
	    }
	} /* End big loop */
	for (j = 0; j < *totevent; j++) eventset[j] = save_eventset[j];
	PutRNGstate();
    }
    /* Restore 'eventset', 'riskset': */
    for (j = 0; j < *totevent; j++) eventset[j]++;
    for (j = 0; j < *totsize; j++) riskset[j]++;

}
