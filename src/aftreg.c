#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>

#include "aftreg.h"
#include "phfun.h"

/*    
ph0S_fun *S0;  
ph0_fun *f0;
ph0_fun *h0;   
ph0_fun *f0_t; 
ph0_fun *h0_t; 
ph0_fun *h0_tt;
*/


static double aft_fun(int n, double *beta, void *vex){
/* Here shape is NOT fixed, i.e. estimated, and part of beta! */
    int i, j, rec;
    Exts *ex;
    double alpha, gamma;
    int nn, indiv, stratum, mb;
    double res1, res2;
    double alphambz;
    double *bz;
    int *n_rec;
    double a_time, b_time;
    int log_p = 1;

    ex = vex;
    mb = *(ex->mb);
    nn = *(ex->nn);

    bz = Calloc(nn, double);

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
    }else if (dist == 3){
	S0 = &S0_ev;
	f0 = &f0_ev;      
	h0 = &h0_ev;      
	f0_t = &f0_t_ev;     
	h0_t = &h0_t_ev;     
	h0_tt = &h0_tt_ev;
    }else if (dist == 4){ /* 4 = Gompertz */    
	error("Gompertz should not go here (internal error)");
    }else{
	error("Unknown distribution");
    }

    indiv = 1;
    for (i = 1; i < nn; i++){
	if (ex->id[i] != ex->id[i-1]) indiv++;
    }
    n_rec = Calloc(indiv, int);
    for (i = 0; i < indiv; i++) n_rec[i] = 1;
    j = 0;
    for (i = 1; i < nn; i++){
	if (ex->id[i] == ex->id[i-1]){
	    n_rec[j]++;
	}else{
	    j++;
	}
    }

    res1 = 0.0;
    res2 = 0.0;

    for (i = 0; i < nn; i++) bz[i] = ex->offset[i];

    if (mb) {
	for (i = 0; i < nn; i++) {
	    /* Kolla detta: Ska senare bytas mot BLAS */  
	    for (j = 0; j < mb; j++) bz[i] += ex->z[i * mb + j] * beta[j];
	}
    }

    rec = 0;

    for (i = 0; i < indiv; i++){
	stratum = ex->strata[rec];
	alpha = beta[mb + 2 * stratum];          /*     Log scale!! */
	gamma = exp(beta[mb + 2 * stratum + 1]); /* NOT Log scale!! */
	alphambz = alpha - bz[rec];

	a_time = ex->time0[rec] * exp(-alphambz); /* Without */
	b_time = ex->time[rec] * exp(-alphambz); /* shape!  */
	if (ex->ind[rec]){ 
	    res1 += log(gamma) - alphambz + 
		(gamma - 1) * (log(ex->time[rec]) - alphambz) +
		log(h0(R_pow(b_time, gamma)));
	}
	res2 += S0(R_pow(a_time, gamma), log_p) - 
	    S0(R_pow(b_time, gamma), log_p);

	if (n_rec[i] >= 2){
	    for (j = 1; j < n_rec[i]; j++){
		rec++;
		stratum = ex->strata[rec];
		alpha = beta[mb + 2 * stratum];           /* See    */
		gamma = exp(beta[mb + 2 * stratum + 1]);  /* above! */
		alphambz = alpha - bz[rec];
		a_time = b_time;
		b_time = a_time + 
		    (ex->time[rec] - ex->time0[rec]) * exp(-alphambz);
		if (ex->ind[rec]){ 
		    res1 += log(gamma) - alphambz + 
			(gamma - 1) * (log(ex->time[rec]) - alphambz) +
			log(h0(R_pow(b_time, gamma)));
		}
		res2 += S0(R_pow(a_time, gamma), log_p) - 
		    S0(R_pow(b_time, gamma), log_p);
	    }
	}
	rec++;
    }

    Free(n_rec);
    Free(bz);

    return( -(res1 - res2) ); /* Minimizing ... */
}
	
void aftsup(int *printlevel,
	    int *ns, int *nn, int *ncov, int *bdim,
	    int *id, int *strata, double *time0, double *time, int *ind,
	    double *covar, double *offset, int *dis, double *beta,
	    double *loglik, int *fail){
/* Here shape is NOT fixed, i.e. estimated, and part of beta!           */
/* *ncov = No. of covariates                                            */
/* *bdim = No. of parameters to setimate                                */
/* *ns   = No. of strata, i.e, *bdim = *ncov + 2 * ns                   */
/* beta  = parameter vector, beta[0:(*mb-1)] for covariates,            */
/*                          beta[*mb, *mb+2, ...] are scale parameters  */
/*                          beta[*mb+1, *mb+3, ...] are shape parameter */
    Exts *ex;
    /* int iok; */

    void *vex;

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
    }else if ((dist == 3) || (dist == 4)){ /* 4 = Gompertz */    
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

    /* iok = 0; */

/* Fill in 'ex': */
    ex->id = id;
    ex->strata = strata;
    ex->ns = ns;
    ex->pfix = NULL; /* Only in aftexpreg, fixed shape */
    ex->mb = ncov;
    ex->nn = nn;
    ex->z = covar;
    ex->time0 = time0;
    ex->time = time;
    ex->ind = ind;
    ex->offset = offset;

/* Done with filling in */

    *loglik = aft_fun(*bdim, beta, vex);
}
