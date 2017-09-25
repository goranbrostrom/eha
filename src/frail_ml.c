#include <stdio.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>

#include "frail_ml.h"
#include "eha_fun.h"
#include "ghq.h"

eha_P_fun *P;
eha_G_fun *G;
eha_Gprim_fun *Gprim;

void frail_ml(int *family,
	      int *method,
	      int *p,
	      int *nn,
	      int *n_rs,
	      int *riskset,
	      double *start_beta,
	      double *start_sigma,
	      double *x, /* Now p x (\sum_1_{n_fam} fam_size[i]) */
	      int *y, /* Looong vector! */
	      int * haz, /* Looong vector! */
	      double *offset,
	      int *fam_size,
	      int *n_fam,
	      int *n_points, /* No. of pts in Gauss-Hermite quadrature */
	      double *epsilon,
	      int *maxit,
	      int *trace,
	      double *beta,
	      double *hazards,
	      double *sigma,
	      double *sigma_sd,
	      double *loglik,
	      double *variance,
	      double *frail,
/*	      double *mu,        */
	      int *convergence,
	      int *fail){
    
#ifndef MATHLIB_STANDALONE
    double abstol;
    double reltol;
    int nREPORT = 1;
    int fncount;
    int grcount;
    int *mask;
#endif
    Exts *ext;
    int i, j;
    
    double Fmin;
    double *b;
    double *gr;
    int bdim;

    int job = 11;
    double *work;
    double *det;
    int lwork;
    double rcond;
    int info;
    int vm_trace;
    int true_bdim;
    int modified = 0; /* Check this!! */

    vm_trace = *trace;
    bdim = *p + *n_rs + 1;

    det = Calloc(2, double);
    lwork = 11 * bdim;
    work = Calloc(lwork, double);

    if (*family == 0){
	P = &eha_P_logit;
	G = &eha_G_logit;
	Gprim = &eha_Gprim_logit;
    }else if (*family == 1){
	P = &eha_P_cloglog;
	G = &eha_G_cloglog;
	Gprim = &eha_Gprim_cloglog;
    }else if (*family == 2){
	P = &eha_P_poisson;
	G = &eha_G_poisson;
	Gprim = &eha_Gprim_poisson;
    }else{
	error("Unknown family\n");
    }
	
    abstol = *epsilon;
    reltol = abstol;

    b = Calloc(bdim, double);
    gr = Calloc(bdim, double);

    ext = Calloc(1, Exts);

    ext->family = *family; /* == 0 for binomial(logit) */

    for (i = 0; i < *n_rs + *p; i++){
	b[i] = start_beta[i];
    }

    if (*n_points == 1){/* Put in log(sigma) */
	b[bdim - 1] = 0.0;
    }else{
	b[bdim - 1] = *start_sigma;
    }
    ext->p = *p;
    ext->n_rs = *n_rs;
    ext->riskset = riskset;
    ext->haz = haz;
    
    ext->offset = offset;
    ext->weights = Calloc(*n_points, double);
    ext->zeros = Calloc(*n_points, double);

    mask = Calloc(bdim, int );

    F77_CALL(ghq)(n_points, ext->zeros, ext->weights, &modified); 

    ext->n = 0;
    for (i = 0; i < *n_fam; i++){
	ext->n += fam_size[i];
    }


    ext->x_beta = Calloc(ext->n, double); 

    ext->gr = gr;
    ext->hessian = variance;
    ext->x = x;
    ext->y = y;
    ext->n_fam = *n_fam;
    ext->fam_size = fam_size;
    ext->n_points = *n_points;

/* First, maximize with regression coefficients fixed at start values: */

    for (i = 0; i < *n_rs; i++) mask[i] = 1;
    for (i = *n_rs; i < *n_rs + *p; i++) mask[i] = 0;

    if (ext->n_points == 1){
	mask[bdim - 1] = 0;
	b[bdim - 1] = 0.0;
    }else{
	mask[bdim - 1] = 1;
    }

 
    *maxit = 100;
    
    vmmin(bdim, b, &Fmin,
	  eha_fun, eha_fun1, *maxit, vm_trace,
	  mask, abstol, reltol, nREPORT,
	      ext, &fncount, &grcount, convergence);
 
    loglik[0] = -Fmin;
    eha_fun2(bdim, b, &Fmin, gr, variance, ext);
    if(*trace){
	Rprintf("Max log likelihood after restricted [vmmin]: %f\n", -Fmin);
	Rprintf("beta: ");
	for (i = *n_rs; i < *n_rs + *p + 1; i++){
	    Rprintf(" %f, ", b[i]);
	}
	Rprintf("\n");
	Rprintf("Gradients: ");
	for (i = *n_rs; i < *n_rs + *p + 1; i++){
	    Rprintf(" %f, ", -ext->gr[i]);
	}
	Rprintf("\n");

	Rprintf("\nhessian:");
	for (i = *n_rs; i < *n_rs + *p + 1; i++){
	    Rprintf("\n");
	    for (j = *n_rs; j < *n_rs + *p + 1; j++){
		Rprintf(" %f ", variance[i + j * bdim]);
	    }
	}
	Rprintf("\n");
    }

    /* Then the full maximization: */
    for (i = 0; i < bdim; i++) mask[i] = 1;
    if (ext->n_points == 1){
	mask[bdim - 1] = 0;
	b[bdim - 1] = 0.0;
    }
    vmmin(bdim, b, &Fmin,
	  eha_fun, eha_fun1, *maxit, vm_trace,
	  mask, abstol, reltol, nREPORT,
	      ext, &fncount, &grcount, convergence);
 
    loglik[1] = -Fmin;

    eha_fun2(bdim, b, &Fmin, gr, variance, ext);
    
    if(*trace){
	Rprintf("Max log likelihood after [vmmin]: %f\n", -Fmin);
	Rprintf("beta: ");
	for (i = *n_rs; i < *n_rs + *p + 1; i++){
	    Rprintf(" %f, ", b[i]);
	}
	Rprintf("\n");
	Rprintf("Gradients: ");
	for (i = *n_rs; i < *n_rs + *p + 1; i++){
	    Rprintf(" %f, ", -ext->gr[i]);
	}
	Rprintf("\n");

	Rprintf("\nhessian:");
	for (i = *n_rs; i < *n_rs + *p + 1; i++){
	    Rprintf("\n");
	    for (j = *n_rs; j < *n_rs + *p + 1; j++){
		Rprintf(" %f ", variance[i + j * bdim]);
	    }
	}
	Rprintf("\n");
    }
/*
    eha_nr_opt(bdim, b, &Fmin, mask, ext, *epsilon, *maxit, *trace);
    loglik[1] = Fmin;
*/

    F77_CALL(dpoco)(variance, &bdim, &bdim, &rcond, work, &info);

    if (info == 0){
	F77_CALL(dpodi)(variance, &bdim, &bdim, det, &job);
	
	for (i = 0; i < bdim; i++){
	    for (j = 0; j < i; j++){
		variance[i + j * bdim] = variance[j + i * bdim];
	    }
	}
    }else{
	if (info == bdim){
	    eha_fun2(bdim, b, &Fmin, gr, variance, ext);
	    true_bdim = bdim - 1;
	    F77_CALL(dpoco)(variance, &bdim, &true_bdim, &rcond, work, &info);
	    if (info == 0){
		F77_CALL(dpodi)(variance, &bdim, &true_bdim, det, &job);
		
		for (i = 0; i < bdim; i++){
		    for (j = 0; j < i; j++){
			variance[i + j * bdim] = variance[j + i * bdim];
		    }
		}
	    }else{
		Rprintf("info from [dpoco] = %d\n", info); 
		Rprintf("No inversion in [frail_ml]\n");
	    }
	}else{
	    Rprintf("info from [dpoco] = %d\n", info); 
	    Rprintf("No inversion in [frail_ml]\n");
	}
    }

    for (i = 0; i < *n_rs; i++){
	hazards[i] = b[i];
    }
    for (i = 0; i < *p; i++){
	beta[i] = b[i + *n_rs];
    }
    *sigma = fabs(b[*p + *n_rs]);
    *sigma_sd = sqrt(variance[bdim * bdim - 1]);
    if(*trace){
	Rprintf("Max log likelihood after [eha_nr_opt]: %f\n", -Fmin);
	Rprintf("beta: ");
	for (i = *n_rs; i < bdim; i++){
	    Rprintf(" %f, ", b[i]);
	}
	Rprintf("\n");
	Rprintf("Gradients: ");
	for (i = *n_rs; i < bdim; i++){
	    Rprintf(" %f, ", -ext->gr[i]);
	}
	Rprintf("\n");

	Rprintf("\nVariance:");
	for (i = *n_rs; i < bdim; i++){
	    Rprintf("\n");
	    for (j = *n_rs; j < bdim; j++){
		Rprintf(" %f ", variance[i + j * bdim]);
	    }
	}
	Rprintf("\n");
    }

    eha_frail_fun(bdim, b, frail, ext);
/*    eha_mu_fun(bdim, b, mu, ext); */

    
    Free(ext->x_beta);
    Free(mask);
    Free(ext->zeros);
    Free(ext->weights);
    Free(ext);
    Free(gr);
    Free(b);
    Free(work);
    Free(det);
 }
