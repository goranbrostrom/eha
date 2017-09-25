#include <stdio.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

#include "glmmml.h"
#include "fun.h"
#include "ghq.h"

int laplace; /* 1 = Laplace, 0 = Gauss-Hermite */
P_fun *P;
G_fun *G;
H_fun *H;
I_fun *I;
K_fun *K;
logpr *logprior;
d_logpr *d_logprior;
d2_logpr *d2_logprior;
d3_logpr *d3_logprior;
d4_logpr *d4_logprior;

void glmm_ml(int *family,
	     int *p, 
	     double *start_beta,
	     int *cluster,
	     double *weights,
	     double *cluster_weights,
	     double *start_sigma,
	     int *fix_sigma,
	     double *x, /* Now p x (\sum_1_{n_fam} fam_size[i]) */
	     double *y, /* NOTE! Change from 'int' !! */
	     double *offset,
	     int *fam_size,
	     int *n_fam,
	     int *method,  /* = method, global variable */
	     int *n_points, /* No. of pts in Gauss-Hermite quadrature */
	     double *epsilon,
	     int *maxit,
	     int *trace,
	     int *boot,
	     int *prior, /* 0 = Normal, 1 = logistic */
	     double *predicted,
	     double *beta,
	     double *sigma,
	     double *loglik,
	     double *variance,
	     double *post_mode,
	     double *post_mean,
	     double *mu,
	     double *boot_p,
	     double *boot_log,
	     int *convergence,
	     int *info){
    
#ifndef MATHLIB_STANDALONE
    double abstol;
    double reltol;
    int nREPORT = 1;
    int fncount;
    int grcount;
    int fail;
    int *mask;
    double *lower, *upper;
    int * nbd;
    char msg[60];
#endif
    Exts *ext;
    int i, j, k, m, indx;
    
    double Fmin;
    double *b;
    double *gr;
    double **hessian;
    double *hess_vec;

    double rcond;
    double *det;
    int lwork;
    double *work;
    int job = 11;
    char *vmax;

    int bdim;

/* New in 0.28; bootstrapping: */
    int upp;

    int modified = 1;

    /* int *conditional; */
    int condi = 1;

    laplace = *method; /* 1 = Laplace, 0 = Gauss-Hermite */
/* This is done to prepare for having conditional as an input parameter */     
    /* conditional = &condi; */

    if (*trace){
	Rprintf("************* Entering [glmmml] **************** \n");
	Rprintf(" p = %d\n\n", *p);
    }

    det = Calloc(2, double);

    if (*family == 0){
	P = &P_logit;
	G = &G_logit;
	H = &H_logit;
	I = &I_logit;
	K = &K_logit;
    }else if (*family == 1){
	P = &P_cloglog;
	G = &G_cloglog;
	H = &H_cloglog;
	I = &I_cloglog;
	K = &K_cloglog;
    }else if (*family == 2){
	P = &P_poisson;
	G = &G_poisson;
	H = &H_poisson;
	I = &I_poisson;
	K = &K_poisson;
    }else{
	error("Unknown family\n");
    }

    if (*prior == 1){
	logprior = &logprior_logistic;
	d_logprior = &d_logprior_logistic;
	d2_logprior = &d2_logprior_logistic;
	d3_logprior = &d3_logprior_logistic;
	d4_logprior = &d4_logprior_logistic;
    }else if (*prior == 2){
	logprior = &logprior_cauchy;
	d_logprior = &d_logprior_cauchy;
	d2_logprior = &d2_logprior_cauchy;
	d3_logprior = &d3_logprior_cauchy;
	d4_logprior = &d4_logprior_cauchy;
    }else{	
	logprior = &logprior_normal;
	d_logprior = &d_logprior_normal;
	d2_logprior = &d2_logprior_normal;
	d3_logprior = &d3_logprior_normal;
	d4_logprior = &d4_logprior_normal;
    }

    abstol = *epsilon;
    reltol = abstol;

    bdim = *p + 1;
    lwork = 11 * bdim;
    work = Calloc(lwork, double);

    lower = Calloc(bdim, double);
    upper = Calloc(bdim, double);
    nbd = Calloc(bdim, int);
    for (j = 0; j < bdim; j++){
	nbd[j] = 0;
	upper[j] = 0.0;
	lower[j] = 0.0;
    }
    nbd[bdim - 1] = 1;
    lower[bdim - 1] = 0.5e-10;
    b = Calloc(bdim, double);
    gr = Calloc(bdim, double);

    hessian = Calloc(bdim, double *);
    hess_vec = Calloc(bdim * bdim, double);
    for (j = 0; j < bdim; j++) hessian[j] = hess_vec + j * bdim;

    /**** Build up 'ext' ********************/
    ext = Calloc(1, Exts);

    ext->family = *family; /* == 0 for binomial(logit) */

    ext->n = 0;
    for (i = 0; i < *n_fam; i++){
	ext->n += fam_size[i];
    }
    ext->p = *p;
    ext->cluster = cluster;
    /* Changed 2006-06-18; may have catastrophic consequences!! */
    ext->x = Calloc(ext->p, double *);
    for (i = 0; i < ext->p; i++){
	ext->x[i] = x + i * (ext->n);
    }
    /*** Note that ext->x is not "filled"; ***/ 
    /*** only points to the right place    ***/
    ext->offset = offset;
    ext->x_beta = Calloc(ext->n, double);
    ext->yw = Calloc(ext->n, double); /* We cannot copy if bootstrapping! */
    for (i = 0; i < ext->n; i++){
	ext->yw[i] = y[i] * weights[i]; /* NOTE !!! */
    }
    /* error("Enough!!!"); */
    ext->weights = weights;
    ext->cluster_weights = cluster_weights;
    ext->n_fam = *n_fam;
    ext->fam_size = fam_size;
    ext->post_mode = Calloc(*n_fam, double); 
    ext->post_mean = Calloc(*n_fam, double); 
    ext->n_points = *n_points;
    ext->wc = Calloc(*n_points, double);
    ext->zeros = Calloc(*n_points, double);
    F77_CALL(ghq)(n_points, ext->zeros, ext->wc, &modified);  

/******* Done with 'ext' ***************/

    for (j = 0; i < *p; j++){
	b[j] = start_beta[j];
    }
/* NOTE; here we do not log sigma!!! */
    b[*p] = *start_sigma;

    mask = Calloc(bdim, int );
    for (i = 0; i < bdim; i++){
        mask[i] = 1;
    }
    if (*fix_sigma) 
	mask[bdim - 1] = 0; /* sigma fixed during optimization */

/*  Note that this performs a minimum: (!!) */
/*
    for (i = 0; i < bdim; i++) Rprintf("b[%d] = %f\n", i, b[i]);
    fun2(bdim, b, &Fmin, gr, hess_vec, ext);
    Rprintf("Fmin = %f\n", Fmin);
    for (i = 0; i < bdim; i++) Rprintf("gr[%d] = %f\n", i, gr[i]);
    for (i = 0; i < bdim * bdim; i++) Rprintf("hess[%d] = %f\n", 
						  i, hess_vec[i]);
*/  
/*    error("Let us stop here"); */

	vmax = vmaxget();
	
	vmmin(bdim, b, &Fmin,
	      fun, fun1, *maxit, *trace,
	      mask, abstol, reltol, nREPORT,
	      ext, &fncount, &grcount, &fail);
	vmaxset(vmax);
/*    }else{
	vmax = vmaxget();
	lbfgsb(bdim, lmm, b, lower, upper, nbd, &Fmin,
	      fun, fun1, &fail, ext, factr, pgtol, 
	       &fncount, &grcount, *maxit, msg, *trace,
	      nREPORT);
	vmaxset(vmax);
	} */
    *convergence = (fail == 0);
    if (fail){
	Rprintf("[glmmml] fail = %d\n", fail);
	if (fail == 1) 
	    Rprintf("Max. No. of iterations reached without convergence");
	warning(msg);
    }
/*
    fun1(bdim, b, gr, ext);
*/
    fun2(bdim, b, &Fmin, gr, hess_vec, ext);

    if (*trace){
	Rprintf("Max log likelihood after vmmin: %f\n", -Fmin);
	Rprintf("beta: ");
	for (i = 0; i < bdim; i++){
	    Rprintf(" %f, ", b[i]);
	}
	Rprintf("\n");
	Rprintf("Gradients: ");
	for (i = 0; i < bdim; i++){
	    Rprintf(" %f, ", -gr[i]);
	}
	Rprintf("\n");
	Rprintf("\n");
	Rprintf("hessian:\n");
	for (i = 0; i < bdim; i++){
	    for (j = 0; j < bdim; j++)
		Rprintf(" %f, ", hessian[i][j]);
	    Rprintf("\n");
	}
    }
    /* Let's avoid nr_opt! Just calculate the hessian and invert! */
/*
  nr_opt(bdim, b, &Fmin, mask, ext, *epsilon, nr_maxit, info, *trace); 
*/
    *loglik = Fmin;
    for (i = 0; i < *p; i++){
	beta[i] = b[i];
    }
    if (!(*fix_sigma)){
	*sigma = b[*p]; /* NOTE!!!!!!!! */
	F77_CALL(dpoco)(*hessian, &bdim, &bdim, &rcond, work, info);
	if (*info == 0){
	    F77_CALL(dpodi)(*hessian, &bdim, &bdim, det, &job);
	    for (m = 0; m < bdim; m++){
		for (k = 0; k < m; k++){
		    hessian[k][m] = hessian[m][k];
		}
	    }
	    
	}else{
	    Rprintf("info = %d\n", *info);
	    warning("Hessian non-positive definite. No variance!");
	}
    }
    
    indx = 0;
    for (m = 0; m < bdim; m++){ /* Hessian if fix_sigma */
	for (k = 0; k < bdim; k++){
	    variance[indx] = hessian[m][k];
	    indx++;
	}
    }
    
    if(*trace){
	Rprintf("Max log likelihood: %f\n", -Fmin);
	Rprintf("Beta: ");
	for (i = 0; i < bdim; i++){
	    Rprintf(" %f, ", b[i]);
	}
	Rprintf("\n");
	Rprintf("Gradients: ");
	for (i = 0; i < bdim; i++){
	    Rprintf(" %f, ", gr[i]);
	}
	if (!(*fix_sigma)){
	    Rprintf("\n");
	    Rprintf("variance:\n");
	    for (i = 0; i < bdim; i++){
		for (j = 0; j < bdim; j++)
		    Rprintf(" %f, ", hessian[i][j]);
		Rprintf("\n");
	    }
	}
	if (laplace == 1){
	    Rprintf("\nMethod is Laplace\n");
	}else{
	    if (laplace == 0){
		Rprintf("Method is Gauss-Hermite\n");
	    }else{
		Rprintf("Method is %d\n", laplace);
	    }
	}
    }

/* Cancelled for the time being...
    frail_fun(bdim, b, ext);
    mu_fun(bdim, b, mu, ext);
*/
    for (i = 0; i < ext->n_fam; i++){
	post_mode[i] = ext->post_mode[i];
	/* post_mean[i] = ext->post_mean[i]; */
    }

    if ((*boot > 0) & !(*fix_sigma)){
/************** Bootstrapping starts *****************************/
	upp = 0;
	GetRNGstate();
	for (i = 0; i < *boot; i++){
	    if (*trace){
		if ((i / 10) * 10 == i)
		    Rprintf("********************* Replicate No. %d\n", i);
	    }
	 
	    if (*family <= 1){ /* Bernoulli */
		for (j = 0; j < ext->n; j++)
		    ext->yw[j] = rbinom((int)weights[j], predicted[j]);
	    }else{
		for (j = 0; j < ext->n; j++) /* Poisson */
		    ext->yw[j] = rpois(weights[j] * predicted[j]);
	    }
	
/* Restore beta as start values: */
	    for ( j = 0; j < *p; j++) b[j] = beta[j];
	    if (*trace){
/*
		Rprintf("Sampled values by cluster:\n");
		for (i = 0; i < ext->n; i++){
		    Rprintf("y[%d] = %d, cluster[%d] = %d\n", 
			    i, ext->y[i], i, ext->cluster[i]);
		}
*/
		Rprintf("Start value to vmmin: %f\n", fun(*p, b, ext)); 

	    }
	    /*    if (*method){ */
		vmax = vmaxget();

		vmmin(*p, b, &Fmin,
		      fun, fun1, *maxit, *trace,
		      mask, abstol, reltol, nREPORT,
		      ext, &fncount, &grcount, &fail);
		vmaxset(vmax);
/* Skip this for the moment! (at least!!)  }else{
		vmax = vmaxget();
		lbfgsb(bdim, lmm, b, lower, upper, nbd, &Fmin,
		       fun, fun1, &fail, ext, factr, pgtol, 
		       &fncount, &grcount, *maxit, msg, *trace,
		       nREPORT);
		
		vmaxset(vmax);
		}*/
	    *convergence = (fail == 0);
	    if (!convergence){
		Rprintf("No convergence...\n");
	    }
	    boot_log[i] = -Fmin;
	    if (*trace){
		Rprintf("boot_log[%d] = %f; loglik = %f\n", i, -Fmin, *loglik);
		Rprintf("beta[0] = %f\n", b[0]);
	    }
	    if (-Fmin >= *loglik) upp++;
	}
	
	if (*boot) *boot_p = (double)upp / (double)*boot;
	else *boot_p = 1.0;
	
	PutRNGstate();
	
    }
    
    Free(mask);

    Free(ext->zeros);
    Free(ext->wc);

    Free(ext->post_mean);
    Free(ext->post_mode);
    Free(ext->x_beta);
    Free(ext->yw);
    Free(ext->x);
    Free(ext);

    Free(hessian);
    Free(hess_vec);
    Free(gr);
    Free(b);
    Free(upper);
    Free(lower);
    Free(nbd);
    Free(work);
    Free(det);
}
