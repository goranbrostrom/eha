#include <stdio.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

#include "glmmboot.h"
#include "bfun.h"
#include "fun.h"

extern P_fun *P;
extern G_fun *G;
extern H_fun *H;

void permute(int n, int *y, int *x)
{
/* Result in y; x is just a "work area" */

    int i, j, k;

    k = n; /* Eg, "wt sampling k-of-n" */

    for (i = 0; i < n; i++)
	x[i] = i;
    for (i = 0; i < k; i++) {
	j = n * unif_rand();
/*	y[i] = x[j] + 1; */
	y[i] = x[j];
	x[j] = x[--n];
    }
}

void glmm_boot(int *family,
	       int *p, 
	       double *start_beta,
	       int *cluster,
	       double *weights,
	       double *x, /* Now p x (\sum_1_{n_fam} fam_size[i]) */
	       double *y,
	       double *offset,
	       int *fam_size,
	       int *n_fam,
	       double *epsilon,
	       int *maxit,
	       int *trace,
	       int *boot,
	       double *beta,
	       double *predicted,
	       double *fitted,
	       double *loglik,
	       double *variance,
	       int *info, 
	       double *frail,
	       double *boot_p,
	       double *boot_log,
	       int *convergence){
    
#ifndef MATHLIB_STANDALONE
    double abstol;
    double reltol;
    int nREPORT = 1;
    int fncount;
    int grcount;
    int fail;
    int *mask;
#endif
    Extb *ext;
    Cluster *clust;
    int i;
    int j;
    int cl;
    int indx;

    double Fmin;
    double *b;
    double *gr;
    int upper;
    char *vmax;
    int bdim, m, k;

    double ** hessian;
    double *hess_vec;

    double *det;
    int lwork;
    double *work;
    double rcond;
    int job = 11;
    int ant_fam_out;

    double tmp;
    double one = 1.0;

    int rep;

    bdim = *p;
    lwork = 11 * (*p);
    work = Calloc(lwork, double);
    det = Calloc(2, double);

    gr = Calloc(bdim, double);
    hessian = Calloc(bdim, double *);
    hess_vec = Calloc(bdim * bdim, double);
    for (j = 0; j < bdim; j++) hessian[j] = hess_vec + j * bdim;

    GetRNGstate(); /* For random number generation */

/*    vmax1 = vmaxget(); */

    if (*family == 0){
	P = &P_logit;
	G = &G_logit;
	H = &H_logit;
    }else if (*family == 1){
	P = &P_cloglog;
	G = &G_cloglog;
	H = &H_cloglog;
    }else if (*family == 2){
	P = &P_poisson;
	G = &G_poisson;
	H = &H_poisson;
    }else{
	error("Unknown family\n");
    }

    abstol = 0.00000001;
    reltol = abstol;

    ext = Calloc(1, Extb);
    clust = Calloc(*n_fam, Cluster);
/************************ Fill in ext: *****************/
    ext->family = *family; /* == 0 for binomial(logit) */

    ext->n = 0;
    for (i = 0; i < *n_fam; i++){
	ext->n += fam_size[i];
    }
    ext->p = *p;
    ext->n_clust = *n_fam;

    ext->clust = clust;

    indx = 0;
    for (cl = 0; cl < ext->n_clust; cl++){
	clust[cl].n = fam_size[cl];
	clust[cl].p = ext->p;
	clust[cl].yw = Calloc(clust[cl].n, double);
	clust[cl].lin = Calloc(clust[cl].n, double);
	clust[cl].weight = weights + indx;
	clust[cl].offset = offset + indx;
	clust[cl].x = Calloc(clust[cl].n, double *); /* KOLLA!!! */
	for (i = 0; i < clust[cl].n; i++){
	    clust[cl].x[i] = x + indx * (ext->p);
	    clust[cl].yw[i] = weights[indx] * y[indx];
	    indx++;
	}
    }
    for (cl = 0; cl < ext->n_clust; cl++){ 
	clust[cl].ytot = 0.0;
	clust[cl].wtot = 0.0;
	for (i = 0; i < clust[cl].n; i++){
	    clust[cl].wtot += clust[cl].weight[i];
	    clust[cl].ytot += clust[cl].yw[i];
	}
    }

    ant_fam_out = 0;
    
    for (i = 0; i < ext->n_clust; i++){
	if (fabs(clust[i].ytot) < 0.001){
	    clust[i].out = -1;
	    clust[i].gamma = -1000.0; /* -inf */
	}else if ((fabs(clust[i].wtot - clust[i].ytot) < 0.001) & 
		  (ext->family <= 1)){
	    clust[i].out = 1;
	    clust[i].gamma = 1000.0; /* +inf */
	}else{
	    ant_fam_out++;
	    clust[i].out = 0;
	}
    }
/**************** Filled in ext  *************************/

    if (!ant_fam_out) 
	error("All clusters are 'trivial' (all zeros or all ones)");

    mask = Calloc(ext->p, int);    

    b = Calloc(ext->p, double);



    for (i = 0; i < *p; i++){
	b[i] = start_beta[i];
	if (*trace) Rprintf("start_beta[%d] = %f", i, b[i]);
    }

    for (i = 0; i < ext->p; i++){
        mask[i] = 1;
    }

/* Note that this searches for a minimum: (!!) */
    vmax = vmaxget();

    vmmin(*p, b, &Fmin,
	  bfun, bfun_gr, *maxit, *trace,
	  mask, abstol, reltol, nREPORT,
	  ext, &fncount, &grcount, &fail);
    *convergence = (fail == 0);
    vmaxset(vmax);
    bfun_gr(*p, b, gr, ext);
    if(*trace){
	Rprintf("Max log likelihood after vmmin: %f\n", -Fmin);
	Rprintf("Gradients: ");
	for (i = 0; i < *p; i++){
	    Rprintf(" %f, ", gr[i]);
	}
	Rprintf("\n");
    }
    
    *loglik = -Fmin;
    for (i = 0; i < *p; i++){
	beta[i] = b[i];
    }
    for (i = 0; i < ext->n_clust; i++){
	frail[i] = clust[i].gamma;
    }

/* Calculate fitted values, _including_ the clustering... */

    indx = -1;
    for (i = 0; i < *n_fam; i++){
	if (clust[i].out){
	    for (j = 0; j < clust[i].n; j++){
		indx++;
		fitted[indx] = clust[i].yw[j];
	    }
	}else{
	    for (j = 0; j < clust[i].n; j++){
		indx++;
		if (ext->family <= 1){
		    fitted[indx] = exp(P(clust[i].lin[j] + 
					 clust[i].gamma, 
					 one, 
					 clust[i].weight[j]));
		}else{
		    fitted[indx] = exp(clust[i].lin[j] +
				       clust[i].gamma);
		}
	    }
	}
    }

/* Done in calling R function.... 
    if (*family <= 1){
	for (j = 0; j < ext->n; j++)
	    predicted[j] = P(ext->x_beta[j], 1);

    }else{
	for (j = 0; j < ext->n; j++)
	    predicted[j] = exp(ext->x_beta[j]);
    }
*/
    hess_vec[0] = 0.0;

    bfun_hess(*p, beta, hess_vec, ext);

    if (*trace){
	Rprintf("Hessian...\n\n");
	for (i = 0; i < *p; i++){
	    for (j = 0; j < *p; j++){
		Rprintf("%f  ", hessian[i][j]);
	    }
	    Rprintf("\n");
	}
    }

    F77_CALL(dpoco)(*hessian, &bdim, &bdim, &rcond, work, info);
    if (*info == 0){
	F77_CALL(dpodi)(*hessian, p, p, det, &job);
	for (m = 0; m < bdim; m++){
	    for (k = 0; k < m; k++){
		hessian[k][m] = hessian[m][k];
	    }
	}
	for (m = 0; m < bdim * bdim; m++) variance[m] = hess_vec[m];
	
    }else{
	Rprintf("info[dpoco] = %d\n", *info);
	warning("[glmmboot:] Information non-positive definite. No variance!");
    }

    upper = 0;

/************** Bootstrapping starts *****************************/
 
    
    for (rep = 0; rep < *boot; rep++){
	if (*trace){
	    if ((rep / 10) * 10 == rep)
		Rprintf("********************* Replicate No. %d\n", rep);
	}
	if (*family <= 1){ /* Bernoulli */
	    indx = -1;
	    for (i = 0; i < ext->n_clust; i++){
		for (j = 0; j < clust[i].n; j++){
		    indx++;
		    clust[i].yw[j] = 
			rbinom((int)weights[indx], predicted[indx]);
		}
	    }
	}else{
	    indx = -1;
	    for (i = 0; i < ext->n_clust; i++){
		for (j = 0; j < clust[i].n; j++){ /* Poisson */
		    indx++;
		    clust[i].yw[j] = 
			rpois(weights[indx] * predicted[indx]);
		}
	    }
	}

	indx = 0;
	ant_fam_out = 0;
	for (i = 0; i < ext->n_clust; i++){
	    clust[i].ytot = 0.0;
	    for (j = 0; j < clust[i].n; j++){
		clust[i].ytot += clust[i].yw[j];
	    }
	    if (fabs(clust[i].ytot) < 0.001){
		clust[i].out = -1;
		clust[i].gamma = -1000.0; /* -inf */
	    }else if ((fabs(clust[i].wtot - clust[i].ytot) < 0.001) & 
		      (ext->family <= 1)){
		clust[i].out = 1;
		clust[i].gamma = 1000.0; /* +inf */
	    }else{
		ant_fam_out++;
		clust[i].out = 0;
	    }
	}

	if (!ant_fam_out){
	    /* Rprintf("Only trivial CLUSTERS!!!\n"); */
	    boot_log[rep] = 0.0;
	    if (0.0 >= *loglik) upper++;
	}else{


/* Restore beta as start values: */
/*	    for ( j = 0; j < *p; j++) b[j] = beta[j]; */

	    for ( j = 0; j < *p; j++) b[j] = 0.0;
	
	    if ( R_FINITE( tmp = bfun(*p, b, ext) ) ){
	        vmax = vmaxget();
		vmmin(*p, b, &Fmin,
		      bfun, bfun_gr, *maxit, *trace,
		      mask, abstol, reltol, nREPORT,
		      ext, &fncount, &grcount, &fail);
		vmaxset(vmax);
		*convergence = (fail == 0);
		boot_log[rep] = -Fmin;
		if (-Fmin >= *loglik) upper++;
	    }else{
		warning("Infinite start value to vmmin");
		Rprintf("\n");
	    }
	}	    
    }

    if (*boot) *boot_p = (double)upper / (double)*boot;
    else *boot_p = 1.0;

    PutRNGstate();

/*    vmaxset(vmax1); */

    Free(gr);
    Free(hessian);
    Free(hess_vec);
    Free(det);
    Free(work);

    for (i = 0; i < ext->n_clust; i++){
	Free(clust[i].yw);
	Free(clust[i].x);
	Free(clust[i].lin);
    }

    Free(clust);
    Free(ext);

    Free(mask);
    Free(b);
}
