#include <stdio.h>
#include <R_ext/Linpack.h>

#include "glmmboot.h"
#include "bfun.h"
#include "fun.h"

extern P_fun *P;
extern G_fun *G;
extern H_fun *H;

void glmm_boot0(int *family,
		int *cluster,
		double *weights,
		double *y,
		double *offset,
		int *fam_size,
		int *n_fam,
		int *trace,
		int *boot,
		double *predicted,
		double *fitted,
		double *loglik,
		double *frail,
		double *boot_p,
		double *boot_log,
		int *convergence){
    
#ifndef MATHLIB_STANDALONE
    double abstol;
    /* double reltol; */
#endif
    Extb *ext;
    Cluster *clust;
    int i;
    int j;
    int cl;
    int rep;
    int indx;
    int ant_fam_out;

    double Fmin;
    double *b = NULL;
    int upper;

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
    /* reltol = abstol; */

    ext = Calloc(1, Extb);
    clust = Calloc(*n_fam, Cluster);
/************************ Fill in ext: *****************/
    ext->family = *family; /* == 0 for binomial(logit) */

    ext->n = 0;
    for (i = 0; i < *n_fam; i++){
	ext->n += fam_size[i];
    }

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
	for (i = 0; i < clust[cl].n; i++){
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

/* Don't need the following in a NULL model! */    
/* Note that this searches for a minimum: (!!) */
/*
    vmax = vmaxget();

    vmmin(*p, b, &Fmin,
	  bfun, bfun_gr, *maxit, *trace,
	  mask, abstol, reltol, nREPORT,
	  ext, &fncount, &grcount, &fail);
    *convergence = (fail == 0);
    vmaxset(vmax);
    bfun_gr(*p, b, gr, ext); 
*/
    Fmin = bfun(ext->p, b, ext);

    *loglik = -Fmin;

    if (*trace) Rprintf("loglik = %f\n", *loglik);

    for (i = 0; i < ext->n_clust; i++)
	frail[i] = clust[i].gamma;

/* Gone for now; predicted comes from calling R function
    if (ext->family <= 1){ 
	for (j = 0; j < ext->n; j++)
	    predicted[j] = ext->pred[j];
*/

    upper = 0;

/************** Bootstrapping starts *****************************/

    for (rep = 0; rep < *boot; rep++){
	if (*trace){
	    if ((rep / 10) * 10 == rep)
		Rprintf("********************* Replicate No. No. %d\n", i);
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
		boot_log[rep] = 0.0;           /* Fixed 090528 (i --> rep)*/
	    if (0.0 >= *loglik) upper++;
	}else{
	
	    Fmin = bfun(ext->p, b, ext);
	    boot_log[rep] = -Fmin;              /* Fixed 090528 (i --> rep)*/
	    if (-Fmin >= *loglik) upper++;
	}
    }

    if (*boot) *boot_p = (double)upper / (double)*boot;
    else *boot_p = 1.0;
    
    PutRNGstate();
    
/*    vmaxset(vmax1); */
    
    for (i = 0; i < ext->n_clust; i++){
	Free(clust[i].yw);
	Free(clust[i].x);
	Free(clust[i].lin);
    }

    Free(ext);
}
