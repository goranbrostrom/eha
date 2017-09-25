#ifndef GLMM_BOOT_H
#define GLMM_BOOT_H

#ifndef MATHLIB_STANDALONE
#include <R.h>
#include <Rmath.h>
/* #include <R_ext/RS.h> */
#else
#include "testa.h"
#endif

typedef struct
{
    int out; /* == 0: is counted in the analysis */
    int n;
    int p;
    double *weight;
    double wtot;
    double *offset;
    double **x;
    double *yw;
    double ytot;
    double *lin; /* = x_beta */
    double gamma; /* depends on beta */
}
Cluster;

typedef struct
{
    /* Fixed data: */
    int family;
    int n;               /* = sum _0^(n_fam - 1) fam_size[i] */
    int p;               /* No. of covariates _excluding_    */ 
                         /* a constant                       */
    int n_clust;         /* No. of clusters.                 */
    Cluster *clust;
}
Extb;


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
	       int *convergence);
    
#endif

