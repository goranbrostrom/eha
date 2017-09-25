#ifndef LOGISTIC_H
#define LOGISTIC_H

#ifndef MATHLIB_STANDALONE
#include <R.h>
#include <Rmath.h>
/* #include <R_ext/RS.h> */
#else
#include "testa.h"
#endif

typedef struct
{
    int family;          /* 0 = Bernoulli, logit link        */ 
                         /* 1 =            cloglog link      */
                         /* 2 = Poisson, log link            */

    int n;               /* = sum _0^(n_fam - 1) fam_size[i] */
    int p;               /* No. of covariates _including_    */ 
                         /* the constant (if any)            */
    int *cluster;        /* n                                */
    double **x;          /* p x n NOTE: check carefully **!! */
    double *offset;      /* n                                */

    double *x_beta;      /* <- x %*% beta                    */
    /* double *gr;          p + 1                            */
    /* double *hessian;     (p+1) x (p+1)                    */
    double *yw;          /* n-vector:  response * weights.   */
    double *weights;     /* n-vector:                        */
    double *cluster_weights;  /* n_fam-vector:                        */
    int n_fam;           /* No. of families.                 */
    int *fam_size;       /* n_fam-vector.                    */
    double *post_mode;   /* n_fam-vector.                    */
    double *post_mean;   /* n_fam-vector.                    */
    int n_points;        /* No. of Gauss-Hermite points.     */
    double *wc;          /* n_points-vector.                 */
    double *zeros;       /* n_points-vector.                 */

}
Exts;

typedef struct
{
    int n; /* fam_size */
    double sigma;
    double *x_beta;
    double *yw;
    double *weights;
    double cluster_weight;
    double **x;
    int p; /* == "Exts->p" */ 
    int m;     /* The actual partial derivative , m = 0, ..., (p - 1)  */
    int k;     /* The actual (second) partial derivative , */
               /* k = 0, ..., (p - 1)  */
}
Family;

void glmm_ml(int *family,
             int *p, 
             double *start_beta,
	     int *cluster,
	     double *weights,
	     double *cluster_weights,
             double *start_sigma,
	     int *fix_sigma,
             double *x, /* Now p x (\sum_1_{n_fam} fam_size[i]) */
             double *y,
             double *offset,
             int *fam_size,
             int *n_fam,
	     int *method,
             int *n_points, /* No. of pts in Gauss-Hermite quadrature */
             double *epsilon,
             int *maxit,
             int *trace,
	     int *boot,
	     int *prior,
	     double *predicted,
	     double *beta,
             double *sigma,
             double *loglik,
             double *variance,
             double *post_mode,
             double *post_mean,
             double *mu,
	     double *boot_p,
	     double* boot_log,
             int *convergence,
	     int *info);

#endif
