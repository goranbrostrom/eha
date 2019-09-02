#ifndef FRAIL_ML_H
#define FRAIL_ML_H

#ifndef MATHLIB_STANDALONE
#include <R.h>
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
    int n_rs;
    int *riskset;
    int nn;
    int *haz;
    double *x;           /* nn x p                            */
    double *offset;      /* nn                                */
    double *x_beta;      /* <- x %*% beta                    */
    double *gr;          /* p + 1                            */
    double *hessian;     /* (p+1) x (p+1)                    */
    int *y;              /* n-vector: binary response (0-1). */
    int *id;             /* n-vector.                        */
    int n_fam;           /* No. of families.                 */
    int *fam_size;       /* n_fam-vector.                    */
    int n_points;        /* No. of Gauss-Hermite points.     */
    double *weights;     /* n_points-vector.                 */
    double *zeros;       /* n_points-vector.                 */

}
Exts;

void eha_frail_ml(int *family,
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
		  double *sigma,
		  double *sigma_sd,
		  double *loglik,
		  double *variance,
		  double *frail,
		  /* double *mu, */
		  int *convergence,
		  int *fail);
    
#endif

