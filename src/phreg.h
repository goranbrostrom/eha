#ifndef PHREG_H
#define PHREG_H

#include "phfun.h"

int dist;         /* Which distribution?           */
ph0S_fun *S0;      /* Survival fcn (standardized!)  */
ph0_fun *f0;      /* Density fcn                   */
ph0_fun *h0;      /* hazard fcn                    */
ph0_fun *f0_t;     /* First order derivativa of f0  */
ph0_fun *h0_t;     /* First order derivative of h0  */
ph0_fun *h0_tt;    /* Second order derivative of h0 */

typedef struct{
    int *ns;
    int *nstra;
    double *pfix;
    int *mb;
    int *nn;
    double *z;
    double *time0;
    double *time;
    int *ind;
    double *offset;
    double *f;
    double *fp;
    double *fpp;
    int *iok;
} Exts;

/*
static double ph_fun(int n, double *beta, void *vex);

static void gph_fun(int n, double *beta, double *dloglik, void *vex);

static void g2ph_fun(int n, double *beta, double *d2loglik, void *vex);

static void ph_nr(int iter, double eps, int printlevel, 
		  int bdim, double *beta, 
		  double *loglik, double *dloglik, double *variance,
		  int *conver, int *fail, Exts *ex);
*/
void phsup(int *iter, double *eps, int *printlevel,
	   int *ns, int *nstra, int *nn, int *ncov, int *bdim,
	   double *time0, double *time, int * ind,
	   double *covar, double *offset, int *dist, /* 'dist' new */
	   double *init, double *beta, double *lambda, double *lambda_sd,
	   double *shape, double *shape_sd,
	   double *loglik, double *dloglik, double *variance, double *sctest,
	   int *conver, int *fail);

#endif
