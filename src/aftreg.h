#ifndef AFTREG_H
#define AFTREG_H

#include "phfun.h"

int dist;         /* Which distribution?           */
ph0S_fun *S0;      /* Survival fcn (standardized!)  */
ph0_fun *f0;      /* Density fcn                   */
ph0_fun *h0;      /* hazard fcn                    */
ph0_fun *f0_t;     /* First order derivativa of f0  */
ph0_fun *h0_t;     /* First order derivative of h0  */
ph0_fun *h0_tt;    /* Second order derivative of h0 */

typedef struct{
    int *id; /* Primary sorting key */
    int *strata; /* Numbered 0, ..., (*ns - 1) */
    int *ns;
    double *pfix;
    int *mb;
    int *nn;
    double *z;
    double *time0; /* Secondary sorting key */
    double *time;
    int *ind;
    double *offset;
} Exts;


void aftsup(int *printlevel,
	    int *ns, int *nn, int *ncov, int *bdim,
	    int *id, int *strata, double *time0, double *time, int *ind,
	    double *covar, double *offset, int *dis, double *beta, 
	    double *loglik, int *fail);

#endif
