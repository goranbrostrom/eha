#ifndef BFUN_H
#define BFUN_H

#include "glmmboot.h"

double bfun(int nvars, double *b, void *ext);

void bfun_gr(int n, double *b, double *gr, void *ext);

/*
double get_gam(int n, 
	       double *lin, 
	       int ytot);
*/
void bnr_opt(int bdim, double *beta, double *loglik, int *mask, 
	     Extb *ext, double epsilon, int maxit, int trace);
/*
double get2_gam(int n, 
	       double *lin,
		int ytot);

double get3_gam(int n,
	       double *lin, 
		int ytot);
*/
void bfun_hess(int p, double *b, double *hessian, Extb *ext);


#endif
