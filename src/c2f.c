/* Makes some C routines (mainly from R) available in Fortran */

#include <R.h>
#include "weibreg.h"

void F77_SUB(randstart)(void) { GetRNGstate();}
void F77_SUB(randend)  (void) { PutRNGstate();}
void F77_SUB(ranf)(double *x) { *x = unif_rand();}
void F77_SUB(rvsort)(double *x, int *index, int *n) {
	revsort(x, index, *n);
}
void F77_SUB(swfun)(int *order,
		    int *bdim, int *mb, double *beta, 
		    int *nn, double *z, double *time0, double *time, int *ind, 
		    double *offset, int *ns, int *nstra,
		    double *f, double *fp, double *fpp, int *iok) {
    sw_fun(order,
	   bdim, mb, beta, 
	   nn, z, time0, time, ind, 
	   offset, ns, nstra,
	   f, fp, fpp, iok);
}
