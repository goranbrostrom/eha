#ifndef LOGLIK_PHEXP_H
#define LOGLIK_PHEXP_H

void loglik_phexp(int *dis,
		  int *mb, 
		  double *b, 
		  double *alpha,
		  double *gamma,
		  int *nn, 
		  double *z, 
		  double *time0, 
		  double *time, 
		  int *ind, 
		  double *offset,
		  double *f);

void d_loglik_phexp(int *dis,
		    int *mb, 
		    double *b,
		    double *alpha,
		    double *gamma,
		    int *nn, 
		    double *z, 
		    double *time0, 
		    double *time, 
		    int *ind, 
		    double *offset, 
		    double *fp);

void d2_loglik_phexp(int *dis,
		     int *mb, 
		     double *b,
		     double *alpha,
		     double *gamma,
		     int *nn, 
		     double *z, 
		     double *time0, 
		     double *time, 
		     int *ind, 
		     double *offset, 
		     double *fpp);

#endif
