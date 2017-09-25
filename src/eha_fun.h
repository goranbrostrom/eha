#ifndef EHA_FUN_H
#define EHA_FUN_H

#ifdef MATHLIB_STANDALONE
#include "testa.h"
#endif

#include "frail_ml.h"



void eha_frail_fun(int pp1, 
		   double *beta,
		   double *frail,
		   void *ex);

void eha_mu_fun(int bdim, 
		double *b, 
		double *mu, 
		void *ex);

double eha_fun(int pp1, 
	       double *beta, 
	       void *ex);

void eha_fun1(int pp1, 
	      double *beta,
	      double *gr,
	      void *ex);

void eha_fun2(int pp1, 
	      double *beta,
	      double *loglik,
	      double *gr,
	      double *hessian,
	      void *ex);

void eha_nr_opt(int bdim, double *beta, double *loglik, int *mask, 
		Exts *ext, double epsilon, int maxit, int trace);

typedef double eha_P_fun(double, int);

typedef double eha_G_fun(double, int);

typedef double eha_Gprim_fun(double, int);

double eha_P_logit(double x, int y); /* logit link */
    
double eha_G_logit(double x, int y);

double eha_Gprim_logit(double x, int y);

double eha_P_cloglog(double x, int y);

double eha_G_cloglog(double x, int y);

double eha_Gprim_cloglog(double x, int y);

double eha_P_poisson(double x, int y);

double eha_G_poisson(double x, int y);

double eha_Gprim_poisson(double x, int y);

void F77_NAME(mlfun)(int *what, int *method,
		     int *totevent, int *totrs, int *ns, 
		     int *antrs, int *antevents, int *size,
		     int *totsize, int *eventset, int *riskset, 
		     int *nn, int *antcov, double *covar, double *offset,
		     double *beta, double *gamma,
		     double *loglik, double *h2, double *h22,
		     double *score);

void F77_NAME(coxfun)(int *what, int *method,
		      int *totevent, int *totrs, int *ns, 
		      int *antrs, int *antevents, int *size,
		      int *totsize, int *eventset, int *riskset, 
		      int *nn, int *antcov, double *covar, double *offset,
		      double *beta,
		      double *loglik, double *h2, double *h22,
		      double *score, double *sumdscore, double *sumd2score);


#endif
