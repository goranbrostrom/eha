#ifndef FUN_H
#define FUN_H

#ifdef MATHLIB_STANDALONE
#include "testa.h"
#endif

#include "glmmml.h"


void frail_fun(int pp1, 
		 double *beta,
	       void *ex);

void mu_fun(int bdim, 
	    double *b, 
	    double *mu, 
	    void *ex);

double fun(int pp1, 
	   double *beta, 
	   void *ex);

void fun1(int pp1, 
	  double *beta,
	  double *gr,
	  void *ex);

void fun2(int pp1, 
	  double *beta,
	  double *loglik,
	  double *gr,
	  double *hessian,
	  void *ex);

void nr_opt(int bdim, double *beta, double *loglik, int *mask, 
	    Exts *ext, double epsilon, int maxit, int *info, int trace);

typedef double P_fun(double, double, double);

typedef double G_fun(double, double, double);

typedef double H_fun(double, double, double);

typedef double I_fun(double, double, double);

typedef double K_fun(double, double, double);

typedef double logpr(double);

typedef double d_logpr(double);

typedef double d2_logpr(double);

typedef double d3_logpr(double);

typedef double d4_logpr(double);

double P_logit(double x, double yw, double weight); /* logit link */
    
double G_logit(double x, double yw, double weight);

double H_logit(double x, double yw, double weight);

double I_logit(double x, double yw, double weight);

double K_logit(double x, double yw, double weight);

double Hbis_logit(double x, double yw, double weight);

double P_cloglog(double x, double yw, double weight);

double G_cloglog(double x, double yw, double weight);

double H_cloglog(double x, double yw, double weight);

double I_cloglog(double x, double yw, double weight);

double K_cloglog(double x, double yw, double weight);

double P_poisson(double x, double yw, double weight);

double G_poisson(double x, double yw, double weight);

double H_poisson(double x, double yw, double weight);

double I_poisson(double x, double yw, double weight);

double K_poisson(double x, double yw, double weight);

double logprior_normal(double u);

double d_logprior_normal(double u);

double d2_logprior_normal(double u);

double d3_logprior_normal(double u);

double d4_logprior_normal(double u);

double logprior_logistic(double u);

double d_logprior_logistic(double u);

double d2_logprior_logistic(double u);

double d3_logprior_logistic(double u);

double d4_logprior_logistic(double u);

double logprior_cauchy(double u);

double d_logprior_cauchy(double u);

double d2_logprior_cauchy(double u);

double d3_logprior_cauchy(double u);

double d4_logprior_cauchy(double u);

#endif
