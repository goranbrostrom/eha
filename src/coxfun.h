#ifndef COXFUN_H
#define COXFUN_H
/*
static double gam1_fun(double gam, void *info);
*/

void get1_gam(RiskSet *risk);


void ml_rs(int what, RiskSet *risk,
	   double *b, double e_frac,
	   double *loglik, double *dloglik, double *d2loglik);

void breslow_rs(int what, RiskSet *risk, 
		double *b,double e_frac,
		double *loglik, double *dloglik, 
		double *d2loglik);

void efron_rs(int what, RiskSet *risk, 
	      double *b, double e_frac,
	     double *loglik, double *dloglik, 
	     double *d2loglik);
/*
static void cox_obs_rs(int what, RiskSet *risk,
		       double *b,
		       double *loglik, double *dloglik);
*/
void mppl_rs(int what, RiskSet *risk,
	     double *b, double e_frac,
	     double *loglik, double *dloglik, double *d2loglik);

void coxfun(int what, int totrs, RiskSet *risks, 
	    /*	    int ml, int method,  */ double e_frac,
	    double *b, 
	    double *loglik, double *dloglik, double *d2loglik);

typedef void RS_fun(int, RiskSet *,
		    double *, double,
		    double *, double *, double *);

#endif
