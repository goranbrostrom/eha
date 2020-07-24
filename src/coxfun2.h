#ifndef COXFUN2_H
#define COXFUN2_H
/*
static double gam1_fun(double gam, void *info);
*/

/* void get1_gam(RiskSet *risk); */

/*
void ml_rs2(int what, RiskSet *risk,
	   double *b, double e_frac,
	   double *loglik, double *dloglik, double *d2loglik);
*/
void breslow_rs2(int *what, /* RiskSet *risk, */
		 int *antevents,
		 int *size,
		 double *weights,
		 double *x,
		 double *lin,
		 int *p,
		 double *b,double *e_frac,
		 /* Return: */
		 double *loglik, double *dloglik, 
		 double *d2loglik);

void efron_rs2(int *what, /* RiskSet *risk, */
	       int *antevents,
	       int *size,
	       double *weights,
	       double *x,
	       double *lin,
	       int *p,
	       double *b, double *e_frac,
	       /* Return: */
	       double *loglik, double *dloglik, 
	       double *d2loglik);
/*
static void cox_obs_rs(int what, RiskSet *risk,
		       double *b,
		       double *loglik, double *dloglik);

void mppl_rs2(int what, RiskSet *risk,
	     double *b, double e_frac,
	     double *loglik, double *dloglik, double *d2loglik);
*/

#endif
