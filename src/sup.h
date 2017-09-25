#ifndef SUP_H
#define SUP_H

/* Global variables */

/* Fixed: */
extern int p;
extern int nn;
extern double *x; /* Covariates: (nn x p) */
/* extern double *offset; */
/* Variable: */ 
extern double *lin;    /* nn */
/* extern double *score;      nn */
extern double *sumdscore;   /* p */
extern double *sumd2score;   /* p x p */

/* End global variables */

typedef struct
{
    int out; /* == 0: is counted in the analysis */
    int stratum;
    double time; /* Failure time for those who fail in this rs */
    int antevents; /* How many events in this risk set? */
    int *eventset; /* Who died in this risk set? */
    int size;      /* No. of members of the risk set */   
    double *weights;
    double *offset; /* Time-varying offset */
    double rs_weight; /* riskset mean of weights; */
    int *riskset;  /* Members of the risk set */
    double gamma;  /* depends on beta and data*/
    double hazard;
    double tot_score; /* sum(score(j)), j in RiskSet; */
}
RiskSet;


void sup(int *meth, 
	 int *iter, 
	 double *eps, 
	 int *prl, 
	 int *totevent, 
	 int *totrs, 
	 int *ns, 
	 int *antrs, 
	 int *antevents, 
	 int *size,
	 double *weights,
	 int *totsize, 
	 int *eventset, 
	 int *riskset, 
	 int *nn_in, 
	 int *p_in, 
	 double *covar, 
	 double *offset_in,
	 double *startbeta,
	 int *boot,
	 double *efrac,
	 double *beta,
	 double *sd_beta,
	 double *loglik, 
	 double *variance, 
	 double *sctest,
	 double *hazard,
	 int *conver, 
	 int *f_conver, 
	 int *fail);

#endif
