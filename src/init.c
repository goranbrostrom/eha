#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void aftexpsup(int *, int *, int *, int *, int *,
		      int *, int *, double *, double *, int *,
		      double *, double *, double *, int *, 
		      double *, double *, int *);
extern void aftregGomp(int *, int *, int *, int *, int *,
		       int *, int *, double *, double *, int *,
		       double *, double *, int *, double *,
		       double *, int *);
extern void aftsup(int *, int *, int *, int *, int *,
		   int *, int *, double *, double *, int *,
		   double *, double *, int *, double *, 
		   double *, int *);
extern void d_loglik_ph(int *, int *, double *, double *,
			double *, int *, double *, double *, 
			double *, int *, double *, double *);
extern void d2_loglik_ph(int *, int *, double *, double *,
			 double *, int *, double *, double *, 
			 double *, int *, double *, double *);
extern void expsup(int *iter, double *eps, int *printlevel,
		   int *nn, int *ncov, int *bdim,
		   double *time0, double *time, int * ind,
		   double *covar, double *offset, double *shape,
		   double *init, double *beta, double *lambda, double *lambda_sd,
		   double *loglik, double *dloglik, double *variance, double *sctest,
		   int *conver, int *fail);
extern void frail_ml(int *family, int *method, int *p, int *nn, int *n_rs,
		     int *riskset, double *start_beta, double *start_sigma,
		     double *x, int *y, int * haz, double *offset, int *fam_size,
		     int *n_fam, int *n_points, double *epsilon, int *maxit,
		     int *trace, double *beta, double *sigma, double *loglik,
		     double *variance, double *frail, /* double *mu, */
		     int *convergence, int *fail);

extern void loglik_ph(int *, int *, double *, double *, double *,
		      int *, double *, double *, double *, 
		      int *, double *, double *);
extern void phexpsup(int *iter, double *eps, int *printlevel,
		     int *ns, int *nstra, int *nn, int *ncov, int *bdim,
		     double *time0, double *time, int * ind,
		     double *covar, double *offset, double *shape, int *dis, 
		     double *init, double *beta, double *lambda, double *lambda_sd,
		     double *loglik, double *dloglik, double *variance, double *sctest,
		     int *conver, int *fail);
extern void phsup(int *iter, double *eps, int *printlevel,
		  int *ns, int *nstra, int *nn, int *ncov, int *bdim,
		  double *time0, double *time, int * ind,
		  double *covar, double *offset, int *dist, /* 'dist' new */
		  double *init, double *beta, double *lambda, double *lambda_sd,
		  double *shape, double *shape_sd,
		  double *loglik, double *dloglik, double *variance, double *sctest,
		  int *conver, int *fail);
extern void risk_get(int *max_s, int *nn, int *ns,
		     double *enter, double *exit, int*event, 
		     int* nstra, int *l_nstra,
		     int *new_totrs,
		     int *antrs, int *n_events, int *size,
		     double *risktimes,
		     int *eventset, int *riskset);
extern void sizes(int *ns, int *nn, double *enter, double *exit, int *event,
		  int *antrs, int *nstra, double *risktimes, 
		  int *n_events,int *size, int *totrs);
extern void sup(int *meth, 
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
extern void weibsup(int *iter, double *eps, int *printlevel,
		    int *ns, int *nstra, int *nn, int *ncov, int *bdim,
		    double *time0, double *time, int * ind,
		    double *covar, double *offset,
		    double *init, double *beta, double *lambda, double *lambda_sd,
		    double *shape, double *shape_sd,
		    double *loglik, double *dloglik, double *variance, double *sctest,
		    int *conver, int *fail);

/* .Fortran calls */
/* Where is 'bootcox'???
extern void F77_NAME(bootcox)(void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *);
*/
extern void F77_NAME(chek)(int *, int *, int *, double *, double *, int *,
			   double *, int *);
extern void F77_NAME(cleanup)(double *, double *, double *, int *, int *,
			      int *, int *, int *, double *, int *,
			      double *, double *, double *, int *, int *);
extern void F77_NAME(geomsup)(void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *, void *);
extern void F77_NAME(ghq)(int *, double *, double *, int *);
extern void F77_NAME(hazards)(int *, int *, int *, int *, int *,
			      int *, int *, int *, double *, double *);
extern void F77_NAME(martres)(int *, int *, int *, int *, int *,
			      int *, int *, int *, double *, double *,
			      double *);
/*
extern void F77_NAME(phfuncnull)(void *, void *, void *, void *, void *,
				 void *, void *, void *, void *, void *,
				 void *, void *, void *);
*/
extern void F77_NAME(split)(double *, int *, int *, double *, int *,
			    int *, int *, int *, double *, int *);
extern void F77_NAME(wfunc)(int *, int *, double *, int *, int *, double *,
			    int *, double *, double *, double *, int *, double *,
			    double *, double *, double *, int *);
extern void F77_NAME(wfuncnull)(int *, int *, double *, int *, double *,
				int *, double *, double *, int *, double *,
				double *, double *, int *);

static const R_CMethodDef CEntries[] = {
    {"aftexpsup",    (DL_FUNC) &aftexpsup,    17},
    {"aftregGomp",   (DL_FUNC) &aftregGomp,   16},
    {"aftsup",       (DL_FUNC) &aftsup,       16},
    {"d_loglik_ph",  (DL_FUNC) &d_loglik_ph,  12},
    {"d2_loglik_ph", (DL_FUNC) &d2_loglik_ph, 12},
    {"expsup",       (DL_FUNC) &expsup,       22},
    {"frail_ml",     (DL_FUNC) &frail_ml,     27},
    {"loglik_ph",    (DL_FUNC) &loglik_ph,    12},
    {"phexpsup",     (DL_FUNC) &phexpsup,     25},
    {"phsup",        (DL_FUNC) &phsup,        26},
    {"risk_get",     (DL_FUNC) &risk_get,     15},
    {"sizes",        (DL_FUNC) &sizes,        11},
    {"sup",          (DL_FUNC) &sup,          30},
    {"weibsup",      (DL_FUNC) &weibsup,      25},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
/*    {"bootcox",    (DL_FUNC) &F77_NAME(bootcox),    32}, */
    {"chek",       (DL_FUNC) &F77_NAME(chek),        8},
    {"cleanup",    (DL_FUNC) &F77_NAME(cleanup),    15},
    {"geomsup",    (DL_FUNC) &F77_NAME(geomsup),    27},
    {"ghq",        (DL_FUNC) &F77_NAME(ghq),         4},
    {"hazards",    (DL_FUNC) &F77_NAME(hazards),    10},
    {"martres",    (DL_FUNC) &F77_NAME(martres),    11},
/*    {"phfuncnull", (DL_FUNC) &F77_NAME(phfuncnull), 13}, */
    {"split",      (DL_FUNC) &F77_NAME(split),      10},
    {"wfunc",      (DL_FUNC) &F77_NAME(wfunc),      16},
    {"wfuncnull",  (DL_FUNC) &F77_NAME(wfuncnull),  13},
    {NULL, NULL, 0}
};

void R_init_eha(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
