#include <stdio.h>
#include <math.h>

#include "eha_fun.h"
#include <R_ext/Linpack.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

extern eha_P_fun *P;
extern eha_G_fun *G;
extern eha_Gprim_fun *Gprim;

/***********************************************************/
/*         Bernoulli distribution, logit link:             */

double eha_P_logit(double x, int y){ /* logit link */
    
    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;
    double res;

    res = plogis(x, location, scale, y, give_log);
    return ( res );
}

double eha_G_logit(double x, int y){
    
    /* Calculates G = P'/P */
    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;
    double res;

    if (y) {
	res = plogis(x, location, scale, 0, give_log);
    }else{
	res = -plogis(x, location, scale, 1, give_log);
    }
    return ( res );
}

double eha_Gprim_logit(double x, int y){ /* Note: Independent of y */
 
    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;

    return ( -dlogis(x, location, scale, give_log) );
}

/******************************************************************/
/*         Bernoulli distribution, cloglog link:                  */

double eha_P_cloglog(double x, int y){
    double res, s;
    
    s = exp(x);

    if (y)
	res = log1p(-exp(-s));
    else
	res = -s;

    return(res); /* Returns "log(P)" */
} 

double eha_G_cloglog(double x, int y){

    double s, res;
    
    s = exp(x);
    if (y)
	res = -s * (1.0 + 1.0 / expm1(-s));
    else
	res = -s;

    return(res);

}

double eha_Gprim_cloglog(double x, int y){

    double q, s;

    s = exp(x);
    q = exp(-s);
    
/*    return ( G_cloglog(x, yw, weight) - 
      yw * R_pow_di(s, 2) * q / R_pow_di(1.0 - q, 2) ); */
    return ( eha_G_cloglog(x, y) - 
	     y * R_pow_di(s, 2) * q / R_pow_di(expm1(-s), 2) ); 

}

/*****************************************************************/
/*       Poisson distribution, log link:                         */

double eha_P_poisson(double x, int y){

    int give_log = 0;
    return ( dpois(y, exp(x), give_log) );
}

double eha_G_poisson(double x, int y){

    return (y - exp(x));
}

double eha_Gprim_poisson(double x, int y){ /* Note: Independent of y */

    return (-exp(x));
}

/*****************************************************************/
/*           Normal distribution, identity link:                 */
/* Needs modification!!! ('int y' no good!)                      */

/*********************
double P_normal(double x, int y){

}

double G_normal(double x, int y){

}

double Gprim_normal(double x, int y){

}
********************/
static void eha_update(int level,
		       int p, 
		       double *beta,
		       double *loglik, 
		       double *score,
		       double *hessian,
		       int n,
		       double *x,
		       double *x_beta,
		       int *y,
		       int *haz,
		       int *riskset,
		       Exts *ext){
    
    double h;
    double *hb = NULL;
    double *hbb = NULL; /* (n_rs + p + 1) x (n_rs + p + 1) */
    double *pip = NULL; /* ext->n_points */
    double *xG = NULL;  /* ext->n_points * (p+1) */
    double sigma;
    double tmp, tmp2;
    int i, j, m, k;
    double factor;
    int count;
    int n_rs;
    int who;
    int bdim;

    n_rs = ext->n_rs;

    bdim = n_rs + p + 1;

    if (level < 0) return;

    factor = 1.0;

    sigma = beta[n_rs + p];

    pip = Calloc(ext->n_points, double);
    if (level > 0){
	xG = Calloc(bdim * ext->n_points, double);
	hb = Calloc(bdim, double);
	if (level > 1)
	    hbb = Calloc(bdim * bdim, double);
    }


/**********************************************************************/    
/* Calculate  h  for the loglik, and pip ("weights"): */
    h = 0.0;
    for (i = 0; i < ext->n_points; i++){
	tmp = 1.0;
	for (j = 0; j < n; j++){
	    tmp *= P(beta[haz[j]] + x_beta[j] + 
		     ext->zeros[i] * sigma, y[j]);
	}
	pip[i] = tmp * ext->weights[i];
	h += pip[i];
/*	Rprintf("pip[%d] = %f\n", i, pip[i]); */
    }

    /* Add into the loglik; We allow *loglik = -Inf!! */
/* But 'vmmin' doesn't allow it! Must do something about it. */

    *loglik += log(h);
    if (isinf(*loglik) == -1) 
    {
	warning("*loglik = -inf");
    }else if (isinf(*loglik) == 1){
	warning("*loglik = inf");
    }
    if (level == 0) {
	Free(pip);
	return;
    }

/********** First derivatives: ***********************************/

    if (h <= 0.0) { /* We don't allow h == 0.0 in the following! */
	Rprintf("h = 0.0; trying to fix...\n");
	h = 0.0;
	count = 0;
	while ( (h <= 0.0) && (count < 10) ){
	    count++;
	    factor = 2.0 * factor;
	    for (i = 0; i < ext->n_points; i++){
		tmp = 1.0;
		for (j = 0; j < n; j++){
		    tmp *= factor * P(beta[haz[j]] + x_beta[j] + 
				      ext->zeros[i] * sigma, y[j]);
		}
		pip[i] = tmp * ext->weights[i];
		h += pip[i];
	    }
	}
	if (h <= 0){
	    error("Unable to get likelihood function POSITIVE!!!!!!!!!\n");
	}
    }


/* Fill xG: */
    for (m = 0; m < n_rs; m++){ /* First for "haz" */

	for (i = 0; i < ext->n_points; i++){
	    tmp = 0.0;
	    for (j = 0; j < n; j++){
		if (haz[j] == m){
		    tmp += factor * /* x[m + who * p] * */  
			G(beta[m] + x_beta[j] + 
			  sigma * ext->zeros[i], y[j]);
		}
	    }
	    xG[i + m * ext->n_points] = tmp;
	}
    }
    for (m = n_rs; m < n_rs + p; m++){ /* Then for  */

	for (i = 0; i < ext->n_points; i++){
	    tmp = 0.0;
	    for (j = 0; j < n; j++){

		who = riskset[j];
		tmp += factor * x[(m - n_rs) + who * p] * 
		    G(beta[haz[j]] + x_beta[j] + 
		      sigma * ext->zeros[i], y[j]);
	    }
	    xG[i + m * ext->n_points] = tmp;
	}
    }
    for (i = 0; i < ext->n_points; i++){ /* Then for w = sigma */
	tmp = 0.0;
	for (j = 0; j < n; j++){
	    tmp += factor * ext->zeros[i] * 
		G(beta[haz[j]] + x_beta[j] + sigma * ext->zeros[i], y[j]);
	}
	xG[i + (n_rs + p) * ext->n_points] = tmp;
    }
    /* Done with xG */

/***********************************************************************/    
/* First derivatives, hb[]. */
    for (m = 0; m <= n_rs + p; m++){ /* hb[m]; note w = log(sigma) INCLUDED! */
	tmp = 0.0;
	for (i = 0; i < ext->n_points; i++){
	    tmp += pip[i] * xG[i + m * ext->n_points];
	    if (!isfinite(tmp)){
		Rprintf("pip[%d] = %f; xG[%d] = %f:: ", i, pip[i], 
			i + m * ext->n_points, xG[i + m * ext->n_points]); 
		Rprintf("tmp = %f\n", tmp);
		Rprintf("\nNumerical problem:\n"); 
		Rprintf("n.points (in 'control') is  %d\n",
			ext->n_points);
		Rprintf("Try a smaller value!\n");
		error("Execution interrupted");
	    }
	}
	hb[m] = tmp;
    }

    /* Add into first derivatives: */
    for (m = 0; m <= n_rs + p; m++){
	score[m] += hb[m] / h;
	if (!isfinite(score[m])){
	    Rprintf("Numerical problem:\n"); 
	    Rprintf("n.points (in 'control') is  %d\n",
		    ext->n_points);
	    Rprintf("Try a smaller value!\n");
	    error("Execution interrupted");
	    error("Non-finite score\n");
	}
    }
	
    if (level == 1){
	Free(xG);
	Free(pip);
	Free(hb);
	return;
    }
	
/***********************************************************************/
/* Done with first derivatives. On to the hessian: */

    /* First the n_rs x n_rs matrix of 'hazards': */
    for (m = 0; m < n_rs; m++){
	for (k = 0; k <= m; k++){
	    tmp = 0.0;
	    for (i = 0; i < ext->n_points; i++){
		tmp += pip[i] * xG[i + m * ext->n_points] *
		    xG[i + k * ext->n_points];
	    }
	    hbb[m + k * bdim] = tmp;
	    if (k == m){
		tmp = 0.0;
		for (i = 0; i < ext->n_points; i++){
		    tmp2 = 0.0;
		    for (j = 0; j < n; j++){
			if (m == haz[j]){
			    tmp2 += factor * 
				Gprim(beta[m] + x_beta[j] + 
				      sigma * ext->zeros[i], y[j]);
			}
		    }
		    tmp += tmp2 * pip[i];
		}		
		    /* And add in: */
		hbb[m + k * bdim] += tmp;
	    }
	}
    }

    /* Second, the n_rs x p matrix of 'haz + coeff' (no symmetry): */
    for (m = n_rs; m < n_rs + p; m++){
	for (k = 0; k < n_rs; k++){
	    tmp = 0.0;
	    for (i = 0; i < ext->n_points; i++){
		tmp += pip[i] * xG[i + m * ext->n_points] *
		    xG[i + k * ext->n_points];
	    }
	    hbb[m + k * bdim] = tmp;
	    tmp = 0.0;
	    for (i = 0; i < ext->n_points; i++){
		tmp2 = 0.0;
		for (j = 0; j < n; j++){
		    if (haz[j] == k){
			who = riskset[j];
			tmp2 += factor * x[(m - n_rs) + who * p] *
			    Gprim(beta[k] + x_beta[j] + 
				  sigma * ext->zeros[i], y[j]);
		    }
		}
		tmp += tmp2 * pip[i];
	    } /* And add in: */
	    hbb[m + k * bdim] += tmp;
	}
    }
    

    /* Third, the p x p matrix of 'coefficients' (symmetry): */
    for (m = n_rs; m < n_rs + p; m++){
	for (k = n_rs; k <= m; k++){
	    tmp = 0.0;
	    for (i = 0; i < ext->n_points; i++){
		tmp += pip[i] * xG[i + m * ext->n_points] *
		    xG[i + k * ext->n_points];
	    }
	    hbb[m + k * bdim] = tmp;
	    tmp = 0.0;
	    for (i = 0; i < ext->n_points; i++){
		tmp2 = 0.0;
		for (j = 0; j < n; j++){
		    who = riskset[j];
		    tmp2 += factor * x[(m - n_rs) + who * p] * 
			x[(k - n_rs) + who * p] *
			Gprim(beta[haz[j]] + x_beta[j] + 
			      sigma * ext->zeros[i], y[j]);
		}
		tmp += tmp2 * pip[i];
	    } /* And add in: */
	    hbb[m + k * bdim] += tmp;
	}
    }


    /* Then the "beta[m] with beta[bdim - 1] = sigma, m = 0,...,(p-1) */
    m = n_rs + p;

    /* 1. The hazards: */
    for (k = 0; k < n_rs; k++){ 
	tmp = 0.0;
	for (i = 0; i < ext->n_points; i++){
	    tmp += pip[i] * xG[i + m * ext->n_points] *
		xG[i + k * ext->n_points];
	}
	hbb[m + k * bdim] = tmp;
	tmp = 0.0; 
	for (i = 0; i < ext->n_points; i++){
	    tmp2 = 0.0;
	    for (j = 0; j < n; j++){ /* Here is the difference: */
		if (haz[j] == k){
		    tmp2 += factor * ext->zeros[i] * 
			Gprim(beta[k] + x_beta[j] + 
			      sigma * ext->zeros[i], y[j]);
		}
	    }
	    tmp += tmp2 * pip[i];
	} /* And add in: */
	hbb[m + k * bdim] += tmp;
    }

    /* 2. Coefficients: */
	
    for (k = n_rs; k < m; k++){
	tmp = 0.0;
	for (i = 0; i < ext->n_points; i++){
	    tmp += pip[i] * xG[i + m * ext->n_points] *
		xG[i + k * ext->n_points];
	}
	hbb[m + k * bdim] = tmp;
	tmp = 0.0; 
	for (i = 0; i < ext->n_points; i++){
	    tmp2 = 0.0;
	    for (j = 0; j < n; j++){ /* Here is the difference: */
		who = riskset[j];
		tmp2 += factor * x[(k - n_rs) + who * p] * 
		    ext->zeros[i] * 
		    Gprim(beta[haz[j]] + x_beta[j] + 
			  sigma * ext->zeros[i], y[j]);
	    }
	    tmp += tmp2 * pip[i];
	} /* And add in: */
	hbb[m + k * bdim] += tmp;
    }

    /* And finally beta[n_rs + p] with beta[n_rs + p] (AKA log(sigma) = w): */
    k = m; /* == n_rs + p */

    tmp = 0.0;
    for (i = 0; i < ext->n_points; i++){ /* Here is a difference: */
	tmp += pip[i] * xG[i + m * ext->n_points] *
	    (1.0 + xG[i + k * ext->n_points]);
    }
    hbb[m + k * bdim] = tmp;
    tmp = 0.0; 
    for (i = 0; i < ext->n_points; i++){
	tmp2 = 0.0;
	for (j = 0; j < n; j++){ /* Here is the difference: */
	    tmp2 += factor * 
		(ext->zeros[i]) * (ext->zeros[i]) *
		Gprim(beta[haz[j]] + x_beta[j] + 
		      sigma * ext->zeros[i], y[j]);
	}
	tmp += tmp2 * pip[i];
    } /* And add in: */

    hbb[m + k * bdim] += tmp;

   /* Now, add it into the hessian (lower triangle): */
    for (m = 0; m < bdim; m++){
	for (k = 0; k <= m; k++){
	    hessian[m + k * bdim] += hbb[m + k * bdim] / h - 
		(hb[m] / h) * (hb[k] / h);
	}
    }
    /* Fill in the upper triangle (symmetry): Necessary?? */
    for (m = 0; m < bdim; m++){
	for (k = (m + 1); k < bdim; k++){
	    hessian[m + k * bdim] = hessian[k + m * bdim];
	}
    }

/*********************************************************************/
/* We are done! Clean up 'the mess'! */

    Free(xG);
    Free(pip);
    Free(hb);
    Free(hbb);
}

static double eha_frail_mean(double sigma,
			     int n,
			     double *x_beta,
			     int *y,
			     Exts *ext){
    
    double h, h_mean;
    double tmp;
    int i, j;
    

/**********************************************************************/    
/* Calculate  h  for the loglik, and pip ("weights"): */
    h = 0.0;
    h_mean = 0.0;
    for (i = 0; i < ext->n_points; i++){
	tmp = 1.0;
	for (j = 0; j < n; j++){
	    tmp *= P(x_beta[j] + ext->zeros[i] * sigma, y[j]);
	}
	h += tmp * ext->weights[i];
	h_mean += tmp * ext->weights[i] * ext->zeros[i];
    }

    return ( h_mean / h);
}

void eha_frail_fun(int pp1, 
		   double *beta,
		   double *frail,
		   void *ex){

    int start;

    double tmp;
    int i, j;

    Exts *ext;
    /* double loglik; */

    int who;

    double sigma;

    ext = ex;

    /* loglik = 0.0; */

    for (i = 0; i < ext->n; i++){
	who = ext->riskset[i];
	tmp = ext->offset[who]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j + who * ext->p];
	}
	ext->x_beta[i] = tmp;
    }

    sigma = beta[ext->n_rs + ext->p];
/*
    F77_CALL(dgemm)(&trans, &trans, &(ex->n), 
    &p, &one, &alpha, x, &n, beta, &p,
		    &alpha, x_beta, &p);

    w = beta[ext->p];
    sigma = exp(w);
*/
    start = 0;
    for (i = 0; i < ext->n_fam; i++){

	frail[i] = eha_frail_mean(sigma,
				  ext->fam_size[i],
				  ext->x_beta + start,
				  ext->y + start,
				  ext);
	start += ext->fam_size[i];
    }
}

/*
void eha_mu_fun(int bdim, double *b, double *mu, void *ex){

    int i;
    Exts *ext;

    ext = ex;

    for (i = 0; i < ext->n; i++){
	mu[i] = 0.0;
    }
}
*/
double eha_fun(int bdim, 
	       double *beta, 
	       void *ex){

/* Dimensions:
   +++++++++++
   beta[p + 1] (0, ... p-1 = regression coefficients; p = log(sigma) = w
   x[p][n]
   y[n]
   fam_size[n_fam] { sum(fam_size) == n! }
   points[n_points][2] : first col: abscissas, second: weights.

   It is assumed that equal values in  'id'  comes together in groups.
*/

    int start;
    double *gr = NULL;

    double tmp;
    int i, j;

    Exts *ext;
    double loglik;

    double *hessian = NULL;

    int level = 0;
    int who;

    ext = ex;

    loglik = 0.0;

    for (i = 0; i < ext->n; i++){
	who = ext->riskset[i];
	tmp = ext->offset[who]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j + ext->n_rs] * ext->x[j + who * ext->p];
	}
	ext->x_beta[i] = tmp;
    }

/*
    F77_CALL(dgemm)(&trans, &trans, &(ex->n), 
    &p, &one, &alpha, x, &n, beta, &p,
		    &alpha, x_beta, &p);

    w = beta[ext->p];
    sigma = exp(w);
*/
    start = 0;
    
    for (i = 0; i < ext->n_fam; i++){

	eha_update(level,
	       ext->p,
	       beta,
	       &loglik,
	       gr,
	       hessian,
	       ext->fam_size[i],
	       ext->x,
	       ext->x_beta + start,
	       ext->y + start,
	       ext->haz + start,
	       ext->riskset + start,
	       ext);
	start += ext->fam_size[i];

    }
    return ( -loglik ); /* Note: minimizing!!! */
}


void eha_fun1(int bdim, 
	      double *beta,
	      double *gr,
	      void *ex){

    int i, j, k;
    int start;

    double loglik;
    double *hessian = NULL;

    Exts *ext;
    double tmp;
    int who;

    int level = 1;

    ext = ex;

    loglik = 0.0;

    for (k = 0; k < bdim; k++){
	gr[k] = 0.0;
    }

    for (i = 0; i < ext->n; i++){
	who = ext->riskset[i];
	tmp = ext->offset[who]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j + ext->n_rs] * ext->x[j + who * ext->p];
	}
	ext->x_beta[i] = tmp;
    }
    
    start = 0;
    for (i = 0; i < ext->n_fam; i++){

	eha_update(level,
	       ext->p,
	       beta,
	       &loglik,
	       gr,
	       hessian,
	       ext->fam_size[i],
	       ext->x,
	       ext->x_beta + start,
	       ext->y + start,
	       ext->haz + start,
	       ext->riskset + start,
	       ext);
	start += ext->fam_size[i];
    }

    for (i = 0; i < bdim; i++){

	gr[i] = -gr[i]; /* Minimization! */
	ext->gr[i] = gr[i];
    }
    /* Rprintf("\n"); */
}

void eha_fun2(int bdim, 
	      double *beta,
	      double *loglik,
	      double *gr,
	      double *hessian,
	      void *ex){

    int i, j, k;
    int start;
    double tmp;

    Exts *ext;

    int who;
    int level = 2;

    ext = ex;

    *loglik = 0.0;

    for (k = 0; k < bdim; k++){
	gr[k] = 0.0;
    }
    
    for (i = 0; i < bdim * bdim; i++){
	hessian[i] = 0.0;
    }

    for (i = 0; i < ext->n; i++){
	who = ext->riskset[i];
	tmp = ext->offset[who]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j + ext->n_rs] * ext->x[j + who * ext->p];
	}
	ext->x_beta[i] = tmp;
    }
    
    start = 0;
    for (i = 0; i < ext->n_fam; i++){

	eha_update(level,
	       ext->p,
	       beta,
	       loglik,
	       gr,
	       hessian,
	       ext->fam_size[i],
	       ext->x,
	       ext->x_beta + start,
	       ext->y + start,
	       ext->haz + start,
	       ext->riskset + start,
	       ext);
        start += ext->fam_size[i];
	}



   for (i = 0; i < bdim * bdim; i++){
	hessian[i] = -hessian[i];
    }

}


void eha_nr_opt(int bdim, double *beta, double *loglik, int *mask, 
		Exts *ext, double epsilon, int maxit, int trace){
    /* Start a Newton-Raphson thing: */
    
    double fstart, fprev;
    int iter;

    double *db = NULL;
    double *gr = NULL;
    double *hess = NULL;
    void *ex;

    int one = 1;
    int conver = 0;

    int info = 1;
    double L1;
    int i;

    int k, m;
    int true_bdim;

    int *ipiv;
    double *work;
    int lwork;

    double rcond;
    double *det;
    int job = 11;
    int pos_def;
    int maxiter = 10;

    det = Calloc(2, double);

    /* p = ext->p; */

    true_bdim = 0;
    for (i = 0; i < bdim; i++){
	true_bdim += mask[i];
    }
    if ( (true_bdim < (bdim - 1)) || (true_bdim > bdim) ) 
	error("Error in [nr_opt]: true dimension wrong.");

    db = Calloc(bdim, double);
    ipiv = Calloc(bdim, int);
    lwork = 11 * bdim;
    work = Calloc(lwork, double);

    gr = ext->gr;
    hess = ext->hessian;
    ex = ext;

    fprev = 0.0;
    fstart = 0.0;
    for (iter = 0; iter < maxiter; iter++){
	eha_fun2(bdim, beta, loglik, gr, hess, ex);
	if (iter ==  0){
	    fstart = *loglik;
	    fprev = fstart;
	}
  
	F77_CALL(dcopy)(&true_bdim, gr, &one, db, &one);

	pos_def = 0;
	while (!pos_def){
	    F77_CALL(dpoco)(hess, &bdim, &true_bdim, &rcond, work, &info);
	    
	    if (info){
		Rprintf("Hessian not positive definite.\n");
		Rprintf("info = %d\n", info);
		if (true_bdim == bdim){
		    eha_fun2(bdim, beta, loglik, gr, hess, ex);
		    Rprintf("We try fixing sigma at %f\n", fabs(beta[bdim-1]));
		    true_bdim--;
		}else{
		    Rprintf("sigma currently = %f", fabs(beta[bdim -1]));
		    error("Try another start value for sigma.\n");
		}
		F77_CALL(dpoco)(hess, &bdim, &true_bdim, 
				&rcond, work, &info);
		if (info) error("Try another start value for sigma.\n");
		pos_def = 1;
	    }else{
		pos_def = 1;
	    }
	}
	F77_CALL(dposl)(hess, &bdim, &true_bdim, db);
	
	L1 = 0.0; /*F77_CALL(dnrm2)(&bdim, db, &one);*/
	for (i = 0; i < true_bdim; i++){
	    L1 += fabs(db[i]);
	    beta[i] += db[i];
	}
	if (trace)
	    Rprintf("*** Iteration %d: L1 = %f, loglik = %f\n", iter, 
		   L1, *loglik);
	conver = (L1 < epsilon);
	if (!conver){ 
	    conver = 
		fabs(*loglik - fprev) < epsilon;/* * (fabs(fprev) + eps); */
	}


	if (conver){
	    if (iter > 0){
		if (trace){
		    Rprintf("Newton-Raphson CONVERGENCE in %d step(s)!!\n", 
			   iter);
		}
		break;
	    }
	}

	fprev = *loglik;

	if (*loglik < fprev){
	    Rprintf("Warning: Decreasing loglik!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
    }

    eha_fun2(bdim, beta, loglik, gr, hess, ex);
    F77_CALL(dpoco)(hess, &bdim, &true_bdim, &rcond, work, &info);
    if (info == 0){
	F77_CALL(dpodi)(hess, &bdim, &true_bdim, det, &job);

	for (m = 0; m < bdim; m++){
	    for (k = 0; k < m; k++){
		hess[m + k * bdim] = hess[k + m * bdim];
	    }
	}
    }else{
	Rprintf("No inversion in [nr_opt]\n");
    }

    Free(work);
    Free(ipiv);
    Free(db);
    Free(det);
}
