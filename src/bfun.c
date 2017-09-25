#include <stdio.h>
#include <math.h>

#include "glmmboot.h"
#include "fun.h"
#include "bfun.h"
#include <Rmath.h>
#include "GB_zeroin.h"
#include <R_ext/Linpack.h>

extern P_fun *P;
extern G_fun *G;
extern H_fun *H;

/* Binomial, logit link: */
static double get0_gam(Cluster *clust);

/* Binomial, cloglog link: */
static double get1_gam(Cluster *clust);

/* Poisson, log link: */
static double get2_gam(Cluster *clust);

/******************************************************************/

static double gam0_fun(double gam, void *info){
    Cluster *clust;
    double dg;
    int i;
    double x;
    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;

    clust = info;

    dg = clust->ytot;
/*    egam = exp(gam); */
    for (i = 0; i < clust->n; i++){
	x = gam + clust->lin[i];
	dg -= clust->weight[i] * plogis(x, location, scale, 1, give_log);
/*	egscore = egam * exp(ex->lin[i]);
	dg -= ex->weights[i] * egscore / ( 1.0 + egscore);
*/
    }
    return(dg);
}

static double get0_gam(Cluster *clust){

    /* Binomial, logit link */
    int i, itmax;

    double gam, eps;

    double gmin, gmax;
    double ax, bx;


    itmax = 35;
    eps = 0.00000001;

    gmin = clust->lin[0];
    gmax = clust->lin[0];
    for (i = 1; i < clust->n; i++){
	if (clust->lin[i] < gmin){ 
	    gmin = clust->lin[i];
	}else{
	    if (clust->lin[i] > gmax) gmax = clust->lin[i];
	}
    }
    gam = log(clust->ytot / (clust->wtot - clust->ytot)); /* start value */
    ax = gam - gmax;
    bx = gam - gmin;
    if (fabs (ax - bx) < eps) return((ax + bx) / 2.0);
    if (gam0_fun(ax, clust) * gam0_fun(bx, clust) > 0.0){
	Rprintf("f(%f) = %f, f(%f) = %f\n", 
		ax, gam0_fun(ax, clust), bx, gam0_fun(bx, clust));
	Rprintf("ytot = %f\n", clust->ytot); 
	Rprintf("wtot = %f\n", clust->wtot); 
	for (i = 0; i < clust->n; i++){
	    Rprintf("lin[%d] = %f\n", i, clust->lin[i]);
	    Rprintf("yw[%d] = %f\n", i, clust->yw[i]);
	    Rprintf("weights[%d] = %f\n", i, clust->weight[i]);
	}
	error("Wrong interval in [get0_gam]");
    }
    gam = GB_zeroin(ax, bx, &gam0_fun, clust, &eps, &itmax);

    return(gam);
}

static double gam1_fun(double gam, void *info){
    Cluster *clust;
    double dg;
    int i;
    clust = info;

    dg = 0.0;
    for (i = 0; i < clust->n; i++){
/*	s = exp(clust->lin[i] + gam);
	dg -= s * (clust->weight[i] + clust->yw[i] / expm1(-s));
*/	
	dg += G_cloglog(clust->lin[i] + gam, 
			clust->yw[i], 
			clust->weight[i]);
	
    }

    return(dg);
}

static double get1_gam(Cluster *clust){

    /* Binomial, cloglog link */
    int i, itmax;

    double gam, eps;

    double gmin, gmax;
    double ax, bx;

    itmax = 35;
    eps = 0.00000001;

    gmin = clust->lin[0];
    gmax = clust->lin[0];
    for (i = 1; i < clust->n; i++){
	if (clust->lin[i] < gmin){ 
	    gmin = clust->lin[i];
	}else{
	    if (clust->lin[i] > gmax) gmax = clust->lin[i];
	}
    }
    gam = log( -log(1.0 - clust->ytot / clust->wtot) ); /* start value */
    ax = gam - gmax;
    bx = gam - gmin;
    if (fabs (ax - bx) < eps) {
	return((ax + bx) / 2.0);
    }
    if (gam1_fun(ax, clust) * gam1_fun(bx, clust) > 0.0){
	Rprintf("f(%f) = %f, f(%f) = %f\n", 
		ax, gam1_fun(ax, clust), bx, gam1_fun(bx, clust));
	Rprintf("ytot = %f\n", clust->ytot); 
	Rprintf("wtot = %f\n", clust->wtot); 
	for (i = 0; i < clust->n; i++){
	    Rprintf("lin[%d] = %f\n", i, clust->lin[i]);
	    Rprintf("yw[%d] = %f\n", i, clust->yw[i]);
	    Rprintf("weights[%d] = %f\n", i, clust->weight[i]);
	}
	error("Wrong interval in [get0_gam]");
    }
    gam = GB_zeroin(ax, bx, &gam1_fun, clust, &eps, &itmax);

/*    if (*trace){
	Rprintf("gam = %f\n", gam);
    } */
    return(gam);
}

static double get2_gam(Cluster *clust){

/* For Poisson family */
    int j;
    double denom;

    denom = 0.0;
    for (j = 0; j < clust->n; j++){
	denom += clust->weight[j] * exp(clust->lin[j]);
    }

    return ( log( clust->ytot / denom ) );
}

double bfun(int p, double *b, void *ex){
    int i;
    int j;
    int cl;

    int indx;
    double loglik;

    Extb *ext;
    Cluster *clust;

    ext = ex;
    clust = ext->clust;
/* Get the "linear predictor": */

    for (cl = 0; cl < ext->n_clust; cl++){
	for (i = 0; i < clust[cl].n; i++){
	    clust[cl].lin[i] = clust[cl].offset[i];
	    for (j = 0; j < p; j++){
		clust[cl].lin[i] += b[j] * clust[cl].x[i][j];
	    }
	}
    }			   

/* Now get the gamma's: */
    
    indx = 0;

    if (ext->family <= 1){ /* binomial family */
	for (i = 0; i < ext->n_clust; i++){  /* NOT Excluding first family!! */
	    if (clust[i].out == 0){
		if (ext->family == 0){ /* logit link */
		    clust[i].gamma = get0_gam( clust + i );
		}else{ /* cloglog link */
		    clust[i].gamma = get1_gam( clust + i );
		} 
	    }
	    indx += clust[i].n;
	}
    }else{ /* Poisson; ext->family == 2 */
	for (i = 0; i < ext->n_clust; i++){ /* Excluding first family!! */
/* Why? ?2007-04-11. */ /* Changed on 2007-12-14 */
	    if (clust[i].out == 0){    
		clust[i].gamma = get2_gam( clust + i );
	    }
	    indx += clust[i].n;
	}
    }

/* Now get the log likelihood: */
    loglik = 0.0;
    /* Rprintf("beta[%d] = %f\n", 0, b[0]);  */
    for (i = 0; i < ext->n_clust; i++){
	if (clust[i].out == 0){
	    for (j = 0; j < clust[i].n; j++){
		loglik += P(clust[i].lin[j] + clust[i].gamma, 
			    clust[i].yw[j], clust[i].weight[j]);
	    }
	}
    }
	
    return(-loglik); /* Return: -loglikelihood */
    }

void bfun_gr(int n, double *b, double *gr, void *ex){
    int i;
    int j;
    int s;

    Extb *ext;
    Cluster *clust;

    ext = ex;
    
    clust = ext->clust;
    
/* Calculate the 'predicted values' at (beta, gamma): */
/* No, we don't! */
    for (s = 0; s < ext->p; s++){
	gr[s] = 0.0;
	for (i = 0; i < ext->n_clust; i++){
	    if (clust[i].out ==  0){
		for (j = 0; j < clust[i].n; j++){
		    gr[s] += clust[i].x[j][s] * 
			G(clust[i].gamma + clust[i].lin[j], 
			  clust[i].yw[j], clust[i].weight[j]);  
		}
	    }
	}
    }


    /* Return minus the gradient! */
    for (s = 0; s < n; s++) gr[s] = -gr[s];
}

void bfun_hess(int p, double *b, double *hessian, Extb *ext){

    int m, s, i, j, indx;
    Cluster *clust;
    
    double *h, *h_fam;
    double **hess;
    double gam, t1, t2;
    
    clust = ext->clust;

    h = Calloc(ext->n, double);
    h_fam = Calloc(ext->n_clust, double);
    hess = Calloc(p, double *);
    for (m = 0; m < p; m++){
	hess[m] = hessian + m * p;
    }

    for (i = 0; i < ext->n; i++) h[i] = 0.0;

    indx = -1;
    for (i = 0; i < ext->n_clust; i++){
	h_fam[i] = 0.0;

	if (clust[i].out == 0){
	    gam = clust[i].gamma;
	    for (j = 0; j < clust[i].n; j++){
		indx++;
		h[indx] = H(clust[i].lin[j] + gam, clust[i].yw[j],
				clust[i].weight[j]);
		h_fam[i] += h[indx];
	    }
	}else{
	    indx += clust[i].n;
	}
    }
    
    for (m = 0; m < p; m++){
	for (s = 0; s <= m; s++){
	    hess[m][s] = 0.0;
	}
    }
    
    for (m = 0; m < p; m++){
	for (s = 0; s <= m; s++){
	    indx = -1;
	    for (i = 0; i < ext->n_clust; i++){
		for (j = 0; j < clust[i].n; j++){
		    indx++;
		    hess[m][s] += clust[i].x[j][m] * 
		    clust[i].x[j][s] * h[indx];
		}
	    }

	    indx = -1;
	    for (i = 0; i < ext->n_clust; i++){
		if (clust[i].out == 0){
		    t1 = 0.0;
		    t2 = 0.0;
		    for (j = 0; j < clust[i].n; j++){
			indx++;
			t1 += clust[i].x[j][m] * h[indx];
			t2 += clust[i].x[j][s] * h[indx];
		    }
		    hess[m][s] -= t1 * t2 / h_fam[i];
		}else{
		    indx += clust[i].n;
		}
	    }
	}
    }
    
    for (m = 0; m < p; m++){
	hess[m][m] = -hess[m][m];     /* Added 0.65-4 */
	for (s = m + 1; s < p; s++){             /* (an old error!) */
	    hess[s][m] = -hess[s][m]; /* Added 0.65-4 */
	    hess[m][s] = hess[s][m];
	}
    }
    Free(hess);
    Free(h_fam);
    Free(h);
}
   
