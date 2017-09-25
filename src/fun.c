#include <stdio.h>
#include <math.h>

#include "glmmml.h"
#include "fun.h"
#include <Rmath.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include "GB_zeroin.h"

extern int laplace; /* 0 = Gauss-Hermite, 1 = Laplace */
extern P_fun *P;
extern G_fun *G;
extern H_fun *H;
extern I_fun *I;
extern K_fun *K;
extern logpr *logprior;
extern d_logpr *d_logprior;
extern d2_logpr *d2_logprior;
extern d3_logpr *d3_logprior;
extern d4_logpr *d4_logprior;


static double g(double u, void *ex);

static double g_u(double u, void *ex);

static double g_uu(double u, void *ex);

static double g_s(double u, void *ex);

static double g_m(double u, int m, void *ex);

static double g_us(double u, void *ex);

static double g_um(double u, int m, void *ex);

static double g_uus(double u, void *ex);

static double g_uum(double u, int m, void *ex);

/***********************************************************/
/*         Bernoulli distribution, logit link:             */

double P_logit(double x, double yw, double weight){ /* logit link */

/* Note: Changed to return log(p) !! (2006-08-21) */

    double res, p;

    res = x * yw - weight * log1p(exp(x));

    if ((yw > 0.001) & ((weight - yw) > 0.001)){ 
	p = yw / weight;
	res = res - (yw * log(p) + (weight - yw) * log(1.0 - p));
    }

    return ( res );
}

double G_logit(double x, double yw, double weight){
    
    /* Calculates G = P'/P */

    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;

/* yw - weight * (exp(x)/(1 + exp(x))) */

    return ( yw - weight * plogis(x, location, scale, 1, give_log) );
}

double H_logit(double x, double yw, double weight){ 
/* Note: Independent of yw */
 
    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;

    return ( -weight * dlogis(x, location, scale, give_log) );
}

double I_logit(double x, double yw, double weight){ 

    double location = 0.0;
    double scale = 1.0;
    int give_log = 0;

/* Note that expm1(x) = exp(x) - 1 = -(1 - exp(x)) !! */
/* -n * exp(x) / (1 + exp(x))^2 * (1 - exp(x)) / (1 + exp(x)) */
    return ( weight * dlogis(x, location, scale, give_log) * 
	     expm1(x) / (1.0 + exp(x)) );
}

double K_logit(double x, double yw, double weight){ 
    
    double s;
    
    s = exp(x);
/*
    n * exp(x) * (4 * exp(x) - exp(2 * x) - 1) / (1 + exp(x))^4
*/ 
    return ( weight * s * (4.0 * s - R_pow_di(s, 2) - 1.0) / 
	     R_pow_di(1 + s, 4) );
	     
}

/******************************************************************/
/*         Bernoulli distribution, cloglog link:                  */

    double P_cloglog(double x, double yw, double weight){

/* Note: Changed to return log(p) !! (2006-08-21) */

    double s, q, p, res;

    s = exp(x);
    q = exp(-s);
    
    res = yw * log1p(-q) - (weight - yw) * s;

    if ((yw > 0.001) & ((weight - yw) > 0.001)){
	/* Rprintf("yw = %f, weight = %f\n", yw, weight); */
	p = yw / weight;
	res = res - (yw * log(p) + (weight - yw) * log(1.0 - p));
    }

    return( res );
}

double G_cloglog(double x, double yw, double weight){

    double s, res;

    s = exp(x);
/*    q = exp(-s); */

/*    return ( s * (yw / (1.0 - q) - weight) ); */
/*    return( s * (yw / (1.0 - exp(-s)) - weight) ); */
    res = -s * (yw / expm1(-s) + weight);

    return(res);

 
}

double H_cloglog(double x, double yw, double weight){
 
    double q, s;

    s = exp(x);
    q = exp(-s);

/*    return ( G_cloglog(x, yw, weight) - 
      yw * R_pow_di(s, 2) * q / R_pow_di(1.0 - q, 2) ); */
    return ( G_cloglog(x, yw, weight) - 
	     yw * R_pow_di(s, 2) * q / R_pow_di(expm1(-s), 2) ); 
}

double I_cloglog(double x, double yw, double weight){
    
    double q, s;

    s = exp(x);
    q = exp(-s);

    return ( H_cloglog(x, yw, weight) -
	     yw * R_pow_di(s, 2) * q * 
	     (2.0 - s * (1 + q)/ R_pow_di((1.0 - q), 3)) ); 
	     /* s * q * ( 2.0 * (1.0 - q) - s * (1 + q)) /
		R_pow_di(1 - q, 3) );*/
	     }
    
double K_cloglog(double x, double yw, double weight){

    double q, s, s2, z;
    double mq3;


    s = exp(x);
    s2 = R_pow_di(s, 2);
    q = exp(-s);
    mq3 = R_pow_di(1-q, 3);

    z = (yw * (2 * s) * q - yw * s2 * q) * 
	(2 - s * (1 + q) / mq3) -
      yw * s2 * q * (((1 + q) - s * q)/ mq3 -
		     s * (1 + q) * 3 * q / R_pow_di((1 - q), 4));

    return ( I_cloglog(x, yw, weight) - s * z );
}
	
/*****************************************************************/
/*       Poisson distribution, log link:                         */

double P_poisson(double x, double yw, double weight){

    double res;

    res = weight * (x * yw - exp(x));
    if (yw > 0.001) res -= weight * yw * (log(yw) - 1.0);
 
    return ( res );
}

double G_poisson(double x, double yw, double weight){

    return ( weight * (yw - exp(x)) );
}

double H_poisson(double x, double yw, double weight){ 
/* Note: Independent of yw */

    return ( -weight * exp(x) );
}

double I_poisson(double x, double yw, double weight){ 
/* Note: Independent of yw */

    return ( -weight * exp(x) );
}

double K_poisson(double x, double yw, double weight){ 
/* Note: Independent of yw */

    return ( -weight * exp(x) );
}

/* Prior distribution and it's derivatives:            */
/* no scale; transformed away! (May be reconsidered... */
/* Normal, Logistic, and Cauchy for the moment...      */

double logprior_normal(double u){

/* \log phi(u) = -0.5 log(2 * pi) - u^2 / 2*/

    return ( -M_LN_SQRT_2PI - R_pow_di(u, 2) / 2.0 );
}

double d_logprior_normal(double u){

    return ( -u );

}

double d2_logprior_normal(double u){

    return ( -1.0 );

}

double d3_logprior_normal(double u){

    return ( 0.0 );
}

double d4_logprior_normal(double u){

    return ( 0.0 );
}

double logprior_logistic(double u){

    return ( u - 2 * log(1+exp(u)) );
/*
 return( u - 2.0 * plogis(u, location, scale, lower_tail, give_log) );
*/
}

double d_logprior_logistic(double u){

    int give_log = 0; /* NOT log scale */
    double location = 0.0;
    double scale = 1.0;
    int lower_tail = 1; /* Yes, lower tail */

    /* 1 - 2 * exp(u) / (1 + exp(u)) */
    return ( 1.0 - 2.0 * plogis(u, location, scale, lower_tail, give_log) );

}

double d2_logprior_logistic(double u){

    int give_log = 0;
    double location = 0.0;
    double scale = 1.0;

/* -2.0 * exp(u) / (1.0 + exp(u))^2 */

    return ( -2.0 * dlogis(u, location, scale, give_log) );

}

double d3_logprior_logistic(double u){

    double s;

    s = exp(u);

/* -2.0 * exp(u) * (1.0 - exp(u)) / (1.0 + exp(u))^3 */

    return ( -2.0 * s * (1.0 - s) / R_pow_di(1 + s, 3) );

}

double d4_logprior_logistic(double u){

    double s, s2;

    s = exp(u);
    s2 = R_pow_di(s, 2);

    return ( 2.0 * s * (2.0 * s - 1.0) / R_pow_di(1 + s, 3) +
	6.0 * s2 * (1.0 - s) / R_pow_di(1 + s, 4) );

}

double logprior_cauchy(double u){

  /* M_PI = pi */

    return( -log(M_PI) -log(1.0 + R_pow_di(u, 2)) );
}

double d_logprior_cauchy(double u){

  /* M_PI = pi */

    return( -2.0 * u / (1.0 + R_pow_di(u, 2)) );
}

double d2_logprior_cauchy(double u){

    double u2;

    u2 = R_pow_di(u, 2);

    return( 2.0 * (u2 - 1) / R_pow_di(u2 + 1, 2) );
}

double d3_logprior_cauchy(double u){

    double u2;

    u2 = R_pow_di(u, 2);

    return( 4.0 * u * (3.0 - u2) / R_pow_di(u2 + 1, 3) );
}

double d4_logprior_cauchy(double u){

    double u2;

    u2 = R_pow_di(u, 2);

    return( (12.0 - 72.0 * u2 + 12.0 * R_pow_di(u, 4)) / 
	     R_pow_di(u2 + 1, 4) );
}


/*****************************************************************/
/*           Normal distribution, identity link:                 */
/* Needs modification!!! ('int y' no good!)                      */

/*********************
double P_normal(double x, double y, double weight){

}

double G_normal(double x, int y){

}

double H_normal(double x, int y){

}
********************/

/****************************************************************/
/* The functions for Laplace comes here. Can also be useful for */
/* Gauss-Hermite (fully adapted).                               */

static double g(double u, void *ex){

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;

    sigu = fam->sigma * u; 
    res = 0.0;
    for (j = 0; j < fam->n; j++){
	res += P(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
    }

    /* Rprintf("[g:] res = %f\n", res); */
    return ( logprior(u) + res );
}

static double g_u(double u, void *ex){

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;

    sigu = fam->sigma * u; 
    
    res = 0.0;
    for (j = 0; j < fam->n; j++){
	res += G(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
    }
    
    return ( d_logprior(u) + res * fam->sigma );
}

static double g_uu(double u, void *ex){

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;

    sigu = fam->sigma * u; 

    res = 0.0;
    for (j = 0; j < fam->n; j++){
	res += H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
    }
    

    return ( d2_logprior(u) +  res * R_pow_di(fam->sigma, 2) );
}

static double g_s(double u, void *ex){

/* Component 'sigma' */

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;

    sigu = fam->sigma * u; 

    res = 0.0;
    for (j = 0; j < fam->n; j++){
	res += G(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
    }
    

    return ( res * u );
}

static double g_m(double u, int m, void *ex){

/* Component 'm' in beta */

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;

    if (fam->p == m){
	res = g_s(u, ex);
    }else{
	
	sigu = fam->sigma * u; 
	
	res = 0.0;
	for (j = 0; j < fam->n; j++){
	    res += fam->x[m][j] * 
		G(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
	}
    }
    return ( res );
}

static double g_ss(double u, void *ex){

/* Component 'sigma' */

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;

    sigu = fam->sigma * u; 

    res = 0.0;
    for (j = 0; j < fam->n; j++){
	res += H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
    }
    

    return ( res * R_pow_di(u, 2) );
}

static double g_sm(double u, int m, void *ex){

/* Component 'm' in beta */

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;
   
    if (fam->p == m){
	res = g_ss(u, ex);
    }else{
	sigu = fam->sigma * u; 
	
	res = 0.0;
	for (j = 0; j < fam->n; j++){
	    res += fam->x[m][j] * 
		H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
	}
	res = res * u;
    }
    
    return ( res );
}

static double g_mk(double u, int m, int k, void *ex){

/* Component 'm' in beta */

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;
    
    if ((m == k) && (m == fam->p)){
	res = g_ss(u, ex);
    }else if (m == fam->p){
	res = g_sm(u, k, ex);
    }else if (k == fam->p){
	res = g_sm(u, m, ex);
    }else{
	sigu = fam->sigma * u; 
	
	res = 0.0;
	for (j = 0; j < fam->n; j++){
	    res += fam->x[m][j] * fam->x[k][j] * 
		H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
	}
    }
    return ( res );
}


static double g_us(double u, void *ex){

/* Component 'm' in beta */

    int j;
    double sum1, sum2, sigu;
    Family *fam;

    fam = ex;

    sigu = fam->sigma * u; 

    sum1 = 0.0;
    sum2 = 0.0;
    for (j = 0; j < fam->n; j++){
	sum1 +=  H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
	sum2 +=  G(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
    }
    
    return ( sum1 * sigu + sum2);
}

static double g_um(double u, int m, void *ex){

/* Component 'm' in beta */

    int j;
    double res, sigu;
    Family *fam;

    fam = ex;

    if (fam->p == m){
	res = g_us(u, ex);
    }else{
	sigu = fam->sigma * u; 
	
	res = 0.0;
	for (j = 0; j < fam->n; j++){
	    res += fam->x[m][j] * 
		H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
	}
	res = fam->sigma * res;
    }
    
    return ( res );
}

static double g_uus(double u, void *ex){

/* Component 'm' in beta */

    int j;
    double sum1, sum2, sigu;
    Family *fam;

    fam = ex;

    sigu = fam->sigma * u; 

    sum1 = 0.0;
    sum2 = 0.0;
    for (j = 0; j < fam->n; j++){
	sum1 +=  H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
	sum2 += I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
    }
    
    return ( 2 * fam->sigma * sum1 + sigu * fam->sigma * sum2 );
}

static double g_uum(double u, int m, void *ex){

/* Component 'm' in beta */

    int j;
    double res, sigu;

    Family *fam;

    fam = ex;

    if (fam->p == m){
	res = g_uus(u, ex);
    }else{
	sigu = fam->sigma * u; 
	
	res = 0.0;
	for (j = 0; j < fam->n; j++){
	    res += fam->x[m][j] * 
		I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]); 
	}
	res = R_pow_di(fam->sigma, 2) * res;
    }
    
    return ( res );
}

static double g_uuu(double u, void *ex){
    int j;
    double sum, sigu;

    Family *fam;

    fam = ex;

    sigu = fam->sigma * u;

    sum = 0.0;
    for (j = 0; j < fam->n; j++){
	sum += I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
    }

    return ( d3_logprior(u) + sum * R_pow_di(fam->sigma, 3) );
} 

static double g_uss(double u, void *ex){

    int j;
    double sum1, sum2, sigu;

    Family *fam;

    fam = ex;

    sigu = fam->sigma * u;

    sum1 = 0.0; sum2 = 0.0;
    for (j = 0; j < fam->n; j++){
	sum1 += I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	sum2 += H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
    }
    
    return (sigu * fam->sigma * sum1 + 2.0 * u * sum2);
}

static double g_usm(double u, int m, void *ex){

    int j;
    double res, sum1, sum2, sigu;

    Family *fam;

    fam = ex;

    if (fam->p == m){
	res = g_uss(u, ex);
    }else{
	sigu = fam->sigma * u;
	
	sum1 = 0.0; sum2 = 0.0;
	for (j = 0; j < fam->n; j++){
	    sum1 += fam->x[m][j] *
		I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	    sum2 += fam->x[m][j] *
		H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	}
	res = sigu * sum1 + sum2;
    }

    return ( res );
}

static double g_umk(double u, int m, int k, void *ex){

    int j;
    double res, sum, sigu;

    Family *fam;

    fam = ex;

    if ((m == k) && (m == fam->p)){
	res = g_uss(u, ex);
    }else if (m == fam->p){
	res = g_usm(u, k, ex);
    }else if (k == fam->p){
	res = g_usm(u, m, ex);
    }else{
	sigu = fam->sigma * u;

	sum = 0.0;
	for (j = 0; j < fam->n; j++){
	    sum += fam->x[m][j] * fam->x[k][j] *
		I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	}
	res = fam->sigma * sum;
    }

    return ( res );
}

static double g_uuss(double u, void *ex){

    int j;
    double sum1, sum2, sum3, sigu;

    Family *fam;

    fam = ex;

    sigu = fam->sigma * u;

    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0;
    for (j = 0; j < fam->n; j++){
	sum1 += K(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	sum2 += I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	sum3 += H(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
    }

    return(R_pow_di(sigu, 2) * sum1 +
	   4.0 * sigu * sum2 +
	   2.0 * sum3);
}

static double g_uusm(double u, int m, void *ex){

    int j;
    double res, sum1, sum2, sigu;

    Family *fam;

    fam = ex;

    if (m == fam->p){
	res = g_uuss(u, ex);
    }else{
	sigu = fam->sigma * u;

	sum1 = 0.0; sum2 = 0.0;
	for (j = 0; j < fam->n; j++){
	    sum1 += fam->x[m][j] * 
		K(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	    sum2 += fam->x[m][j] *
		I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	}
	res = sigu * fam->sigma * sum1 + 2.0 * fam->sigma * sum2;
    }

    return( res );
}

static double g_uumk(double u, int m, int k, void *ex){

    int j;
    double res, sum, sigu;

    Family *fam;

    fam = ex;

    if ((m == k) && (m == fam->p)){
	res = g_uuss(u, ex);
    }else if (m == fam->p){
	res = g_uusm(u, k, ex);
    }else if (k == fam->p){
	res = g_uusm(u, m, ex);
    }else{

	sigu = fam->sigma * u;

	sum = 0.0;
	for (j = 0; j < fam->n; j++){
	    sum += fam->x[m][j] * fam->x[k][j] *
		K(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	}
	
	res = R_pow_di(fam->sigma, 2) * sum;
    }

    return( res );
}

static double g_uuuu(double u, void *ex){

    int j;
    double sum, sigu;

    Family *fam;

    fam = ex;

    sigu = fam->sigma * u;

    sum = 0.0;
    for (j = 0; j < fam->n; j++){
	sum += K(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
    }

    return(d4_logprior(u) + R_pow_di(fam->sigma, 4) * sum);
}

static double g_uuus(double u, void *ex){

    int j;
    double sum1, sum2, s2, sigu;

    Family *fam;

    fam = ex;

    s2 = R_pow_di(fam->sigma, 2);
    sigu = fam->sigma * u;

    sum1 = 0.0; sum2 = 0.0;
    for (j = 0; j < fam->n; j++){
	sum1 += K(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	sum2 += I(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
    }

    return(sigu * s2 * sum1 + 3.0 * s2 * sum2);
}

static double g_uuum(double u, int m, void *ex){

    int j;
    double res, sum, s3, sigu;

    Family *fam;

    fam = ex;

    if (fam->p == m){
	res = g_uuus(u, ex);
    }else{
	s3 = R_pow_di(fam->sigma, 3);
	sigu = fam->sigma * u;
	
	sum = 0.0;
	for (j = 0; j < fam->n; j++){
	    sum += fam->x[m][j] *
		K(fam->x_beta[j] + sigu, fam->yw[j], fam->weights[j]);
	}
	
	res = s3 * sum;
    }

    return(res);
}

/*************** End "Laplace" (and Gauss-hermite) fun's ***************/


static void update(int level,
		   int p, 
		   double *beta,
		   double *loglik, 
		   double *score,
		   double *hessian,
		   double *post_mode,
		   Family *fam,
		   int n_points,
		   double *wc,
		   double *zeros){
    
    double h = 0.0;
    double *hb = NULL;
    double *hbb = NULL; /* (p+1) x (p+1) */
    /* double sigma; */
    double tmp, tmp1, tmp2;
    int i, m, k;

    /* For integrate */
    int limit = 100;
    int lenw;
    int *iwork;
    double *work;

    /* End For integrate */
    
    void *ex;
    double ax, bx;
    /* int n; */
    /* double **x; */
    /* double *x_beta; */
    /* double *yw; */

    double u_hat, sigma_hat, sigma2_hat;
    double reltol = 0.000001;
    int maxit = 100;

    /* temporary gig: */
    double *u;
    double *wght;

    double u_s, u_ss;
    double *u_m, *u_sm;
    double *u_mkvec;
    double **u_mk;

    double sh_s, sh_ss; /* sh: short for 'sigma_hat' ! */
    double *sh_m, *sh_sm;
    double *sh_mkvec;
    double **sh_mk;

    double gu, guu, guuu, guuuu; /* Counterparts to g_u, g_uu, etc */

    if (level < 0) return;

    u = Calloc(n_points, double);
    wght = Calloc(n_points, double);

    for (i = 0; i < n_points; i++){
/* Should be moved to 'ghq.f' after testing! DONE now!!*/
	wght[i] = wc[i]; /* * exp(R_pow_di(zeros[i], 2)); */
    }

    ex = fam;

    /* n = fam->n; */
    /* x = fam->x; */
    /* x_beta = fam->x_beta; */
    /* yw = fam->yw; */

    lenw = 4 * limit;
    iwork = Calloc(limit, int);
    work = Calloc(lenw, double);

    /* sigma = beta[p]; */

    /* First, find the u value that maximizes the integrand */
    /* (or, equivalently, the log of it:                    */ 

    tmp = 0.0;
    for (i = 0; i < fam->n; i++){
	tmp += fam->x_beta[i];
    }

    tmp /= (double)fam->n;
    if (fabs(tmp) < 1.e-16 || fabs(fam->sigma) < 1.e-16) u_hat = 0;
    else u_hat = -tmp / fam->sigma; /* Start value for u_hat */

    /* tmp = fun_fam(one, &u_hat, ex); */
    tmp = g(u_hat, ex);

    ax = -1.0;
    bx = -ax;
    while (g_u(ax, ex) * g_u(bx, ex) >= 0){
	ax = 2.0 * ax;
	bx = 2.0 * bx;
    }
    u_hat = GB_zeroin(ax, bx, g_u, ex, &reltol, &maxit);
    
    gu = g_u(u_hat, ex);
    guu = g_uu(u_hat, ex);
    if (guu >= 0.0) 
	error("[update:] Second derivative non-negative at max!");
    sigma_hat = sqrt(-1.0 / guu);
    
    sigma2_hat = R_pow_di(sigma_hat, 2);
    
    for (i = 0; i < n_points; i++){
	u[i] = M_SQRT2 * sigma_hat * zeros[i] + u_hat;
    }
    *post_mode = u_hat;
    /* *post_mode = u_hat * sigma_hat;  Changed back 2009-02-11, Version 0.81-4. */

/**********************************************************************/    
/* Calculate  h  for the loglik: */
    
/* We choose GHQ for the loglik....

    Rdqagi(int_fun, fam, &bound, &inf, &epsabs, &epsrel, &h,
	   &abserr, &neval, &ier, &limit, &lenw, &last,
	   iwork, work);
*/

/* GHQ: */
    
    if (!laplace){
	
	h = 0.0;
	for (i = 0; i < n_points; i++){
	    h +=  wght[i] * exp(g(u[i], ex));
	}
	/* M_LN2 = log(2) */

	*loglik += 
	    fam->cluster_weight * 
	    (0.5 * M_LN2 + log(sigma_hat) + log(h));

    }else{ /* Laplace */
	 *loglik += 
	     fam->cluster_weight * 
	     (g(u_hat, ex) + M_LN_SQRT_2PI + log(sigma_hat));
    }

    if (level == 0) {
	Free(u);
	Free(wght);
	Free(iwork);
	Free(work);
	return;
    }


/***********************************************************************/    
/* First order derivatives, hb[]. */

/* level >= 1 */

    u_m = Calloc(p + 1, double); /* must be changed for GHQ (n_points > 1)!! */

    for (i = 0; i < p; i++) 
	u_m[i] = sigma2_hat * g_um(u_hat, i, ex);
    u_s = sigma2_hat * g_us(u_hat, ex);
    u_m[p] = u_s;

    sh_m = Calloc(p + 1, double); /* must be changed for GHQ (n_points > 1)!! */

    guuu = g_uuu(u_hat, ex);
    for (i = 0; i < p; i++) 
	sh_m[i] = 0.5 * R_pow_di(sigma_hat, 3) * 
	    (u_m[i] * guuu + g_uum(u_hat, i, ex)); 
    sh_s = 0.5 * R_pow_di(sigma_hat, 3) * 
	(u_s * guuu + g_uus(u_hat, ex));
    sh_m[p] = sh_s;
/*******************************************************************/

    hb = Calloc((p + 1), double);


    for (m = 0; m <= p; m++){ /* hb[m]; note w = log(sigma) INCLUDED! */
	if (m < p){

/*	    if (n_points >= 2){ */
	    if (!laplace){
		tmp = 0.0;
		for (i = 0; i < n_points; i++){
		    tmp2 = (zeros[i] * sh_m[m] * M_SQRT2 + u_m[m]) *
			g_u(u[i], ex) + g_m(u[i], m, ex);
		    tmp += tmp2 * wght[i] * exp(g(u[i], ex));
		}
		hb[m] = tmp;
	    }else{ /* Laplace */
		hb[m] = u_m[m] * gu + g_m(u_hat, m, ex);
	    }  

	}else{ /* p == m */

	    /* GHQ: */
	    /*   if (n_points >= 2){ */
	    if (!laplace){
		tmp = 0.0;
		for (i = 0; i < n_points; i++){
		    tmp2 = (zeros[i] * sh_s * M_SQRT2 + u_s) *
			g_u(u[i], ex) + g_s(u[i], ex);
		    tmp += tmp2 * wght[i] * exp(g(u[i], ex));
		}
		hb[m] = tmp;
		
	    }else{ /* Laplace */
		hb[m] = u_s * gu + g_s(u_hat, ex);
	    }
	}
    }

    
    /* Add into first derivatives: KOLLA HERE!??!!*/

    if (laplace){
	for (m = 0; m < p; m++){
	    score[m] += 
		fam->cluster_weight * 
		(sh_m[m] / sigma_hat + hb[m]);
	}
	score[p] += fam->cluster_weight * (sh_s / sigma_hat + hb[p]);
    }else{
	for (m = 0; m < p; m++){
	    score[m] += 
		fam->cluster_weight * 
		(sh_m[m] / sigma_hat + hb[m] / h);
	}
	score[p] += fam->cluster_weight * (sh_s / sigma_hat + hb[p] / h);
    }

    if (level == 1){
	Free(u);
	Free(wght);
	Free(hb);
	Free(iwork);
	Free(work);
	Free(u_m);
	Free(sh_m);

	return;
    }
	
/***********************************************************************/
/* Done with first derivatives. On to the hessian: */

/* level >= 2 */

    u_sm = Calloc(p + 1, double);
    u_mkvec = Calloc((p + 1) * (p + 1), double);
    u_mk = Calloc(p + 1, double *);
    for (m = 0; m <= p; m++){
	u_mk[m] = u_mkvec + m * (p + 1);
    }

    u_ss = sigma2_hat * (R_pow_di(u_s, 2) * guuu + 
			 2.0 * u_s * g_uus(u_hat, ex) +
			 g_uss(u_hat, ex));
    u_mk[p][p] = u_ss;
    /* Rprintf("u_ss = %f\n", u_ss); */
    for (m = 0; m < p; m++){
	u_sm[m] = sigma2_hat * (u_s * u_m[m] * guuu + 
				u_s * g_uum(u_hat, m, ex) +
				u_m[m] * g_uus(u_hat, ex) +
				g_usm(u_hat, m, ex));
	u_mk[m][p] = u_sm[m];
	u_mk[p][m] = u_sm[m];
	for (k = 0; k < p; k++){
	    u_mk[m][k] = sigma2_hat * (u_m[k] * u_m[m] * guuu + 
				       u_m[k] * g_uum(u_hat, m, ex) +
				       u_m[m] * g_uum(u_hat, k, ex) +
				       g_umk(u_hat, m, k, ex));
	    /*  Rprintf("u_mk[%d][%d] = %f\n", m, k, u_mk[m, k]); */
	}
    }
/*    for (m = 0; m < p; m++)Rprintf("u_sm[%d] = %f\n", m, u_sm[m]); */

    sh_sm = Calloc(p + 1, double);
    sh_mkvec = Calloc((p + 1) * (p + 1), double);
    sh_mk = Calloc(p + 1, double *);
    for (m = 0; m <= p; m++){
	sh_mk[m] = sh_mkvec + m * (p + 1);
    }

    guuuu = g_uuuu(u_hat, ex);
    sh_ss = 0.75 * R_pow(sigma_hat, 5) * 
	(u_s * guuu + g_uus(u_hat, ex)) *
	(u_s * guuu + g_uus(u_hat, ex))
	+ 0.5 * R_pow_di(sigma_hat, 3) *
	(u_ss * guuu + u_s *u_s * guuuu + 
	 u_s * g_uuus(u_hat, ex) + u_s * g_uuus(u_hat, ex) + 
	 g_uuss(u_hat, ex));
    sh_mk[p][p] = sh_ss;

    for (m = 0; m < p; m++){
	sh_sm[m] = 0.75 * R_pow(sigma_hat, 5) * 
	    (u_m[m] * guuu + g_uum(u_hat, m, ex)) *
	    (u_s * guuu + g_uus(u_hat, ex))
	    + 0.5 * R_pow_di(sigma_hat, 3) *
	    (u_sm[m] * guuu + u_s *u_m[m] * guuuu + 
	     u_s * g_uuum(u_hat, m, ex) + u_m[m] * g_uuus(u_hat, ex) + 
	     g_uusm(u_hat, m, ex));
	sh_mk[m][p] = sh_sm[m];
	sh_mk[p][m] = sh_sm[m];
	for (k = 0; k < p; k++){
	    sh_mk[m][k] = 0.75 * R_pow(sigma_hat, 5) * 
	    (u_m[m] * guuu + g_uum(u_hat, m, ex)) *
	    (u_m[k] * guuu + g_uum(u_hat, k, ex))
	    + 0.5 * R_pow_di(sigma_hat, 3) *
	    (u_mk[m][k] * guuu + 
	     u_m[m] *u_m[k] * guuuu + 
	     u_m[k] * g_uuum(u_hat, m, ex) + u_m[m] * g_uuum(u_hat, k, ex) + 
	     g_uumk(u_hat, m, k, ex));
	}
    }

/********************************************************************/

    hbb = Calloc((p + 1) * (p + 1), double);

    /* First the pxp matrix of 'coefficients': */
    /* Now: Everything at once!! (2008-05-24) */
    if (!laplace){
	for (m = 0; m <= p; m++){
	    for (k = 0; k <= m; k++){
		
		hbb[m + k * (p + 1)] = 0.0;
		fam->m = m;
		fam->k = k;
		for(i = 0; i < n_points; i++){

		    tmp = u_mk[m][k] + M_SQRT2 * sh_mk[m][k] * zeros[i];
		    tmp1 = u_m[m] + M_SQRT2 * sh_m[m] * zeros[i];
		    tmp2 = u_m[k] + M_SQRT2 * sh_m[k] * zeros[i];

		    hbb[m + k * (p + 1)] += wght[i] * (
			(tmp * g_u(u[i], ex) + 
			tmp1 * (g_um(u[i], k, ex) + tmp2 * g_uu(u[i], ex)) +
			tmp2 * g_um(u[i], m, ex) + g_mk(u[i], m, k, ex) 
			) * exp(g(u[i], ex)) +
			(tmp1 * g_u(u[i], ex) + g_m(u[i], m, ex)) *
			(tmp2 * g_u(u[i], ex) + g_m(u[i], k, ex)) *
			exp(g(u[i], ex))
			);
						     
		}
		
	    }
	}
    }else{ /* Laplace 'complete' */
	for (m = 0; m <= p; m++){
	    for (k = 0; k <= m; k++){
		fam->m = m;
		fam->k = k;
		hbb[m + k * (p + 1)] =
		    (sh_mk[m][k] * sigma_hat - sh_m[m] * sh_m[k]) /
		    sigma2_hat +
		    u_mk[m][k] * gu +
		    u_m[m] * (u_m[k] * guu + g_um(u_hat, k, ex)) +
		    u_m[k] * g_um(u_hat, m, ex) + 
		    g_mk(u_hat, m, k, ex);		    
		    /*    l_mk(u_hat, sigma2_hat, m, k, ex) +
		    l_u(u_hat, sigma2_hat, ex) * 
		    u_mk[m][k] +
		    l_uu(u_hat, sigma2_hat, ex) * u_m[m] * u_m[k] +
		    l_um(u_hat, sigma2_hat, m, ex) * u_m[k] +
		    l_um(u_hat, sigma2_hat, k, ex) * u_m[m];
		    */
	    }
	}
    }

   /* Now, add it into the hessian (lower triangle): */
    /* if (n_points >= 2){ */
    if (!laplace){
	for (m = 0; m <= p; m++){
	    for (k = 0; k <= m; k++){
		hessian[m + k * (p+1)] += 
		    fam->cluster_weight * (
		    (sh_mk[m][k] * sigma_hat - sh_m[m] * sh_m[k]) /
		    sigma2_hat +

		    hbb[m + k * (p+1)] / h  -   
		    (hb[m] / h) * (hb[k] / h)); /* 0 at the solution? No!!! */
	    }
	}
    }else{ /* Laplace */
	for (m = 0; m <= p; m++){
	    for (k = 0; k <= m; k++){
		hessian[m + k * (p+1)] += 
		    fam->cluster_weight * (hbb[m + k * (p+1)]); /* - */     
		/*  hb[m] * hb[k]);  */
/* 0 at the solution? No!!! */
	    }
	}
    }
    /* Fill in the upper triangle (symmetry): Necessary?? */
    for (m = 0; m <= p; m++){
	for (k = (m + 1); k <= p; k++){
	    hessian[m + k * (p+1)] = hessian[k + m * (p+1)];
	}
    }

/*********************************************************************/
/* We are done! Clean up 'the mess'! */
    Free(u);
    Free(wght);
    Free(hbb);
    Free(hb);
    Free(iwork);
    Free(work);

    Free(u_m);
    Free(u_sm);
    Free(u_mkvec);
    Free(u_mk);

    Free(sh_m);
    Free(sh_sm);
    Free(sh_mkvec);
    Free(sh_mk);

}

static double frail_mean(int level,
			 int p, 
			 double *beta,
			 double *loglik, 
			 double *score,
			 double *hessian,
			 Family *fam,
			 int n_points,
			 double *wc,
			 double *zeros){
    
    double h, h_mean;
    double sigma;
    double tmp;
    int i, j;
    
    sigma = beta[p];

/**********************************************************************/    
/* Calculate  h  for the loglik, and pip ("wc"): */
    h = 0.0;
    h_mean = 0.0;
    for (i = 0; i < n_points; i++){
	tmp = 1.0;
	for (j = 0; j < fam->n; j++){
	    tmp *= exp(P(fam->x_beta[j] + zeros[i] * sigma, 
		     fam->yw[j], fam->weights[j]));
	}
	h += tmp * wc[i];
	h_mean += tmp * wc[i] * zeros[i];
/*
	Rprintf("zeros[%d] = %f\n", i, zeros[i]);
	Rprintf("h_mean = %f, h = %f\n", h_mean, h);
*/
    }

    return ( h_mean / h);
}

void frail_fun(int pp1, 
		 double *beta,
		 void *ex){

    int start;
    double *gr = NULL;

    double tmp;
    int i, j;

    Exts *ext;
    Family *fam;

    double loglik;

    double *hessian = NULL;

    int level = 0;

    ext = ex;

    fam = Calloc(1, Family);
    fam->x = Calloc(ext->p, double *);
    fam->p = ext->p;
    fam->m = 0;
    fam->k = 0;
    fam->sigma = beta[ext->p];

    loglik = 0.0;

    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j][i];
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
	fam->n = ext->fam_size[i];
	fam->cluster_weight = ext->cluster_weights[i];
	fam->x_beta = ext->x_beta + start;
	fam->yw = ext->yw + start;
	fam->weights = ext->weights + start;
	for (j = 0; j < ext->p; j++)	
	    fam->x[j] = ext->x[j] + start;  /* check this! */

	ext->post_mean[i] = frail_mean(level,
				       ext->p,
				       beta,
				       &loglik,
				       gr,
				       hessian,
				       fam,
				       ext->n_points,
				       ext->wc,
				       ext->zeros);
	
	start += ext->fam_size[i];
    }
    Free(fam->x);
    Free(fam);
}

void mu_fun(int bdim, double *b, double *mu, void *ex){

    int i;
    Exts *ext;

    ext = ex;

    for (i = 0; i < ext->n; i++){
	mu[i] = 0.0; /* More to come later */
    }
}

double fun(int pp1, 
	   double *beta, 
	   void *ex){

/* Dimensions:
   +++++++++++
   beta[p + 1] (0, ... p-1 = regression coefficients; p = sigma
   x[p][n]
   y[n]
   fam_size[n_fam] { sum(fam_size) == n! }
   points[n_points][2] : first col: abscissas, second: wc.

   It is assumed that equal values in  'id'  comes together in groups.
*/

    int start;
    double *gr = NULL;
    double post_mode;

    char trans = 'N';
    double alpha = 1.0;
    int one = 1;

    int i, j;

    Exts *ext;
    Family *fam;

    double loglik;

    double *hessian = NULL;

    int level = 0;

    ext = ex;

    fam = Calloc(1, Family);
    fam->x = Calloc(ext->p, double *);
    fam->p = ext->p;
    fam->m = 0;
    fam->k = 0;
    fam->sigma = beta[ext->p];

    loglik = 0.0;
/*
    for (i = 0; i < ext->p + 1; i++) 
	Rprintf("[fun] beta[%d] = %f\n", i, beta[i]);
*/
/*
    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j][i];
	}
	ext->x_beta[i] = tmp;
    }
*/

    F77_CALL(dcopy)(&(ext->n), ext->offset, &one, ext->x_beta, &one);
    F77_CALL(dgemv)(&trans, &(ext->n), &(ext->p), 
    &alpha, ext->x[0], &(ext->n), beta, &one,
		    &alpha, ext->x_beta, &one);


    start = 0;
    for (i = 0; i < ext->n_fam; i++){
	fam->n = ext->fam_size[i];
	fam->cluster_weight = ext->cluster_weights[i];
	fam->x_beta = ext->x_beta + start;
	fam->yw = ext->yw + start;
	fam->weights = ext->weights + start;
	for (j = 0; j < ext->p; j++)	
	    fam->x[j] = ext->x[j] + start;  /* check this! */
	update(level,
	       ext->p,
	       beta,
	       &loglik,
	       gr,
	       hessian,
	       &post_mode,
	       fam,
	       ext->n_points,
	       ext->wc,
	       ext->zeros);
	ext->post_mode[i] = post_mode;
	start += ext->fam_size[i];
    }
    Free(fam->x);
    Free(fam);
    /* Rprintf("[fun]; -loglik = %f\n", -loglik); */ 
    return ( -loglik ); /* Note: minimizing!!! */
}


void fun1(int pp1, 
	  double *beta,
	  double *gr,
	  void *ex){

    int i, j, k;
    int start;

    double loglik;
    double *hessian = NULL;

    Exts *ext;
    Family *fam;

    double post_mode;

    int level = 1;

    char trans = 'N';
    double alpha = 1.0;
    int one = 1;


/* Note that we here trust 'ext->x_beta' to be properly updated !! */
/* In 'fun' */
 
    ext = ex;

    fam = Calloc(1, Family);
    fam->x = Calloc(ext->p, double *);
    fam->p = ext->p;
    fam->m = 0;
    fam->k = 0;
    fam->sigma = beta[ext->p];

    loglik = 0.0;

    for (k = 0; k < pp1; k++){
	gr[k] = 0.0;
    }

/*
    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j][i]; 
	}
	ext->x_beta[i] = tmp;
    }
*/

    F77_CALL(dcopy)(&(ext->n), ext->offset, &one, ext->x_beta, &one);
    F77_CALL(dgemv)(&trans, &(ext->n), &(ext->p), 
    &alpha, ext->x[0], &(ext->n), beta, &one,
		    &alpha, ext->x_beta, &one);

  
    start = 0;

    for (i = 0; i < ext->n_fam; i++){
	fam->n = ext->fam_size[i];
	fam->cluster_weight = ext->cluster_weights[i];
	fam->x_beta = ext->x_beta + start;
	fam->yw = ext->yw + start;
	fam->weights = ext->weights + start;
	for (j = 0; j < ext->p; j++)
	    fam->x[j] = ext->x[j] + start; 
	update(level,
	       ext->p,
	       beta,
	       &loglik,
	       gr,
	       hessian,
	       &post_mode,
	       fam,
	       ext->n_points,
	       ext->wc,
	       ext->zeros);
	ext->post_mode[i] = post_mode;
	start += ext->fam_size[i];
    }

    for (i = 0; i < pp1; i++){
	gr[i] = -gr[i]; /* Minimization! */
/*
	Rprintf("gr[%d] = %f   ", i, gr[i]);
	Rprintf("beta[%d] = %f\n", i, beta[i]);
*/
    }

    Free(fam->x);
    Free(fam);
}

void fun2(int pp1, 
	  double *beta,
	  double *loglik,
	  double *gr,
	  double *hessian,
	  void *ex){

    int i, j, k;
    int start;
    double post_mode;

    Exts *ext;
    Family *fam;

    int level = 2;

    char trans = 'N';
    double alpha = 1.0;
    int one = 1;

    ext = ex;

    fam = Calloc(1, Family);
    fam->x = Calloc(ext->p, double *);
    fam->p = ext->p;
    fam->m = 0;
    fam->k = 0;
    fam->sigma = beta[ext->p];

    *loglik = 0.0;

    for (k = 0; k < pp1; k++){
	gr[k] = 0.0;
    }
    
    for (i = 0; i < pp1 * pp1; i++){
	hessian[i] = 0.0;
    }
/*
    for (i = 0; i < ext->n; i++){
	tmp = ext->offset[i]; 
	for (j = 0; j < ext->p; j++){
	    tmp += beta[j] * ext->x[j][i];
	}
	ext->x_beta[i] = tmp;
    }
*/
    F77_CALL(dcopy)(&(ext->n), ext->offset, &one, ext->x_beta, &one);
    F77_CALL(dgemv)(&trans, &(ext->n), &(ext->p), 
    &alpha, ext->x[0], &(ext->n), beta, &one,
		    &alpha, ext->x_beta, &one);

    
    start = 0;
    for (i = 0; i < ext->n_fam; i++){
	fam->n = ext->fam_size[i];
	fam->cluster_weight = ext->cluster_weights[i];
	fam->x_beta = ext->x_beta + start;
	fam->yw = ext->yw + start;
	fam->weights = ext->weights + start;
	for (j = 0; j < ext->p; j++)
	    fam->x[j] = ext->x[j] + start; /* check this! */

	update(level,
	       ext->p,
	       beta,
	       loglik,
	       gr,
	       hessian,
	       &post_mode,
	       fam,
	       ext->n_points,
	       ext->wc,
	       ext->zeros);
	ext->post_mode[i] = post_mode;
        start += ext->fam_size[i];
    }

    for (i = 0; i < pp1 * pp1; i++){
	hessian[i] = -hessian[i];
    }
    Free(fam->x);
    Free(fam);
}
