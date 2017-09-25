#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>

#include "phfun.h"

extern ph0S_fun *S0;
extern ph0_fun *f0;
extern ph0_fun *h0;
extern ph0_fun *f0_t;
extern ph0_fun *h0_t;
extern ph0_fun *h0_tt;

/* Weibull (exponential, really) distribution */

double S0_weibull(double x, int log_p){
    double scale = 1.0;
    int lower_tail = 0;
    
    double ret;

/* NOTE: Exponential distribution with rate = scale = 1 */

    ret = pexp(x, scale, lower_tail, log_p);
    
    return(ret);
}


double f0_weibull(double x){
    double scale = 1.0;
    int give_log = 0;
    double ret;

    ret = dexp(x, scale, give_log);
    
    return(ret);
}

double h0_weibull(double x){
    double ret;

    if (x < 0.0){
	ret = 0.0;
    }else{
	ret = 1.0;
    }

    return(ret);
}

double f0_t_weibull(double x){
    double scale = 1.0;
    int give_log = 0;
    double ret;

    ret = -dexp(x, scale, give_log);
    
    return(ret);
}

double h0_t_weibull(double x){
    double ret;

    ret = 0.0;
    
    return(ret);
}

double h0_tt_weibull(double x){
	double ret;

    ret = 0.0;
    
    return(ret);
}

/* Log-logistic distribution */

double S0_loglogistic(double x, int log_p){
	double ret;

	if (log_p){
		if (x <= 0){
			ret = 0.0;
		}else{ 
			ret = -log1p(x);
		}
	}else{
		if (x <= 0.0){
			ret = 1.0;
		}else{
			ret = 1.0 / (1.0 + x);
		}

	}


    return(ret);
}

double f0_loglogistic(double x){
    double ret;

    if (x < 0){
	ret = 0;
    }else{
	ret = 1.0 / R_pow_di((1 + x), 2);
    }    
    
    return(ret);
}

double h0_loglogistic(double x){
    double ret;

    if (x < 0.0){
	ret = 0.0;
    }else{
	ret = 1.0 / (1.0 + x);
    }

    return(ret);
}

double f0_t_loglogistic(double x){
    double ret;

    if(x < 0){
	ret = 0.0;
    }else{
	ret = -2.0 / R_pow_di((1.0 + x), 3);
    }
    
    return(ret);
}

double h0_t_loglogistic(double x){
    double ret;

    if (x < 0.0){
	ret = 0.0;
    }else{
	ret = -1.0 / R_pow_di((1.0 + x), 2);
    }
    
    return(ret);
}

double h0_tt_loglogistic(double x){
    double ret;

    if (x < 0.0){
	ret = 0.0;
    }else{
	ret = 2.0 / R_pow_di((1.0 + x), 3);
    }
    
    return(ret);
}

/* Lognormal distribution */

double S0_lognormal(double x, int log_p){
    double logmean = 0.0;
    double logsd = 1.0;
    int lower_tail = 0;

    double ret;

    if (x <= 0.0){
      if (log_p){
	ret = 0.0;
      }else{
	ret = 1.0;
      }
    }else{
      ret = plnorm(x, logmean, logsd, lower_tail, log_p);
    }

    return(ret);
}

double f0_lognormal(double x){
    double logmean = 0.0;
    double logsd = 1.0;
    int give_log = 0;
    double ret;

    ret = dlnorm(x, logmean, logsd, give_log);
    
    return(ret);
}

double h0_lognormal(double x){
    double logmean = 0.0;
    double logsd = 1.0;
    int lower_tail = 0;
    int give_log = 0;
    double ret;

    if (x < 0.0){
	ret = 0.0;
    }else{
	ret = dlnorm(x, logmean, logsd, give_log) / 
	plnorm(x, logmean, logsd, lower_tail, give_log);
    }

    return(ret);
}

double f0_t_lognormal(double x){
    double ret;

    if (x <= 0.0){
	ret = 0.0;
    }else{
	ret = -f0_lognormal(x) * (1.0 + log(x)) / x;
    }
    
    return(ret);
}

/* NOT needed!? Oh yes!! See h0_tt_lognormal!! */
double f0_tt_lognormal(double x){
    double ret;
    
    if (x <= 0.0){
	ret = 0.0;
    }else{
	ret = -f0_t_lognormal(x) * (1.0 + log(x))/ x + 
	    f0_lognormal(x) * log(x) / R_pow_di(x, 2);
    }
    
    return(ret);
}


double h0_t_lognormal(double x){
    double ret;
    int log_p = 0;

    if (x <= 0){
	ret = 0.0;
    }else{
	ret = f0_t_lognormal(x) / S0_lognormal(x, log_p) +
	    R_pow_di(h0_lognormal(x), 2);
    }
    
    return(ret);
}

double h0_tt_lognormal(double x){
    double ret;
    int log_p = 0;

    if (x <= 0){
	ret = 0.0;
    }else{
	ret = f0_tt_lognormal(x) / S0_lognormal(x, log_p) +
	    f0_lognormal(x) * f0_t_lognormal(x) / 
	    R_pow_di(S0_lognormal(x, log_p), 2) +
	    2.0 * h0_t_lognormal(x) * h0_lognormal(x);
    }
    
    return(ret);
}

/* Ev distribution */
/* Note that a constant column in the design matrix is needed */

double S0_ev( double x, int log_p){
    double ret;

    if (x <= 0){
	ret = 0.0;
    }else{
	ret = -expm1(x);
    }
    if (!log_p) ret = exp(ret);
    return ( ret );
}

double f0_ev(double x){
    double ret;
    
    if (x < 0){
	ret = 0.0;
    }else{
	ret = exp(x - expm1(x));
    }

    return ( ret );
}

double h0_ev(double x){
    double ret;

    if (x < 0){
	ret = 0.0;
    }else{
	ret = exp(x);
    }

    return ( ret );
}

double f0_t_ev(double(x)){
    double ret, tmp;

    if (x < 0){
	ret = 0.0;
    }else{
	tmp = expm1(x);
	ret = -tmp * exp(x - tmp);
    }

    return ( ret );
}

/* Not needed!?
double f0_tt_ev(double x){
    double ret, tmp, expx;

    if (x < 0){
	ret = 0.0;
    }else{
	tmp = expm1(x);
	expx = exp(x);
	ret = expx * exp(-tmp) * (R_pow_di(tmp, 2) - expx);
    }

    return ( ret );
}
*/

double h0_t_ev(double x){
    double ret;

    if (x < 0){
	ret = 0.0;
    }else{
	ret = exp(x);
    }

    return ( ret );
}

double h0_tt_ev(double x){
    double ret;

    if (x < 0){
	ret = 0.0;
    }else{
	ret = exp(x);
    }

    return ( ret );
}

    

/********* The 'g' function and its partial derivatives: ***********/

double g(double time, double gam, double alpha){
    /* The "shape/scale" family: */
    double ret;
 
    if (time < 0.0){
	ret = -99.0;
	error("Negative 'time' to 'g' not allowed");
    }else if (time == 0.0){
	ret = 0.0;
    }else{
	ret = R_pow((time / exp(alpha)), exp(gam));
    }

    return(ret);
}

double g_t(double time, double gam, double alpha){
    /* The "shape/scale" family: */
    double ret;
 
    if (time <= 0.0){
	ret = -99.0;
	error("Non-positive 'time' to 'g_t' not allowed");
    }else{
	ret = (exp(gam) / time) * g(time, gam, alpha);
    }

    return(ret);
}

double g_gamma(double time, double gam, double alpha){
    /* The "shape/scale" family: */
    double ret;
 
    ret = g(time, gam, alpha);

    if (time > 0.0){
	ret = ret * log(ret);
    }

    return(ret);
}

double g_alpha(double time, double gam, double alpha){
    double ret;

    ret = -exp(gam) * g(time, gam, alpha);

    return(ret);
}

double g_t_gamma(double time, double gam, double alpha){
    double ret; 

    if (time <= 0.0){
	ret = -99.0;
	error("'time' must be strictly positive in 'g_t_gamma");
    }else{
	ret = g_t(time, gam, alpha) +
	    (exp(gam) / time) * g_gamma(time, gam, alpha);
    }

    return(ret);
}

double g_t_alpha(double time, double gam, double alpha){
    double ret;

    if (time <= 0.0){
	ret = -99.0;
	error("'time' must be strictly positive in 'g_t_gam");
    }else{
	ret = -exp(gam) * g_t(time, gam, alpha);
    }

    return(ret);
}

double g_gamma2(double time, double gam, double alpha){
    double ret;

    ret = g_gamma(time, gam, alpha);
    if (time > 0.0){
	ret = ret * (1.0 + log(g(time, gam, alpha)));
    }

    return(ret);
}

double g_gamma_alpha(double time, double gam, double alpha){
    double ret;

    ret = g_alpha(time, gam, alpha);

    if (time > 0.0){
	ret = ret * (1.0 + log(g(time, gam, alpha)));
    }

    return(ret);
}

double g_alpha2(double time, double gam, double alpha){
    double ret;

    ret = -exp(gam) * g_alpha(time, gam, alpha);

    return(ret);
}

double g_t_gamma2(double time, double gam, double alpha){
    double ret;

    if (time <= 0){
	ret = -99.0;
	error("'time' must be strictly positive in 'g_t_gamma2'");
    }else{
	ret = g_t_gamma(time, gam, alpha) + (exp(gam) / time) *
	    g_gamma(time, gam, alpha) * (2.0 + log(g(time, gam, alpha)));
    }

    return(ret);
}

double g_t_gamma_alpha(double time, double gam, double alpha){
    double ret;

    if (time <= 0.0){
	ret = -99.0;
	error("'time' must be strictly positive in 'g_t_gamma_alpha'");
    }else{
	ret =  -exp(gam) * (g_t(time, gam, alpha) + 
			    g_t_gamma(time, gam, alpha));
    }

    return(ret);
}

double g_t_alpha2(double time, double gam, double alpha){
    double ret;

    if (time <= 0.0){
	ret = -99.0;
	error("'time' must be strictly positive in 'g_t_alpha2'");
    }else{
	ret =  -exp(gam) * g_t_alpha(time, gam, alpha);
    }

    return(ret);
}

/************ And now the derived functions: **************/

double S(double x, double gam, double alpha, int log_p){
    double ret;

    ret = S0(g(x, gam, alpha), log_p);
    return(ret);
}

double f(double x, double gam, double alpha){
    double ret;
    
    ret = g_t(x, gam, alpha) * f0(g(x, gam, alpha));
    return(ret);
}

double S_gamma(double x, double gam, double alpha){
    return(-g_gamma(x, gam, alpha) * f0(g(x, gam, alpha)));
}

double S_alpha(double x, double gam, double alpha){
    double ret;

    ret = -g_alpha(x, gam, alpha) * f0( g(x, gam, alpha));

    return(ret);
}

double S_gamma2(double x, double gam, double alpha){
    return(-(g_gamma2(x, gam, alpha) * f0(g(x, gam, alpha)) +
	     R_pow_di(g_gamma(x, gam, alpha), 2) * 
	     f0_t(g(x, gam, alpha))));
}

double S_gamma_alpha(double x, double gam, double alpha){
    return(-(g_gamma_alpha(x, gam, alpha) * f0(g(x, gam, alpha)) +
	     g_gamma(x, gam, alpha) * g_alpha(x, gam, alpha) *
	     f0_t(g(x, gam, alpha))));
}

double S_alpha2(double x, double gam, double alpha){
    return(-(g_alpha2(x, gam, alpha) * f0(g(x, gam, alpha)) +
	     R_pow_di(g_alpha(x, gam, alpha), 2) * 
	     f0_t(g(x, gam, alpha))));
}

double h(double x, double gam, double alpha){
    /* The hazard function */
    double ret;

    ret = g_t(x, gam, alpha) * h0(g(x, gam, alpha));
    /* return(f(x, gam, alpha) / S(x, gam, alpha)); */
    return(ret);
}

double h_gamma(double x, double gam, double alpha){
    /* d/dgam h */
    double y;

    y = g(x, gam, alpha);
    return(g_t_gamma(x, gam, alpha) * h0(y) +
	   g_t(x, gam, alpha) * g_gamma(x, gam, alpha) * h0_t(y));
}
    
double h_alpha(double x, double gam, double alpha){
    /* d/dgam h */
    double y;

    y = g(x, gam, alpha);
    return(g_t_alpha(x, gam, alpha) * h0(y) +
	   g_t(x, gam, alpha) * g_alpha(x, gam, alpha) * h0_t(y));
}

double h_gamma2(double x, double gam, double alpha){
    double y, r1, r2, r3;

    y = g(x, gam, alpha);

    r1 = h0(y) * g_t_gamma2(x, gam, alpha);
    r2 = h0_t(y) * (g_t(x, gam, alpha) * g_gamma2(x, gam, alpha) +
		    g_t_gamma(x, gam, alpha) * g_gamma(x, gam, alpha) +
		    g_t_gamma(x, gam, alpha) * g_gamma(x, gam, alpha));
    r3 = h0_tt(y) * g_t(x, gam, alpha) *
	g_gamma(x, gam, alpha) * g_gamma(x, gam, alpha);

    return(r1 + r2 + r3);
}

double h_gamma_alpha(double x, double gam, double alpha){
    double y, r1, r2, r3;

    y = g(x, gam, alpha);
    r1 = h0(y) * g_t_gamma_alpha(x, gam, alpha);
    r2 = h0_t(y) * (g_t(x, gam, alpha) * g_gamma_alpha(x, gam, alpha) +
		    g_t_gamma(x, gam, alpha) * g_alpha(x, gam, alpha) +
		    g_t_alpha(x, gam, alpha) * g_gamma(x, gam, alpha));
    r3 = h0_tt(y) * g_t(x, gam, alpha) *
	g_gamma(x, gam, alpha) * g_alpha(x, gam, alpha);

    return(r1 + r2 + r3);
}

double h_alpha2(double x, double gam, double alpha){
    double y, r1, r2, r3;

    y = g(x, gam, alpha);
    r1 = h0(y) * g_t_alpha2(x, gam, alpha);
    r2 = h0_t(y) * (g_t(x, gam, alpha) * g_alpha2(x, gam, alpha) +
		    g_t_alpha(x, gam, alpha) * g_alpha(x, gam, alpha) +
		    g_t_alpha(x, gam, alpha) * g_alpha(x, gam, alpha));
    r3 = h0_tt(y) * g_t(x, gam, alpha) *
	g_alpha(x, gam, alpha) * g_alpha(x, gam, alpha);

    return(r1 + r2 + r3);
}
