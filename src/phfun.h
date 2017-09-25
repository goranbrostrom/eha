#ifndef PHFUN_H
#define PHFUN_H

typedef double ph0S_fun(double, int);
typedef double ph0_fun(double);

/* Weibull (exponential, really) distribution */

double S0_weibull(double x, int log_p);

/* NOTE: Exponential distribution with rate = scale = 1 */

double f0_weibull(double x);

double h0_weibull(double x);

double f0_t_weibull(double x);

double h0_t_weibull(double x);

double h0_tt_weibull(double x);

/* Log-logistic distribution */

double S0_loglogistic(double x, int log_p);

double f0_loglogistic(double x);

double h0_loglogistic(double x);

double f0_t_loglogistic(double x);

double h0_t_loglogistic(double x);

double h0_tt_loglogistic(double x);

/* Lognormal distribution */

double S0_lognormal(double x, int log_p);

double f0_lognormal(double x);

double h0_lognormal(double x);

double f0_t_lognormal(double x);

double f0_tt_lognormal(double x); /* Unique for the lognormal! */

double h0_t_lognormal(double x);

double h0_tt_lognormal(double x);

/* Ev distribution */

double S0_ev(double x, int log_p);

double f0_ev(double x);

double h0_ev(double x);

double f0_t_ev(double x);

double h0_t_ev(double x);

double h0_tt_ev(double x);

/********* The 'g' function and its partial derivatives: ***********/

double g(double time, double gam, double alpha);

double g_t(double time, double gam, double alpha);

double g_gamma(double time, double gam, double alpha);

double g_alpha(double time, double gam, double alpha);

double g_t_gamma(double time, double gam, double alpha);

double g_t_alpha(double time, double gam, double alpha);

double g_gamma2(double time, double gam, double alpha);

double g_gamma_alpha(double time, double gam, double alpha);

double g_alpha2(double time, double gam, double alpha);

double g_t_gamma2(double time, double gam, double alpha);

double g_t_gamma_alpha(double time, double gam, double alpha);

double g_t_alpha2(double time, double gam, double alpha);


/************ And now the derived functions: **************/

double S(double x, double gam, double alpha, int log_p);

double f(double x, double gam, double alpha);

double S_gamma(double x, double gam, double alpha);

double S_alpha(double x, double gam, double alpha);

double S_gamma2(double x, double gam, double alpha);

double S_gamma_alpha(double x, double gam, double alpha);

double S_alpha2(double x, double gam, double alpha);

double h(double x, double gam, double alpha);

double h_gamma(double x, double gam, double alpha);
    
double h_alpha(double x, double gam, double alpha);

double h_gamma2(double x, double gam, double alpha);

double h_gamma_alpha(double x, double gam, double alpha);

double h_alpha2(double x, double gam, double alpha);

#endif
