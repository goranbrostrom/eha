#ifndef GB_ZEROIN
#define GB_ZEROIN

double GB_zeroin(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit);				/* Max # of iterations */

#endif
