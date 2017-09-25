#include <R.h>

static void Sample_wtr(int k, int n, int *y)
{
    /* k-out-of-n without replacement. Result in  y. */
  
    int *x;
    int i, j;
  
  
    x = Calloc(n, int);
    
    for (i = 0; i < n; i++)
	x[i] = y[i]; /* Note change!! */
    
    for (i = 0; i < k; i++) {
	j = n * unif_rand();
	/* Rprintf("n = %d::j = %d\n",n, j); */
	y[i] = x[j]; /* Note: '+ 1' removed */ 
	x[j] = x[--n];
    }
    Free(x);
}

static void riskset_fill(int p_start, 
			 int *nn, double *enter, double *exit, int *event,
			 int *antrs, double *risktimes, 
			 int *n_events,int *size, int *eventset, 
			 int *riskset, int *max_s){
    /** nn = stratum size,
	enter[nn], exit[nn], event[nn] as usual
	antrs = No. of risksets in this stratum.
	risktimes[nn] (e.g. risksets[antrs])
	n.events[nn], size[nn] (eg [antrs] )!
    **/

    /* Data sorted ascending wrt exit, descending wrt event (for tied exit)

    */

    int i, j, start, nextstart, eindx, sindx, survs;
    double th;
    int *rs;

    start = 0;
    eindx = 0;
    sindx = 0;
    j = -1;
    while (start < *nn){
	/* Changed order of conditions in 2.2-3: */
	/* for (nextstart = start; (nextstart < *nn) & (event[nextstart] == 0); 
	   nextstart++); */
	nextstart = start;
	while (nextstart < *nn){
	    if (event[nextstart] == 1) break;
	    nextstart++;
	}
	if (nextstart >= *nn) return; /* Done in this stratum! */
    	j++;
	th = exit[nextstart];
    /* Changed order of conditions in 2.2-3: */
	/* for (start = nextstart; (start < *nn) & (exit[start] == th) & 
	   (event[start] == 1); start++){ */
	start = nextstart;
	while (start < *nn){
	    if ((exit[start] == th) & (event[start] == 1)){
		eventset[eindx] = start + p_start + 1; 
                /* +1 because C vs. R... */
		riskset[sindx] = start + p_start + 1;
		eindx++;
		sindx++;
	    }else{
		break;
	    }
	    start++;
	}
    
	for (i = start; i < *nn; i++){
	    if (enter[i] < th){
		riskset[sindx] = i + p_start + 1;
		sindx++;
	    }
	}
	survs = size[j] - n_events[j];
    
	if (survs > *max_s){
	    GetRNGstate();
	    rs = riskset + (sindx - survs); /* Back to first survivor. */
	    Sample_wtr(*max_s, survs, rs);
	    sindx += *max_s - survs;
	    size[j] = *max_s + n_events[j];
	    PutRNGstate();
	}
    }
}

void risk_get(int *max_s, int *nn, int *ns,
	      double *enter, double *exit, int*event, 
	      int* nstra, int *l_nstra,
	      int *new_totrs,
	      int *antrs, int *n_events, int *size,
	      double *risktimes,
	      int *eventset, int *riskset){

    int j, stra, start, stopp, rsindx, eindx, sindx, l_nn;
  
    rsindx = 0;
    eindx = 0;
    sindx = 0;

    for (stra = 0; stra < *l_nstra - 1; stra++){
	start = nstra[stra];
	stopp = nstra[stra + 1];
	l_nn = stopp - start;
	if (antrs[stra] > 0){
	    riskset_fill(start, &l_nn, 
			 enter + start, exit + start, event + start,
			 &(antrs[stra]), risktimes + rsindx,
			 n_events + rsindx, size + rsindx, eventset + eindx, 
			 riskset + sindx, max_s);
	    for (j = 0; j < antrs[stra]; j++){
		eindx += n_events[rsindx + j];
		sindx += size[rsindx + j];
	    }
	    rsindx += antrs[stra];
	}
    }
    *new_totrs = sindx;
}
