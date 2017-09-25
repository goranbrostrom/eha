#include <R.h>

void strat_sizes(int *nn, double *enter, double *exit, int *event,
		 int *antrs, double *risktimes, 
		 int *n_events,int *size){
    /** nn = stratum size,
	enter[nn], exit[nn], event[nn] as usual
	antrs = No. of risksets in this stratum.
	risktimes[nn] (e.g. risksets[antrs])
	n.events[nn], size[nn] (eg [antrs] )!
    **/

    /* Data sorted ascending wrt exit, descending wrt event (for tied exit)

    */

    int i, start, nextstart;
    double th;

    for (i = 0; i < *nn; i++){
	n_events[i] = 0;
	size[i] = 0;
    }

    *antrs = 0;

    start = 0;

    while (start < *nn){
/* Reordered conditions in 2.2-3: */ 
	/* for (nextstart = start; (nextstart < *nn) & (event[nextstart] == 0); 
	   nextstart++); */
	nextstart = start;
	while (nextstart < *nn){
	    if (event[nextstart] == 1) break;
	    nextstart++;
	}
	if (nextstart >= *nn) return; /* Done in this stratum! */

	if (*antrs >= *nn) Rprintf("Error antrs in [sizes]\n");
	th = exit[nextstart];
	risktimes[*antrs] = th;
    
/* Reordered conditions in 2.2-3: */ 
	/* for (start = nextstart; (start < *nn) & (exit[start] == th) & 
	   (event[start] == 1); start++){ */
	start = nextstart;
	while (start < *nn){
	    if ((exit[start] == th) & (event[start] == 1)){
		n_events[*antrs]++;
		size[*antrs]++;
	    }else{
		break;
	    }
	    start++;
	}    
    
	for (i = start; i < *nn; i++){
	    if (enter[i] < th) size[*antrs]++;
	}
	(*antrs)++;
    }
}


void sizes(int *ns, int *nn, double *enter, double *exit, int *event,
	   int *antrs, int *nstra, double *risktimes, 
	   int *n_events,int *size, int *totrs){

    int j, start, s_nn, cum_antrs;

    cum_antrs = 0;
    for (j = 0; j < *ns; j++){
	start = nstra[j];
	s_nn = nstra[j + 1] - start;

	strat_sizes(&s_nn, enter + start, exit + start, event + start, 
		    antrs + j, risktimes + cum_antrs, 
		    n_events + cum_antrs, size + cum_antrs);
	cum_antrs += antrs[j];
    }
    *totrs = cum_antrs;

}
    
