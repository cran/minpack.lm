#include <R.h>
#include <Rdefines.h>
#include "minpack_lm.h"


void fcn_lmdif(int *m, int *n, double *par, double *fvec, int *iflag)
{
    int i;
    double sumf;
    SEXP sexp_fvec;

    /* Rprintf("fcn-lmdif calling...\n"); */
    if (IS_NUMERIC(OS->par))
      for (i = 0; i < *n; i++) {
	if(par[i] < NUMERIC_POINTER(OS->lower)[i])
	  par[i] = NUMERIC_POINTER(OS->lower)[i];
	if(par[i] > NUMERIC_POINTER(OS->upper)[i])
	  par[i] = NUMERIC_POINTER(OS->upper)[i];
	NUMERIC_POINTER(OS->par)[i] = par[i];
      }
    else
      for (i = 0; i < *n; i++) {
	if(par[i] < NUMERIC_POINTER(OS->lower)[i])
	  par[i] = NUMERIC_POINTER(OS->lower)[i];
	if(par[i] > NUMERIC_POINTER(OS->upper)[i])
	  par[i] = NUMERIC_POINTER(OS->upper)[i];
	NUMERIC_POINTER(VECTOR_ELT(OS->par, i))[0] = par[i];
      }
    if(*iflag == 0){ 
      if(OS->nprint > 0) {
	Rprintf("It. %4d, RSS = %10g, Par. =", OS->niter, 
		OS->rsstrace[OS->niter]);
	for (i = 0; i < *n; i++) 
	  Rprintf(" % 10g", par[i]);
	Rprintf("\n");
	
      }
      OS->niter++;
    }
    else if (*iflag == 1 || *iflag == 2) {
      SETCADR(OS->fcall, OS->par);
      PROTECT(sexp_fvec = eval(OS->fcall, OS->env));
      
      for (i = 0; i < *m; i++) 
	fvec[i] = NUMERIC_POINTER(sexp_fvec)[i];
      UNPROTECT(1);	
      sumf = 0;
      for (i = 0; i < *m; i++) 
	sumf += (fvec[i]*fvec[i]);
     
      OS->rsstrace[OS->niter] = sumf;
      
    }
    if(OS->niter == OS->maxiter)
      *iflag = -1;
    
}
