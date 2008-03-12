#include <R.h>
#include <Rdefines.h>
#include "minpack_lm.h"


void fcn_lmdif(int *m, int *n, double *par, double *fvec, 
	       int *iflag, double *rss)
{
    int i;
    SEXP sexp_fvec;

    /* Rprintf("fcn-lmdif calling...\n"); */
    if (IS_NUMERIC(OS->par))
        for (i = 0; i < *n; i++) {
            if (!R_FINITE(par[i]))
                error("non-finite value supplied by lmdif!");
            NUMERIC_POINTER(OS->par)[i] = par[i];
        }
    else
        for (i = 0; i < *n; i++) {
            if (!R_FINITE(par[i]))
                error("non-finite value supplied by lmdif!");
            NUMERIC_POINTER(VECTOR_ELT(OS->par, i))[0] = par[i];
        }

    if      (*iflag == 0) {
      Rprintf("It. %4d, RSS = %10g, Par. =", niter, *rss);
      for (i = 0; i < *n; i++)
	Rprintf(" % 10g", par[i]);
      Rprintf("\n");
    }
    else if (*iflag == 1 || *iflag == 2) {
        SETCADR(OS->fcall, OS->par);
        PROTECT(sexp_fvec = eval(OS->fcall, OS->env));

        for (i = 0; i < *m; i++)
            fvec[i] = NUMERIC_POINTER(sexp_fvec)[i];

        UNPROTECT(1);

    }
    else if(*iflag == 3) niter++;
    
}
