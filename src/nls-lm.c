#include <R.h>
#include <Rdefines.h>


typedef struct opt_struct {
    SEXP par;
    SEXP fcall;
    SEXP jcall;
    SEXP env;
    double ftol;
    double xtol;
    double gtol;
    double epsfcn;
    double *diag;
    double factor;
} opt_struct, *OptStruct;

extern void F77_NAME(lmdif)();
extern void F77_NAME(lmder)();
extern void fcn_lmdif(int *m, int *n, double *par, double *fvec,
                      int *iflag);
extern void fcn_lmder(int *m, int *n, double *par, double *fvec, double *fjac, int *ldfjac,
                      int *iflag);
extern char *fcn_message(char*, int, int);

int niter;
OptStruct OS;


static SEXP getListElement(SEXP list, char *str)
{
    SEXP elmt = NULL_USER_OBJECT, names = GET_NAMES(list);
    int i;

    for (i = 0; i < length(list); i++)
        if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    return elmt;
}

SEXP nls_lm(SEXP par_arg, SEXP fn, SEXP jac, SEXP control, SEXP rho)
{
    int     i, j, k;
    int     n, m, ldfjac;
    int     info, nfev, njev;

    double  *par, *fvec, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4,
              *r, *r2, *r2_x_perm, *hess;
    int     *ipvt, *perm;

    SEXP    eval_test;
    SEXP    sexp_diag, sexp_hess, sexp_fvec, sexp_info, sexp_message;
    SEXP    out, out_names;

    char    message[256];

    int     maxfev, nprint;
    int     mode;
    PROTECT_INDEX ipx;


    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));

    PROTECT(OS->par = AS_VECTOR(duplicate(par_arg)));
    n = length(OS->par);
    for (i = 0; i < n; i++)
        SET_VECTOR_ELT(OS->par, i, AS_NUMERIC(VECTOR_ELT(OS->par, i)));

    if (!isFunction(fn)) error("fn is not a function!");
    PROTECT(OS->fcall = lang2(fn, OS->par));

    if (!isEnvironment(rho)) error("rho is not an environment!");
    OS->env = rho;

    PROTECT(eval_test = eval(OS->fcall, OS->env));
    if (!IS_NUMERIC(eval_test))
        error("evaluation of fn function returns non-numeric vector!");
    UNPROTECT(1);

    m = length(eval_test);
    ldfjac = (m > n ? m : n);

    par       = (double*) R_alloc(n,          sizeof(double));
    fvec      = (double*) R_alloc(m,          sizeof(double));
    fjac      = (double*) R_alloc(ldfjac * n, sizeof(double));
    qtf       = (double*) R_alloc(n,          sizeof(double));
    wa1       = (double*) R_alloc(n,          sizeof(double));
    wa2       = (double*) R_alloc(n,          sizeof(double));
    wa3       = (double*) R_alloc(n,          sizeof(double));
    wa4       = (double*) R_alloc(m,          sizeof(double));
    ipvt      = (int   *) R_alloc(n,          sizeof(int   ));
    perm      = (int   *) R_alloc(n * n,      sizeof(int   ));
    r         = (double*) R_alloc(n * n,      sizeof(double));
    r2        = (double*) R_alloc(n * n,      sizeof(double));
    r2_x_perm = (double*) R_alloc(n * n,      sizeof(double));
    hess      = (double*) R_alloc(n * n,      sizeof(double));

    OS->ftol   = NUMERIC_VALUE(getListElement(control, "ftol"));
    OS->xtol   = NUMERIC_VALUE(getListElement(control, "ptol"));
    OS->gtol   = NUMERIC_VALUE(getListElement(control, "gtol"));
    OS->epsfcn = NUMERIC_VALUE(getListElement(control, "epsfcn"));
    OS->factor = NUMERIC_VALUE(getListElement(control, "factor"));
    OS->diag   = (double*) R_alloc(n, sizeof(double));

    PROTECT_WITH_INDEX(sexp_diag = getListElement(control, "diag"), &ipx);
    if (length(sexp_diag) == n) {
        REPROTECT(sexp_diag = AS_VECTOR(duplicate(sexp_diag)), ipx);
        for (i = 0; i < n; i++) {
            SET_VECTOR_ELT(sexp_diag, i, AS_NUMERIC(VECTOR_ELT(sexp_diag, i)));
            OS->diag[i] = NUMERIC_VALUE(VECTOR_ELT(sexp_diag, i));
        }
        mode = 2;
    }
    else {
        REPROTECT(sexp_diag = NEW_LIST(n), ipx);
        for (i = 0; i < n; i++)
            SET_VECTOR_ELT(sexp_diag, i, NEW_NUMERIC(1));
        mode = 1;
    }

    maxfev = INTEGER_VALUE(getListElement(control, "maxfev"));
    nprint = INTEGER_VALUE(getListElement(control, "nprint"));

    for (i = 0; i < n; i++)
        par[i] = NUMERIC_VALUE(VECTOR_ELT(OS->par, i));

    niter = -1;

/*========================================================================*/

    if (isNull(jac)) {
        F77_CALL(lmdif)(&fcn_lmdif, &m, &n, par, fvec,
                        &OS->ftol, &OS->xtol, &OS->gtol,
                        &maxfev, &OS->epsfcn, OS->diag, &mode,
                        &OS->factor, &nprint, &info, &nfev, fjac, &ldfjac,
                         ipvt, qtf, wa1, wa2, wa3, wa4);
    }
    else {
        if (!isFunction(jac))
            error("jac is not a function!");
        PROTECT(OS->jcall = lang2(jac, OS->par));
        PROTECT(eval_test = eval(OS->jcall, OS->env));
        if (!IS_NUMERIC(eval_test))
            error("evaluation of jac function returns non-numeric vector!");
        if (length(eval_test) != n*m)
            error("jac function must return numeric vector with length"
                  " == length(par) * length(fn(par, ...)). Your function"
                  " returns one with length %d while %d expected.",
                  length(eval_test), n*m);
        UNPROTECT(1);
        F77_CALL(lmder)(&fcn_lmder, &m, &n, par, fvec,
                         fjac, &ldfjac,
                        &OS->ftol, &OS->xtol, &OS->gtol,
                        &maxfev, OS->diag, &mode,
                        &OS->factor, &nprint, &info, &nfev, &njev,
                         ipvt, qtf, wa1, wa2, wa3, wa4);
        UNPROTECT(1);
    }

/*========================================================================*/

    fcn_message(message, info, maxfev);

    if (info < 1 || 8 < info)
        error("nls-lm: info = %d. %s\n\n", info, message);

    PROTECT(sexp_hess = NEW_NUMERIC(n*n));
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++) {
            perm[j*n + i] = (i + 1 == ipvt[j]) ? 1 : 0;
            r   [j*n + i] = (i <= j) ? fjac[j*ldfjac + i] : 0;
        }

    for (i = 0; i < n*n; i++) 
        r2[i] = r2_x_perm[i] = hess[i] = 0;

   /*  sexp_hess = sexp_perm %*% T(sexp_r) %*% sexp_r %*% T(sexp_perm)  *
    *              |             |________ r2_______|                |  *
    *              |             |_________r2_x_perm_________________|  *
    *              |_______________________hess______________________|  */

    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            for (k = 0; k < n; k++)
                r2[j*n + i] += r[i*n + k] * r[j*n + k];

    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            for (k = 0; k < n; k++)
                r2_x_perm[j*n + i] += r2[k*n + i] * perm[k*n + j];

    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            for (k = 0; k < n; k++)
                hess[j*n + i] += perm[k*n + i] * r2_x_perm[j*n + k];

    for (i = 0; i < n*n; i++)
        NUMERIC_POINTER(sexp_hess)[i] = hess[i];

    PROTECT(sexp_fvec = NEW_NUMERIC(m));
    for (i = 0; i < m; i++)
        NUMERIC_POINTER(sexp_fvec)[i] = fvec[i];

    PROTECT(sexp_info = NEW_INTEGER(1));
    INTEGER_POINTER(sexp_info)[0] = info;

    PROTECT(sexp_message = NEW_STRING(1));
    SET_VECTOR_ELT(sexp_message, 0, mkChar(message));

    for (i = 0; i < n; i++)
        NUMERIC_POINTER(VECTOR_ELT(sexp_diag, i))[0] = OS->diag[i];

    PROTECT(out = NEW_LIST(6));
    SET_VECTOR_ELT(out, 0, OS->par);
    SET_VECTOR_ELT(out, 1, sexp_hess);
    SET_VECTOR_ELT(out, 2, sexp_fvec);
    SET_VECTOR_ELT(out, 3, sexp_info);
    SET_VECTOR_ELT(out, 4, sexp_message);
    SET_VECTOR_ELT(out, 5, sexp_diag);

    PROTECT(out_names = NEW_STRING(6));
    SET_VECTOR_ELT(out_names, 0, mkChar("par"));
    SET_VECTOR_ELT(out_names, 1, mkChar("hessian"));
    SET_VECTOR_ELT(out_names, 2, mkChar("fvec"));
    SET_VECTOR_ELT(out_names, 3, mkChar("info"));
    SET_VECTOR_ELT(out_names, 4, mkChar("message"));
    SET_VECTOR_ELT(out_names, 5, mkChar("diag"));

    SET_NAMES(out, out_names);

    UNPROTECT(9);

    return out;
}
