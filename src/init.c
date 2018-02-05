#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _SMED_SMED_selectC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SMED_SMED_selectYC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_SMED_SMED_selectC",  (DL_FUNC) &_SMED_SMED_selectC,  5},
  {"_SMED_SMED_selectYC", (DL_FUNC) &_SMED_SMED_selectYC, 6},
  {NULL, NULL, 0}
};

void R_init_SMED(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
