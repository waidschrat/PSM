#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(lti_kalman_fulla_withinput)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(matexpfortransub)(void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"lti_kalman_fulla_withinput", (DL_FUNC) &F77_NAME(lti_kalman_fulla_withinput), 27},
  {"matexpfortransub",           (DL_FUNC) &F77_NAME(matexpfortransub),            5},
  {NULL, NULL, 0}
};

void R_init_PSM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}