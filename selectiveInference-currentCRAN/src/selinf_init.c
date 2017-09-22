#include <R.h>
#include <Rinternals.h>

/* . entry points */
extern void update1(double *Q2, double *w, int *mp, int *kp);
static R_NativePrimitiveArgType update1_t[] = {
  REALSXP, REALSXP, INTSXP, INTSXP
};

extern void downdate1(double *Q1, double *R, int *j0p, int *mp, int *np);
static R_NativePrimitiveArgType downdate1_t[] = {
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

static const R_CMethodDef CEntries[] = {
  {"update1", (DL_FUNC) &update1, 4},
  {"downdate1", (DL_FUNC) &downdate1, 5},
  {NULL, NULL, 0}
};

void R_init_cubature(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
