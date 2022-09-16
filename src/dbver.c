#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <ctype.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("lotri", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#include "dbver.h"

SEXP _nlmixr2libMd5() {
  SEXP ret = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(ret, 0, Rf_mkChar(__MD5__));
  UNPROTECT(1);
  return ret;
}


void R_init_nlmixr2lib(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_nlmixr2libMd5", (DL_FUNC) &_nlmixr2libMd5, 0},
    {NULL, NULL, 0}
  };
  static const R_CMethodDef cMethods[] = {
    {NULL, NULL, 0, NULL}
  };
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_RegisterCCallable("nlmixr2lib", "_nlmixr2libMd5", (DL_FUNC) _nlmixr2libMd5);
  R_useDynamicSymbols(info, FALSE);
}
