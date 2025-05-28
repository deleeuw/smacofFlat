#include "aplinline.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define CDEF(name, n)  {#name, (DL_FUNC) &name, n, {INTSXP, INTSXP, INTSXP, INTSXP}}

R_CallMethodDef callMethods[]  = {
	CALLDEF(APLDECODE, 2),
	CALLDEF(APLENCODE, 2),
	CALLDEF(APLSELECT, 3),
	CALLDEF(APLTRANSPOSE, 5),
	CALLDEF(APLSCAN, 5),
	CALLDEF(APLREDUCE, 6),
	CALLDEF(APLINNERPRODUCT, 9),
	{NULL, NULL, 0}
};
/*
R_CMethodDef cMethods[] = {
	CDEF(apldecode, 4),
	CDEF(aplencode, 4),
	{NULL, NULL, 0}
};
*/
//{"apldecode", (DL_FUNC) &apldecode, 4, {INTSXP, INTSXP, INTSXP, INTSXP}},
//{"aplencode", (DL_FUNC) &aplencode, 4, {INTSXP, INTSXP, INTSXP, INTSXP}},

/* Initializer for apl called upon loading the package.*/
void
R_init_apl(DllInfo *info)
{
	//R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	//R_RegisterCCallable("apl", "APLDECODE",  APLDECODE);
	//R_RegisterCCallable("apl", "APLENCODE", APLENCODE);
	//R_RegisterCCallable("apl", "APLSELECT",  APLSELECT);
	//R_RegisterCCallable("apl", "APLTRANSPOSE",  APLTRANSPOSE);
	//R_RegisterCCallable("apl", "APLSCAN", APLSCAN);
	//R_RegisterCCallable("apl", "APLREDUCE", APLREDUCE);
	//R_RegisterCCallable("apl", "APLINNERPRODUCT", APLINNERPRODUCT);
	//R_RegisterCCallable("apl", "apldecode", apldecode);
	//R_RegisterCCallable("apl", "aplencode", aplencode);
}



/* Finalizer for apl called upon unloading the package.
void R_unload_apl(DllInfo *info){
    M_cholmod_finish(&c);
}
 */
