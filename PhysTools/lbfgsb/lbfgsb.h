
#ifndef PHYSTOOLS_LBFGSB_H_INCLUDED
#define PHYSTOOLS_LBFGSB_H_INCLUDED

namespace phys_tools{
namespace lbfgsb{

///\cond 0
//hide from doxygen
	
typedef int integer;
typedef float real;
typedef double doublereal;
typedef int logical;
typedef short ftnlen;

/* Subroutine */ int setulb_(integer *n, integer *m, doublereal *x, 
	doublereal *l, doublereal *u, integer *nbd, doublereal *scale,
	doublereal *f, doublereal *g, doublereal *factr, doublereal *pgtol,
	doublereal *wa, integer *iwa, char *task, integer *iprint,
	char *csave, logical *lsave, integer *isave, doublereal *dsave);

///\endcond

} //namespace lbfgsb
} //namespace phys_tools

#include "interface.h"

#endif
