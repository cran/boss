#include "R.h"
#include "Rinternals.h"
#include "Rmath.h" 	
#include "math.h"
#include "R_ext/Lapack.h"
#include "R_ext/BLAS.h"

SEXP thatnew(SEXP list) {

R_len_t n=LENGTH(list);
R_len_t i;

double  one  =  1.0; 
double  zero  =  0.0; 
SEXP aat;
SEXP A;
int *dimA;

if(!isNewList(list)) error("`list' must be a list"); 

A = VECTOR_ELT(list,0);
//PrintValue(A);
//PrintValue(list);
dimA = INTEGER(getAttrib(A, R_DimSymbol));

PROTECT(aat = allocMatrix(REALSXP, dimA[0],dimA[0]));


F77_CALL(dgemm)("N", "T", dimA+0, dimA+0, dimA+1, &one, REAL(A), dimA+0, REAL(A), dimA+0, &zero, REAL(aat), dimA+0);

for(i = 1; i < n; i++) {
	F77_CALL(dgemm)("N", "T", dimA+0, dimA+0, dimA+1, &one, REAL(VECTOR_ELT(list,i)), dimA+0, REAL(VECTOR_ELT(list,i)), dimA+0, &one, REAL(aat), dimA+0);
}


UNPROTECT(1); 

return aat;

}





