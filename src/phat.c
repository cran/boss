#include "R.h"
#include "Rinternals.h"
#include "Rmath.h" 	
#include "math.h"
#include "R_ext/Lapack.h"
#include "R_ext/BLAS.h"

SEXP phatsum(SEXP listX, SEXP listr) {

R_len_t n=LENGTH(listX);
R_len_t i;

double  one  =  1.0; 
double  zero  =  0.0; 
SEXP xtrxtrtsum, xtr, xtrxtrt;

SEXP* xtrxtrti;
xtrxtrti=(SEXP *) malloc(n*sizeof(SEXP));



SEXP X, r;
int *dimX, *dimr;

X = VECTOR_ELT(listX,0);
r = VECTOR_ELT(listr,0);

//PrintValue(listX);
//PrintValue(listr);

dimX = INTEGER(getAttrib(X, R_DimSymbol));
dimr = INTEGER(getAttrib(r, R_DimSymbol));

PROTECT(xtrxtrtsum = allocMatrix(REALSXP, dimX[1],dimX[1]));
PROTECT(xtr = allocMatrix(REALSXP, dimX[1],dimr[1]));
PROTECT(xtrxtrt = allocVector(VECSXP, n+1));

F77_CALL(dgemm)("T", "N", dimX+1, dimr+1, dimX+0, &one, REAL(X), dimX+0, REAL(r), dimr+0, &zero, REAL(xtr), dimX+1);

PROTECT(xtrxtrti[0] = allocMatrix(REALSXP, dimX[1],dimX[1]));	
F77_CALL(dgemm)("N", "T", dimX+1, dimX+1, dimr+1, &one, REAL(xtr), dimX+1, REAL(xtr), dimX+1, &zero, REAL(xtrxtrti[0]), dimX+1);
SET_VECTOR_ELT(xtrxtrt, 0, xtrxtrti[0]);

F77_CALL(dgemm)("N", "T", dimX+1, dimX+1, dimr+1, &one, REAL(xtr), dimX+1, REAL(xtr), dimX+1, &zero, REAL(xtrxtrtsum), dimX+1);



for(i = 1; i < n; i++) {

X = VECTOR_ELT(listX,i);
r = VECTOR_ELT(listr,i);

dimX = INTEGER(getAttrib(X, R_DimSymbol));
dimr = INTEGER(getAttrib(r, R_DimSymbol));

F77_CALL(dgemm)("T", "N", dimX+1, dimr+1, dimX+0, &one, REAL(X), dimX+0, REAL(r), dimr+0, &zero, REAL(xtr), dimX+1);

PROTECT(xtrxtrti[i] = allocMatrix(REALSXP, dimX[1],dimX[1]));	
F77_CALL(dgemm)("N", "T", dimX+1, dimX+1, dimr+1, &one, REAL(xtr), dimX+1, REAL(xtr), dimX+1, &zero, REAL(xtrxtrti[i]), dimX+1);
SET_VECTOR_ELT(xtrxtrt, i, xtrxtrti[i]);

F77_CALL(dgemm)("N", "T", dimX+1, dimX+1, dimr+1, &one, REAL(xtr), dimX+1, REAL(xtr), dimX+1, &one, REAL(xtrxtrtsum), dimX+1);
}


SET_VECTOR_ELT(xtrxtrt, n, xtrxtrtsum);



UNPROTECT(3+n); 
free(xtrxtrti);

return xtrxtrt;

}





