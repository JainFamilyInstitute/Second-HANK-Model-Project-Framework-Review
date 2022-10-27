#include "mex.h"
#include "blas.h"
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )


/*
DIAGPROD.C
This function calculates diagonal elements of matrix
product C=diag(A*B') (A,B - has same dimensions) and up to 2.5
times faster, then Matlab expression X=sum(A,*B,2).
(c) D.Tikhonov
*/

{
  double ddot();
  double *A, *B, *C, res;
  int m,n,m1,n1,i, one=1;

  A = mxGetPr(prhs[0]);
  B = mxGetPr(prhs[1]);
  m = mxGetM(prhs[0]); /* number of rows */
  n = mxGetN(prhs[0]); /* number of columns */


  n1 = mxGetN(prhs[1]);
  m1 = mxGetM(prhs[1]);

if (n != n1) {
    mexErrMsgTxt("Matrixs must have same dimensions");
  }


  plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
  C = mxGetPr(plhs[0]);


  /* Pass all arguments to Fortran by reference */

/* function ddot(n,cx,incx,cy,incy) */
for (i = 0; i < m; i++)
{
 C[i]=ddot(&n,A+i,&m,B+i,&m);
 }


}
