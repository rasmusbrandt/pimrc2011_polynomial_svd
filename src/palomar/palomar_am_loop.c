// pimrc2011_polynomial_svd
// Copyright (C) 2011, Rasmus Brandt

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#if !defined(_WIN32)
#define ddot ddot_
#endif

#include "mex.h"
#include "blas.h"
#include <math.h>
                    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* matrix variables */
    mxArray *Q, *X, *Xin;
    double *prQ, *prX, *piQ, *piX;
    
    /* dimension variables */
    mwSize N;
    
    /* blas variables */
    char *chn = "N"; int one = 1;
    
    /* other variables */
    double acc_r, acc_i, phi;
    unsigned int n, nNp1;
    
    /* consistency checks */
    if (nrhs != 2) {
        mexErrMsgTxt("Incorrect number of inputs.");
    }
    
    /* set up matrices */
    Q = prhs[0];
    Xin = prhs[1];
    
    /* get input matrix dimensions */
    N = mxGetN(Q);
    
    /* consistency check */
    if (mxGetM(Q) != N) {
        mexErrMsgTxt("Q not square.");
    }
    if (N != mxGetM(Xin)) {
        mexErrMsgTxt("Incompatible dimensions for multiplication.");
    }
    if (mxGetN(Xin) != 1) {
        mexErrMsgTxt("X not a vector.");
    }
    if (!mxIsComplex(Q)) {
        mexErrMsgTxt("Q should be complex");
    }
    if (!mxIsComplex(Xin)) {
        mexErrMsgTxt("Please supply complex initial guess");
    }
    
    /* create output vector */
    X = mxCreateDoubleMatrix(N,1,mxCOMPLEX);
    
    /* get the pointers */
    prQ = mxGetPr(Q); piQ = mxGetPi(Q);
    prX = mxGetPr(X); piX = mxGetPi(X);
    
    /* copy */
    memcpy(prX, mxGetPr(Xin), N*sizeof(double));
    memcpy(piX, mxGetPi(Xin), N*sizeof(double));
    
    /* loop! */
    for(n = 0; n < N; n++) {
        /* Qtkk = Q(k,:)*x - Q(k,k)*x(k); */
        
        /* prQ*prX */
        acc_r = ddot(&N, prQ + n, &N, prX, &one);
        
        /* - piQ*piX */
        acc_r = acc_r - ddot(&N, piQ + n, &N, piX, &one);
        
        /* j*prQ*piX */
        acc_i = ddot(&N, prQ + n, &N, piX, &one);
        
        /* j*piQ*prX */
        acc_i = acc_i + ddot(&N, piQ + n, &N, prX, &one);
        
        /* - Q(k,k)*x(k) */
        nNp1 = n*(N+1);
        acc_r = acc_r - (*( prQ + nNp1 )) * (*( prX + n ))
                      + (*( piQ + nNp1 )) * (*( piX + n ));
        acc_i = acc_i - (*( prQ + nNp1 )) * (*( piX + n ))
                      - (*( piQ + nNp1 )) * (*( prX + n ));

        /* phi_k = angle(Qtkk) + pi; */
        phi = atan2(acc_i, acc_r) + M_PI;
        
        /* x(k) = exp(1i*phi_k); */
        *(prX + n) = cos(phi);
        *(piX + n) = sin(phi);
    }
    
    /* output variables */
    plhs[0] = X;
}
