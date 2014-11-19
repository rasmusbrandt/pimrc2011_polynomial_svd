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
#define dgemm dgemm_
#endif

#include "mex.h"
#include "blas.h"

#define TYP_DIMS 3

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* matrix variables */
    mxArray *A, *B, *C;
    double *prA, *prB, *prC, *piA, *piB, *piC;
    unsigned int ciA, ciB, ciC;
    
    /* dimension variables */
    mwSize ndimsA, ndimsB;
    const mwSize *tmp_dims;
    mwSize dimsA[TYP_DIMS], dimsB[TYP_DIMS], dimsC[TYP_DIMS];
    mwSize el2A, el2B, el2C;
    
    /* blas variables */
    char *chn = "N"; double one = 1.0, minusone = -1.0;
    
    /* other variables */
    int i,j;
    double *dtmp;
    
    /* consistency checks */
    if (nrhs != 4) {
        mexErrMsgTxt("Pass two matrices, with const_ind's.");
    }
    
    /* set up matrices nicely (-1 for matlab index style) */
    A = prhs[0];
    dtmp = mxGetPr(prhs[1]); ciA = (unsigned int)dtmp[0] - 1; 
    B = prhs[2];
    dtmp = mxGetPr(prhs[3]); ciB = (unsigned int)dtmp[0] - 1;
    
    /* get input matrix dimensions */
    ndimsA = mxGetNumberOfDimensions(A);
    tmp_dims = mxGetDimensions(A);
    dimsA[0] = tmp_dims[0]; dimsA[1] = tmp_dims[1];
    if (ndimsA == 2) {
        dimsA[2] = 1;
    } else {
        dimsA[2] = tmp_dims[2];
    }
    el2A = dimsA[0] * dimsA[1];
    
    ndimsB = mxGetNumberOfDimensions(B);
    tmp_dims = mxGetDimensions(B);
    dimsB[0] = tmp_dims[0]; dimsB[1] = tmp_dims[1];
    if (ndimsB == 2) {
        dimsB[2] = 1;
    } else {
        dimsB[2] = tmp_dims[2];
    }
    el2B = dimsB[0] * dimsB[1];
    
    /* consistency check */
    if (dimsA[1] != dimsB[0]) {
        mexErrMsgTxt("Bad matrix dimensions for multiplication.");
    }
            
    /* set up output matrix dimensions */
    dimsC[0] = dimsA[0]; dimsC[1] = dimsB[1];
    dimsC[2] = dimsA[2] + dimsB[2] - 1;
    el2C = dimsC[0] * dimsC[1];
    
    /* set up output matrix */
    C = mxCreateNumericArray(TYP_DIMS, dimsC, mxDOUBLE_CLASS, mxCOMPLEX);
    ciC = ciA + ciB;
    
    /* get the pointers */
    prA = mxGetPr(A); piA = mxGetPi(A);
    prB = mxGetPr(B); piB = mxGetPi(B);
    prC = mxGetPr(C); piC = mxGetPi(C);
    
    /* multipy! */
    for(i = 0; i < dimsA[2]; i++) {
        for(j = 0; j < dimsB[2]; j++) {
            /* prA*prB */
            if ((prA != NULL) && (prB != NULL)) {
                dgemm(chn, chn,                         /* operators */
                      &dimsA[0], &dimsB[1], &dimsA[1],  /* dimensions */
                      &one,                             /* constant */
                      prA + i*el2A, &dimsA[0],          /* A matrix */
                      prB + j*el2B, &dimsB[0],          /* B matrix */
                      &one,                             /* constant */
                      prC + (i+j)*el2C, &dimsA[0]);     /* C matrix */
            }
            
            /* - piA*piB */
            if ((piA != NULL) && (piB != NULL)) {
                dgemm(chn, chn,                         /* operators */
                      &dimsA[0], &dimsB[1], &dimsA[1],  /* dimensions */
                      &minusone,                        /* constant */
                      piA + i*el2A, &dimsA[0],          /* A matrix */
                      piB + j*el2B, &dimsB[0],          /* B matrix */
                      &one,                             /* constant */
                      prC + (i+j)*el2C, &dimsA[0]);     /* C matrix */
            }
            
            /* j*prA*piB */
            if ((prA != NULL) && (piB != NULL)) {
                dgemm(chn, chn,                         /* operators */
                      &dimsA[0], &dimsB[1], &dimsA[1],  /* dimensions */
                      &one,                             /* constant */
                      prA + i*el2A, &dimsA[0],          /* A matrix */
                      piB + j*el2B, &dimsB[0],          /* B matrix */
                      &one,                             /* constant */
                      piC + (i+j)*el2C, &dimsA[0]);     /* C matrix */
            }
            
            /* j*piA*prB */
            if ((piA != NULL) && (prB != NULL)) {
                dgemm(chn, chn,                         /* operators */
                      &dimsA[0], &dimsB[1], &dimsA[1],  /* dimensions */
                      &one,                             /* constant */
                      piA + i*el2A, &dimsA[0],          /* A matrix */
                      prB + j*el2B, &dimsB[0],          /* B matrix */
                      &one,                             /* constant */
                      piC + (i+j)*el2C, &dimsA[0]);     /* C matrix */
            }
        }
    }
    
    /* output variables */
    plhs[0] = C;
    plhs[1] = mxCreateNumericMatrix(1,1, mxDOUBLE_CLASS, mxREAL);
    dtmp = mxGetPr(plhs[1]);
    dtmp[0] = ciC + 1;          /* +1 for MATLAB index... */
}
