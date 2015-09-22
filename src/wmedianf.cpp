// The software is made available under the MIT license:

//Copyright (c) 2015 Andre Rauh
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in
//all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//THE SOFTWARE.

// To create the DLL for MATLAB use.
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmex.lib")

#include <algorithm>
#include <cmath>
#include <limits>
#include "mex.h"
#include "wmedianf_impl.hpp"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  switch( nrhs ) {
    /* Note: Fall through is intentional */
    /* Argument 3: Specifies 1st pivot */
    case 3:
      if( mxIsComplex(prhs[2]) )
        mexErrMsgTxt("Third argument must be real.");
      if( mxGetNumberOfElements(prhs[2]) != 1 )
        mexErrMsgTxt("Third argument must be a scalar.");
      if( mxGetClassID(prhs[2]) != mxGetClassID(prhs[1]) )
        mexErrMsgTxt("Third argument must be of same type as second argument.");


    /* Argument 2: Either weights or first pivot (for median) */
    case 2:
      if( mxIsComplex(prhs[1]) )
        mexErrMsgTxt("Seond argument must be real.");
      if( mxGetClassID(prhs[1]) != mxGetClassID(prhs[0]) )
        mexErrMsgTxt("Second argument must be of same type as first argument.");
      {
        const mwSize* dim1 = mxGetDimensions(prhs[1]);
        if( ((dim1[0] != 1) && (dim1[1] != 1)) && (mxGetNumberOfElements(prhs[1]) != 1) ) {
          mexErrMsgTxt("Second argument must be a vector (weights) or scalar (pivot).");
        }
      }

    /* Argument 1: Data (N elements) */
    case 1:
      if( mxIsComplex(prhs[0]) )
        mexErrMsgTxt("wmedianf does only support real arguments");
      {
        const mwSize* dim0 = mxGetDimensions(prhs[0]);
        if( (dim0[0] != 1) && (dim0[1] != 1) ) {
          mexErrMsgTxt("First argument must be a vector");
        }
      }

      break;
    default:
      mexErrMsgTxt("wmedianf needs at least 1 but no more than 3 arguments.");
  }

  if( nlhs <= 0 )
    mexErrMsgTxt("wmedianf must have at least one return element");

#if _DEBUG
  // Reset counter of comparisons
  int num_comparisons = 0;
#endif
  mwSize N = mxGetNumberOfElements(prhs[0]);
  mwSize NW = 1; // Number weights
  if( nrhs >= 2 )
    NW = mxGetNumberOfElements(prhs[1]);

  if( (NW != N) && (NW != 1) )
    mexErrMsgTxt("Number of weights must either be a scalar or same size as x.");

  switch( mxGetClassID(prhs[0]) )
  {
    case mxDOUBLE_CLASS: {
      /* Pointer for the real part of the input */
      mxArray* x_mx = mxDuplicateArray(prhs[0]);
      double* x = (double*)mxGetData(x_mx);
      double* pivot = 0;

      // If there are 3 arguments the last on is the pivot
      if( nrhs >= 3 )
        pivot = (double*)mxGetData(prhs[2]);
      // If there are 2 arguments and the second one is a scalar this is the pivot then.
      if( (nrhs == 2) && (mxGetNumberOfElements(prhs[1]) == 1) )
        pivot = (double*)mxGetData(prhs[1]);

      plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
      double* y = (double*)mxGetData(plhs[0]);

      if( NW > 1 ) {
        // Compute weighted median
        double* w = (double*)mxGetData(mxDuplicateArray(prhs[1]));
#if _DEBUG
        *y = wmedianf(x, w, pivot, (int)N, -1.0, num_comparisons);
#else
        *y = wmedianf(x, w, pivot, (int)N, -1.0);
#endif
      } else {
        // Compute standard median.
#if _DEBUG
        *y = medianf(x, pivot, (int)N, num_comparisons);
#else
        *y = medianf(x, pivot, (int)N);
#endif
      }

#if _DEBUG
      if( nlhs >= 2 ) {
        plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        int* c = (int*)mxGetData(plhs[1]);
        *c = num_comparisons;
      }
#else
      if( nlhs >= 2 ) {
        mexErrMsgTxt("This is the release build and hence doesn't calculate number of comparisions.");
      }
#endif

      break;
    }
    case mxSINGLE_CLASS: // TODO
    case mxUINT32_CLASS: // TODO
    default:
      mexErrMsgTxt("Unsuporrted type.");
  }
}
