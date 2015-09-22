// The software is made available under the MIT license:

//Copyright (c) 2010 Andre Rauh
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


#pragma once

#if _DEBUG
template <typename Ty> Ty medianf(Ty* x, Ty* pivot, int N, int& num_comparisons);
template <typename Ty> Ty quickSelect(Ty* x, int delta, int k, int N, int& num_comparisons);
double select_floyd_rivest(double *a, int l, int r, int k, int& num_comparison);

template <typename Ty> Ty wmedianf(Ty* x, Ty* w, Ty* pivot, int N, Ty W0, int& num_comparisons);
template <typename Ty> Ty wquickSelect(Ty* x, Ty* w, int N, Ty W0, int& num_comparisons, Ty Wl = 0.0, Ty Wr = 0.0);
template <typename Ty> Ty wmSort(Ty* x, Ty* w, int N, Ty W0);
#else
template <typename Ty> Ty medianf(Ty* x, Ty* pivot, int N);
template <typename Ty> Ty quickSelect(Ty* x, int delta, int k, int N);
double select_floyd_rivest(double *a, int l, int r, int k);

template <typename Ty> Ty wmedianf(Ty* x, Ty* w, Ty* pivot, int N, Ty W0);
template <typename Ty> Ty wquickSelect(Ty* x, Ty* w, int N, Ty W0, Ty Wl = 0.0, Ty Wr = 0.0);
template <typename Ty> Ty wmSort(Ty* x, Ty* w, int N, Ty W0);
#endif
