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

#pragma once

#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include "wmedianf_impl.hpp"
#include "mathfunctions.hpp"

// Factor by which the left side is multiplied to yield the right side when
// using the bisection method to get the optimal k for a given M
#define BISECTION_DELTA 1.1
// Tolerance at which to stop when using the bisection method
#define BISECTION_TOL 1e-2
// If alph*M is smaller than this value the algorithm will use poisson
// distribution to model the cost function. Otherwise normal distribution will hold.
#define NORMAL_POISSON_LIMIT 10


#if _DEBUG
// Partitions the vector x around piv and sets the values for W0, Wl, Wr accordingly. If W0 is negative
// it will be computed. lrb specifies from which side of the array to count.
template <typename Ty> int wpartition(Ty* x, Ty* w, Ty piv, int N, int& lrb, Ty& W0, Ty& Wl, Ty& Wr, int& num_comparisons);
// A standard 2-way partitioning function
template <typename Ty> int  partition(Ty* x,        Ty piv, int N,          int& num_comparisons);
// Floyd & Rivest's select implementation with a delta to specify the inter element distance
double select_floyd_rivest_delta(double *a, int l, int r, int k, int delta, int& num_comparisons);
// Helper function to get the k th OS of a subset (equi-distant) of samples
template <typename Ty> Ty medianf_equidist(Ty* a, Ty* w, int l, int r, int k, int dist, int& num_comparisons);
#else
// Same prototypes as above without the num_comparisons argument
template <typename Ty> int wpartition(Ty* x, Ty* w, Ty piv, int N, int& lrb, Ty& W0, Ty& Wl, Ty& Wr);
template <typename Ty> int  partition(Ty* x,        Ty piv, int N);
double select_floyd_rivest_delta(double *a, int l, int r, int k, int delta);
template <typename Ty> Ty medianf_equidist(Ty* a, Ty* w, int l, int r, int k, int dist);
#endif

// Returns the optimal OS (k) for a given alpha & M
int getK(int M, double alph);
// Returns the optimal M to a given N
int getM(int N);
// Function used by getK to return f(alph,M)/f'(alph,M)
double fdivfp(double xold, double alph, int M);
// Interpolation function used to get the optimal beta to a given alpha (used for 3+ iteration)
double approx_fct(double a);
// The expected cost at x & alpha & M
double cost_fct(double x, double alph, int M);

// A simple ShellSort routine working on x
template <typename Ty> void ShellSort(Ty* x, unsigned size);
// A simple ShellSort routine working on x
// Will also work with NaN (Not A Number)
template <typename Ty> void ShellSortNanAware(Ty* x, unsigned size);
// A simple ShellSort routine working on x and swapping w with each x
template <typename Ty> void ShellSortPair(Ty* x, Ty* w, unsigned size);
// A simple ShellSort routine working on x and swapping w with each x
// Will also work with NaN (Not A Number)
template <typename Ty> void ShellSortPairNanAware(Ty* x, Ty* w, unsigned size);

// Little debugger function to check if Wl,Wr,x,w etc are in agreement
template <typename Ty> void checkConsistency(Ty* x, Ty* w, Ty W0_ref, Ty Wl_ref, Ty Wr_ref, int left, int right, int N);

// Helper function which checks if a given value is NaN (Not A Number)
template<typename T> inline bool isnan(T value)
{
  return value != value;
}

// Helper function which check if a given value is infinite
template<typename T> inline bool isinf(T value)
{
  return std::numeric_limits<T>::has_infinity && (value == std::numeric_limits<T>::infinity() || value == -std::numeric_limits<T>::infinity());
}


/**
 * Returns the median of the vector x[0..N0-1]. First pivot can be specified. If not
 * desired a null pointer needs to be passed.
 * Note: If at any point the pivot happens to be NaN then NaN is returned. However the
 * algorithm will not loop over the array to search for NaNs.
 */
#if _DEBUG
template <typename Ty> Ty medianf(Ty* x, Ty* pivot_in, int N0, int& num_comparisons)
#else
template <typename Ty> Ty medianf(Ty* x, Ty* pivot_in, int N0)
#endif
{
  const int K = 32; // If problem size is smaller than this abort and find median using quickselect

  int M0; // Number of samples to get the first pivot
  int left  = 0; // Always points to the first element of (reduced) problem
  int right = N0-1; // Always points to the last element of (reduced) problem
  int median = (N0-1)/2;
  Ty pivot = x[0];
  Ty xmax = x[0];
  Ty xmin = x[0];

  // Stop recursion
  if( N0 < 256 ) {
#if _DEBUG
    return select_floyd_rivest(x, left, right, median, num_comparisons);
#else
    return select_floyd_rivest(x, left, right, median);
#endif
  }

  // Calculate the number of samples to get the first pivot
  M0  = getM(N0);

  // If the passed pivot is null pointer then compute the first pivot
  if( pivot_in == 0 || isnan(*pivot_in) ) {
    // Need to choose the first pivot
#if _DEBUG
    pivot = medianf_equidist(&x[left], 0, M0-1, (M0-1)/2, N0/M0, num_comparisons);
#else
    pivot = medianf_equidist(&x[left], 0, M0-1, (M0-1)/2, N0/M0);
#endif
    if( isnan(pivot) ) {
      return pivot;
    }
  } else {
    pivot = *pivot_in;
  }

  // partition will return the index of the first element >(=) than the pivot
#if _DEBUG
  int idx0 = partition(x, pivot, N0, num_comparisons);
#else
  int idx0 = partition(x, pivot, N0);
#endif

  if( idx0 > median ) {
    // Pivot was larger than the median
    right = idx0-1;
    xmax = pivot;
  } else if( idx0 < median ) {
    left = idx0;
    xmin = pivot;
  } else {
    return pivot; // OK since pivot was selection type
  }
  // Note: In case the pivot was the WM we just continue until the end where
  // the standard Quickselect solves the remaining problem.
  int N1;

  // While the problem isn't bounded above and below:
  while( (left == 0) || (right == N0-1) )
  {
    double alph;
    N1 = right - left + 1; // Number of samples left

    if( N1 < K ) {
#if _DEBUG
      return quickSelect(&x[left], 1, median-left, N1, num_comparisons);
#else
      return quickSelect(&x[left], 1, median-left, N1);
#endif
    }

    if( left == 0 ) {
      alph = ((double)N0/2) / (double)(N1);
    } else {
      alph = 1 - ((double)N0/2) / (double)(N1);
    }
    if( alph <= 0 || alph >= 1 ) {
      // Found the median
      return pivot;
    }
    int M1 = (int)((double)M0*pow((double)N1/N0,0.618764));

    // Let M1 never be less than 5:
    M1 = std::max(M1,5);
    if( M1 > N1 ) {
      M1 = N1;
    }

    int k = getK(M1, alph);

    // Never choose min/max as pivot
    k = std::max(3,k); // Never select the min
    k = std::min(M1-2,k); // Never select the max
    k = std::min(M1,k);
    k = std::max(1,k);

#if _DEBUG
    pivot = medianf_equidist(&x[left], 0, M1-1, k-1, N1/M1, num_comparisons);
#else
    pivot = medianf_equidist(&x[left], 0, M1-1, k-1, N1/M1);
#endif
    // Partitioning will fail if pivot happend to be NaN
    if( isnan(pivot) ) {
      return pivot;
    }

    // partition will return the index of the first element >(=) to the pivot
#if _DEBUG
    int idx1 = partition(&x[left], pivot, N1, num_comparisons);
#else
    int idx1 = partition(&x[left], pivot, N1);
#endif

    // idx1 can still be 0 even though k >= 3 if there are many eq. elements
    if( idx1 == 0 ) {
      idx1 = 2;
    }
    idx1 += left; // partition() will count from x[left]

    if( idx1 > median ) {
      right = idx1-1;
      xmax = pivot;
    } else if( idx1 < median ) {
      left = idx1;
      xmin = pivot;
    } else {
      return pivot;
    }  
  }

  int N2plus = right - left + 1;
  int elements_removed = 0;

  // Start reducing the set using the approximation of the median
  for(;;)
  {
    if( N2plus <= K )
      break;
    // If any of the prev pivots were inf then we cannot take a linear
    // combination of them and have to fall back.
    if( isinf(xmin) || isinf(xmax) )
      break;

    // 0 < a,b < 1
    double a = (double)(N0/2 - left)  / (double)(right - left);
    double c = approx_fct(a);
    // This pivot is problably not in the set.
    pivot = (Ty)(c*xmax + (1-c)*xmin);

    // partition will return the index of the first element >(=) to the pivot
#if _DEBUG
    int idx2plus = partition(&x[left], pivot, N2plus, num_comparisons);
#else
    int idx2plus = partition(&x[left], pivot, N2plus);
#endif
    idx2plus += left; // partition() will count from x[left]
    if( idx2plus > median ) {
      // Pivot was larger than the median
      right = idx2plus-1;
      xmax = pivot;
    } else {
      // Pivot was smaller than the median
      left = idx2plus;
      xmin = pivot;
    }
    elements_removed = N2plus - (right - left + 1);
    N2plus = right - left + 1;
    if( elements_removed <= 2 ) {
      // We're stuck -> fallback to old method
      break;
    }
  }
  // variable median is one-based index
#if _DEBUG
  return quickSelect(&x[left], 1, median-left, N2plus, num_comparisons);
#else
  return quickSelect(&x[left], 1, median-left, N2plus);
#endif
}

// Finds the median of the elements x[0], x[dist], x[2dist] ... x[(N-1)*dist]
#if _DEBUG
template <typename Ty> Ty medianf_equidist(Ty* a, Ty* w, int l, int r, int k, int dist, int& num_comparisons)
#else
template <typename Ty> Ty medianf_equidist(Ty* a, Ty* w, int l, int r, int k, int dist)
#endif
{
  int n, i, j, s, sd, ll, rr;
  double z, t;

  if( dist <= 0 )
    dist = 1;

  while( r > l )
  {
     if( (r-l) > 600 )
     {
       n = r-l+1;
       i = k-l+1;
       z = log((double)n);
       s = (int)(0.5 * exp(2*z/3));
       sd = (int)(0.5 * sqrt((double)z*s*(n-s)/n));
       if( i-n/2.0 < 0 ) {
         sd = -sd;
       }
       ll = (int)std::max((double)l, k-(double)i*s/n+sd);
       rr = (int)std::min((double)r, k+(double)(n-i)*s/n+sd);
#if _DEBUG
       medianf_equidist(a, w, ll, rr, k, dist, num_comparisons);
#else
       medianf_equidist(a, w, ll, rr, k, dist);
#endif
     }
     t = a[k*dist]; // Pivot
     i = l;
     j = r;
     swap(a[l*dist], a[k*dist]);
     swap(w[l*dist], w[k*dist]);
     if (a[r*dist] > t) {
       swap(a[r*dist], a[l*dist]);
       swap(w[r*dist], w[l*dist]);
     }
     while( i < j ) {
       swap(a[i*dist], a[j*dist]);
       swap(w[i*dist], w[j*dist]);
       i++; j--;
#if _DEBUG
       num_comparisons += 2;
#endif
       while( a[i*dist] < t ) {
         i++;
#if _DEBUG
         num_comparisons++;
#endif
       }
       while( t < a[j*dist] ) {
         j--;
#if _DEBUG
         num_comparisons++;
#endif
       }
     }
     if( a[l*dist] == t ) {
       swap(a[l*dist], a[j*dist]);
       swap(w[l*dist], w[j*dist]);
     } else {
       j++;
       swap(a[j*dist], a[r*dist]);
       swap(w[j*dist], w[r*dist]);
     }
     if( j <= k ) l = j+1;
     if( k <= j ) r = j-1;
#if _DEBUG
     num_comparisons += 2;
#endif
  }
  return a[k*dist];
}

// Finds the median of the elements x[0], x[dist], x[2dist] ... x[(N-1)*dist]
#if _DEBUG
template <typename Ty> Ty medianf_equidist(Ty* a, int l, int r, int k, int dist, int& num_comparisons)
#else
template <typename Ty> Ty medianf_equidist(Ty* a, int l, int r, int k, int dist)
#endif
{
  int n, i, j, s, sd, ll, rr;
  double z, t;

  if( dist <= 0 )
    dist = 1;

  while( r > l )
  {
     if( (r-l) > 600 )
     {
       n = r-l+1;
       i = k-l+1;
       z = log((double)n);
       s = (int)(0.5 * exp(2*z/3));
       sd = (int)(0.5 * sqrt((double)z*s*(n-s)/n));
       if( i-n/2.0 < 0 ) {
         sd = -sd;
       }
       ll = (int)std::max((double)l, k-(double)i*s/n+sd);
       rr = (int)std::min((double)r, k+(double)(n-i)*s/n+sd);
#if _DEBUG
       medianf_equidist(a, ll, rr, k, dist, num_comparisons);
#else
       medianf_equidist(a, ll, rr, k, dist);
#endif
     }
     t = a[k*dist]; // Pivot
     i = l;
     j = r;
     swap(a[l*dist], a[k*dist]);
     if (a[r*dist] > t) {
       swap(a[r*dist], a[l*dist]);
     }
     while( i < j ) {
       swap(a[i*dist], a[j*dist]);
       i++; j--;
#if _DEBUG
       num_comparisons += 2;
#endif
       while( a[i*dist] < t ) {
         i++;
#if _DEBUG
         num_comparisons++;
#endif
       }
       while( t < a[j*dist] ) {
         j--;
#if _DEBUG
         num_comparisons++;
#endif
       }
     }
     if( a[l*dist] == t ) {
       swap(a[l*dist], a[j*dist]);
     } else {
       j++;
       swap(a[j*dist], a[r*dist]);
     }
     if( j <= k ) l = j+1;
     if( k <= j ) r = j-1;
#if _DEBUG
     num_comparisons += 2;
#endif
  }
  return a[k*dist];
}

// Note: pivot is not necessarily in the set
// lrb:
//   - (-1): count from left
//   - (+1): count from right
//   - all other:  count both
// WO is only written if lrb is not +-1.
// lrb is also a return value and is set to either -1 or +1:
//   o -1 indicates that the pivot was less than the WM
//   o +1 indicates that the pivot was less than the WM
// Wl and Wr will be set to the sum of the weigts of the elements removed.
#if _DEBUG
template <typename Ty> int partition(Ty* x, Ty pivot, int N, int& num_comparisons)
#else
template <typename Ty> int partition(Ty* x, Ty pivot, int N)
#endif
{
#if 0
  int left   = 0;
  int right  = N-1;

  while( left <= right ) {
    while( x[left] < pivot ) {
      left++;
      num_comparisons++;
    }
    while( x[right] > pivot ) {
      right--;
      num_comparisons++;
    }
    // +2 since the very last comp. wasnt counted
    num_comparisons += 2;
    if( left <= right ) {
      swap(x[left], x[right]);
      left++;
      right--;
    }
  }
  return std::min(left,N);
#endif
#if 0
  Ty* left  = &x[0];
  Ty* right = &x[N-1];
  while (left <= right) {
    while (*left < pivot) {
      ++left; num_comparisons++;
    }
    while (*right >= pivot) {
      --right; num_comparisons++;
    }
    num_comparisons += 2;
    if (left < right) {
      swap(*left, *right);
      ++left;
      --right;
    }
  }
  return (int)(left - &x[0]);
#endif
#if 0
  int i = 0;
  int j = N-1;
  while( i < j ) {
    while( x[i] < pivot ) { i++; num_comparisons++; }
    while( pivot < x[j] ) { j--; num_comparisons++; }
    swap(x[i], x[j]);
    i++; j--;
    num_comparisons += 2;
  }
  i--; j++;
  if( i > j && x[i] < x[j] ) {
    swap(x[i], x[j]);
    return i;
  }
  return j;
#endif
#if 1
  int low = -1;
  int high = N;
  /* Nibble from each end towards middle, swapping items when stuck */
  for(;;) {
    do {
      low++;
#if _DEBUG
      num_comparisons++;
#endif
    } while( x[low] < pivot );
    do {
      high--;
#if _DEBUG
      num_comparisons++;
#endif
    } while( x[high]  >= pivot );

    if( high <= low )
      break;

    swap(x[low], x[high]);
  }
  return low; // Index of first element which is >(=)
#endif
}



#if _DEBUG
template <typename Ty> Ty wmedianf(Ty* x, Ty* w, Ty* pivot_in, int N0, Ty W0, int& num_comparisons)
#else
template <typename Ty> Ty wmedianf(Ty* x, Ty* w, Ty* pivot_in, int N0, Ty W0)
#endif
{
  const int K = 32; // If problem size is smaller than this abort and find median using quickselect
  int M0; // Number of samples to get the first pivot
  int left  = 0; // Always points to the first element of (reduced) problem
  int right = N0-1; // Always points to the last element of (reduced) problem
  Ty pivot = x[0];
  Ty xmax = x[0];
  Ty xmin = x[0];
  Ty Wl = 0.0; // Sum of weight discarded to the left (0...left-1)
  Ty Wr = 0.0; // Sum of weight discarded to the right (right-1...END)

  // Stop recursion
  if( N0 < 256 ) {
#if _DEBUG
    return wquickSelect(x, w, N0, W0, num_comparisons);
#else
    return wquickSelect(x, w, N0, W0);
#endif
  }

  // Calculate the number of samples to get the first pivot
  M0 = getM(N0);

  // If the passed pivot is null pointer then compute the first pivot
  if( pivot_in == 0 || isnan(*pivot_in) ) {
    // Need to choose the first pivot
#if _DEBUG
    pivot = medianf_equidist(x, w, 0, M0-1, (M0-1)/2, N0/M0, num_comparisons);
#else
    pivot = medianf_equidist(x, w, 0, M0-1, (M0-1)/2, N0/M0);
#endif
    // Partitioning will fail if pivot happend to be NaN
    if( isnan(pivot) ) {
#if _DEBUG
      return wquickSelect(x, w, N0, W0, num_comparisons);
#else
      return wquickSelect(x, w, N0, W0);
#endif
    }
  } else {
    pivot = *pivot_in;
  }

  // partition will return the index of the first element >(=) than the pivot
  int lrb = 0 ;// Indicates to count both and compute W0
#if _DEBUG
  int idx0 = wpartition(x, w, pivot, N0, lrb, W0, Wl, Wr, num_comparisons);
#else
  int idx0 = wpartition(x, w, pivot, N0, lrb, W0, Wl, Wr);
#endif

  if( lrb > 0 ) {
    // Pivot was larger than the median
    right = idx0;
    xmax = pivot;
  } else {
    // Pivot was smaller than the median
    left = idx0;
    xmin = pivot;
  }
  checkConsistency(x, w, W0, Wl, Wr, left, right, N0);

  // Note: In case the pivot was the WM we just continue until the end where
  // the standard Quickselect solves the remaining problem.

  int N1 = right - left + 1; // Number of samples left

  // While the problem isn't bounded above and below:
  while( (left == 0) || (right == N0-1) )
  {
    double alph;
    int M1;

    if( N1 < K ) {
#if _DEBUG
      return wquickSelect(&x[left], &w[left], N1, W0, num_comparisons, Wl, Wr);
#else
      return wquickSelect(&x[left], &w[left], N1, W0, Wl, Wr);
#endif
    }

    // This can occur if negative weights are given.
    if( (Wr < 0) || (Wl < 0) )
      return std::numeric_limits<Ty>::max();

    if( left == 0 ) {
      alph = W0 / (2*W0 - Wr);
      M1 = getScaledM(M0,(W0 - 0.5*Wr)/W0);
    } else {
      alph = 1 - W0 / (2*W0 - Wl);
      M1 = getScaledM(M0,(W0 - 0.5*Wl)/W0);
    }
    if( alph <= 0 || alph >= 1 ) {
      //We found the median
      return pivot;
    }
    // Let M1 never be less than 5:
    M1 = std::max(M1,5);
    if( M1 > N1 )
      M1 = N1;

    int k = getK(M1, alph); 

    // Never choose min/max as pivot
    k = std::max(3,k); // Never select the min
    k = std::min(M1-2,k); // Never select the max

#if _DEBUG
    pivot = medianf_equidist(&x[left], &w[left], 0, M1-1, k-1, N1/M1, num_comparisons);
#else
    pivot = medianf_equidist(&x[left], &w[left], 0, M1-1, k-1, N1/M1);
#endif
    // Partitioning will fail if pivot happend to be NaN
    if( isnan(pivot) )
      return wmSortNanAware(&x[left], &w[left], N1, W0 - Wl);

    // partition will return the index of the first element >(=) to the pivot
    lrb = lrb; // lrb stays the same so to count the hopefully smaller side
    checkConsistency(x, w, W0, Wl, Wr, left, right, N0);
#if _DEBUG
    int idx1 = wpartition(&x[left], &w[left], pivot, N1, lrb, W0, Wl, Wr, num_comparisons);
#else
    int idx1 = wpartition(&x[left], &w[left], pivot, N1, lrb, W0, Wl, Wr);
#endif
    idx1 += left; // partition() will count from x[left]
    if( lrb > 0 ) {
      // Pivot was larger than the median
      right = idx1;
      xmax = pivot;
    } else {
      // Pivot was smaller than the median
      left = idx1;
      xmin = pivot;
    }
    checkConsistency(x, w, W0, Wl, Wr, left, right, N0);
    N1 = right - left + 1; // Update samples left
  }
  int N2plus = N1;
  int elements_removed = 0;

  // Start reducing the set using the approximation of the median
  for(;;)
  {
    if( N2plus <= K )
      break;
    // If any of the prev pivots were inf/NaN then we cannot take a lin
    // combination of them and have to fall back.
    if( isinf(xmin) || isinf(xmax) )
      break;

    // 0 < a < 1
    double a = (W0 - Wl) / (2*W0 - Wl - Wr);
    double c = approx_fct(a);
    // This pivot is problably not in the set.
    pivot = (Ty)(c*xmax + (1-c)*xmin);

    if( c < 0.5 ) {
      lrb = -1;
    } else {
      lrb = +1;
    }
#if _DEBUG
    int idx2plus = wpartition(&x[left], &w[left], pivot, N2plus, lrb, W0, Wl, Wr, num_comparisons);
#else
    int idx2plus = wpartition(&x[left], &w[left], pivot, N2plus, lrb, W0, Wl, Wr);
#endif
    idx2plus += left; // wpartition() will count from x[left]
    if( lrb > 0 ) {
      // Pivot was larger than the median
      right = idx2plus;
      xmax = pivot;
    } else {
      // Pivot was smaller than the median
      left = idx2plus;
      xmin = pivot;
    }
    checkConsistency(x, w, W0, Wl, Wr, left, right, N0);
    elements_removed = N2plus - (right - left + 1);
    N2plus = right - left + 1;
    if( elements_removed <= 2 ) {
      // We're stuck -> fallback to old method
      break;
    }
  }
  // variable median is one-based index
  checkConsistency(x, w, W0, Wl, Wr, left, right, N0);
#if _DEBUG
  return wquickSelect(&x[left], &w[left], N2plus, W0, num_comparisons, Wl, Wr);
#else
  return wquickSelect(&x[left], &w[left], N2plus, W0, Wl, Wr);
#endif
}



// Note: pivot is not necessarily in the set
// lrb:
//   - (-1): count from left
//   - (+1): count from right
//   - all other:  count both
// WO is only written if lrb is not +-1.
// lrb is also a return value and is set to either -1 or +1:
//   o -1 indicates that the pivot was less than the WM
//   o +1 indicates that the pivot was less than the WM
// Wl and Wr will be set to the sum of the weigts of the elements removed.
#if _DEBUG
template <typename Ty> int wpartition(Ty* x, Ty* w, Ty pivot, int N, int& lrb, Ty& W0, Ty& Wl, Ty& Wr, int& num_comparisons)
#else
template <typename Ty> int wpartition(Ty* x, Ty* w, Ty pivot, int N, int& lrb, Ty& W0, Ty& Wl, Ty& Wr)
#endif
{
  // summed up concomitant weights which are smaller/larger than the pivot
  Ty wleftSum  = 0;
  Ty wrightSum = 0;

  int left   = 0;
  int right  = N-1;

  if( lrb == -1 )
  {
    // W0 must be provided if lrb is -1
    if( W0 <= 0.0 )
      throw "No W0 provided!";

    // Sum up the weights left of the pivot
    while( left <= right ) {
      while( x[left] < pivot ) {
        wleftSum += w[left];
        left++;
#if _DEBUG
        if( left > N )
          throw "Ran out of bounds";
        num_comparisons++;
#endif
      }
      while( x[right] > pivot ) {
        right--;
#if _DEBUG
        if( right < -1 )
          throw "Ran out of bounds";
        num_comparisons++;
#endif
      }
#if _DEBUG
      // +2 since the very last comp. wasnt counted
      num_comparisons += 2;
#endif
      if( left <= right ) {
        swap(x[left], x[right]);
        swap(w[left], w[right]);
        wleftSum += w[left];
        left++;
        right--;
      }
    }
    // In case there were elements equal to the pivot, we need to adjust the sum of weights
    while( left > right+1 ) {
      left--;
      wleftSum -= w[left];
    }
    if( wleftSum + Wl < W0 ) {
      Wl += wleftSum;
      lrb = -1;
      return left;
    } else {
      Wr = 2*W0 - wleftSum - Wl;
      lrb = +1;
      return right;
    }
  }
  else if( lrb == 1 )
  {
    if( W0 <= 0.0 )
      throw "No W0 provided!";

    // Sum up the weights right of the pivot
    while( left <= right ) {
      while( x[left] < pivot ) {
        left++;
#if _DEBUG
        if( left > N )
          throw "Ran out of bounds";
        num_comparisons++;
#endif
      }
      while( x[right] > pivot ) {
        wrightSum += w[right];
        right--;
#if _DEBUG
        if( right < -1 )
          throw "Ran out of bounds";
        num_comparisons++;
#endif
      }
#if _DEBUG
      // +2 since the very last comp. wasnt counted
      num_comparisons += 2;
#endif
      if( left <= right ) {
        swap(x[left], x[right]);
        swap(w[left], w[right]);
        wrightSum += w[right];
        ++left;
        --right;
      }
    }
    // In case there were elements equal to the pivot, we need to adjust the sum of weights
    while( left > right+1 ) {
      left--;
      wleftSum -= w[left];
    }
    if( wrightSum + Wr < W0 ) {
      Wr += wrightSum;
      lrb = +1;
      return right;
    } else {
      Wl = 2*W0 - wrightSum - Wr;
      lrb = -1;
      return left;
    }
  }
  else
  {
    // Sum up the weights right of the pivot
    while( left <= right ) {
      while( x[left] < pivot ) {
        wleftSum += w[left];
        left++;
#if _DEBUG
        if( left > N )
          throw "Ran out of bounds";
        num_comparisons++;
#endif
      }
      while( x[right] > pivot ) {
        wrightSum += w[right];
        right--;
#if _DEBUG
        if( right < -1 )
          throw "Ran out of bounds";
        num_comparisons++;
#endif
      }
      // +2 since the very last comp. wasnt counted
#if _DEBUG
      num_comparisons += 2;
#endif
      if( left <= right ) {
        swap(x[left], x[right]);
        swap(w[left], w[right]);
        wleftSum += w[left];
        wrightSum += w[right];
        ++left;
        --right;
      }
    }
    // In case there were elements equal to the pivot, we need to adjust the sum of weights
    while( left > right+1 ) {
      left--;
      wleftSum -= w[left];
    }
    W0 = 0.5*(wleftSum + wrightSum);
    if( wleftSum < W0 ) {
      Wl = wleftSum;
      lrb = -1;
      return left;
    } else {
      Wr = wrightSum;
      lrb = +1;
      return right;
    }
  }
}


/*
 * This Quickselect routine is based on the algorithm described in
 * "Numerical recipes in C", Second Edition,
 * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 * This code by Nicolas Devillard - 1998. Public domain.
 */
#if _DEBUG
template <typename Ty> Ty quickSelect(Ty* arr, int delta, int k, int n, int& num_comparisons)
#else
template <typename Ty> Ty quickSelect(Ty* arr, int delta, int k, int n)
#endif
{
  int low = 0;
  int high = n-1;
  int middle, ll, hh;

  for(;;)
  {
    if( high <= low ) { /* One element only */
      return arr[k*delta];
    }
    if( high == low + 1 ) { /* Two elements only */
      if (arr[low*delta] > arr[high*delta])
        swap(arr[low*delta], arr[high*delta]) ;
      return arr[k*delta] ;
    }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high)/2;
    // Take 3-sample median of low,middle & high and put it in arr[low]
    if( arr[middle*delta] > arr[high*delta])    swap(arr[middle*delta], arr[high*delta]) ;
    if( arr[low*delta] > arr[high*delta])       swap(arr[low*delta], arr[high*delta]) ;
    if( arr[middle*delta] > arr[low*delta])     swap(arr[middle*delta], arr[low*delta]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    swap(arr[middle*delta], arr[(low+1)*delta]);

    if( isnan(arr[low]) )
      return arr[low];

#if _DEBUG
    num_comparisons += 4;
#endif
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do {
        ll++;
#if _DEBUG
        num_comparisons++;
#endif
      } while (arr[low*delta] > arr[ll*delta]);
      do {
        hh--;
#if _DEBUG
        num_comparisons++;
#endif
      } while (arr[hh*delta]  > arr[low*delta]);

#if _DEBUG
      num_comparisons += 2;
#endif

      if (hh < ll)
        break;

      swap(arr[ll*delta], arr[hh*delta]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    swap(arr[low*delta], arr[hh*delta]) ;

    /* Re-set active partition */
    if (hh <= k)
      low = ll;
    if (hh >= k)
      high = hh - 1;
  }
}


/*
 * This Quickselect routine is based on the algorithm described in
 * "Numerical recipes in C", Second Edition,
 * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 * This code by Nicolas Devillard - 1998. Public domain.
 */
#if _DEBUG
template <typename Ty> Ty wquickSelect(Ty* x, Ty* w, int N, Ty W0, int& num_comparisons, Ty Wl, Ty Wr)
#else
template <typename Ty> Ty wquickSelect(Ty* x, Ty* w, int N, Ty W0, Ty Wl, Ty Wr)
#endif
{
  int low = 0;
  int high = N-1;
  int middle, ll, hh;
  Ty wSumLeft = 0.0;
  Ty wLeft  = Wl; // The sum of weights of discarded element to the left
  Ty wRight = Wr; // ... to the right

  if( W0 < 0 ) {
    // W0 is not provided -> compute
    W0 = 0.0;
    for( int i = 0; i < N; i++ )
      W0 += w[i];
    W0 *= 0.5;
  }

  for (;;)
  {
    if( high <= low ) /* One element only */
    {
      return x[low];
    }

    if( high == low + 1 ) { /* Two elements only */
      if (x[low] > x[high]) {
        swap(x[low], x[high]);
        swap(w[low], w[high]);
      }
      if( wLeft + w[low] >= W0 ) {
        return x[low];
      } else {
        return x[high];
      }
    }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high)/2;
    // Take 3-sample median of low,middle & high and put it in x[low]
    if( x[middle] > x[high] ) {
      swap(x[middle], x[high]);
      swap(w[middle], w[high]);
    }
    if( x[low] > x[high] ) {
      swap(x[low], x[high]);
      swap(w[low], w[high]);
    }
    // Now the largest sample is in x[high]
    if( x[middle] > x[low] ) {
      swap(x[middle], x[low]);
      swap(w[middle], w[low]);
    }
    // x[middle] < x[low] < x[high]

    /* Swap low item (now in position middle) into position (low+1) */
    swap(x[middle], x[low+1]);
    swap(w[middle], w[low+1]);

    // We will use x[low] as the pivot
    Ty piv = x[low];
    if( isnan(x[low]) )
      return wmSortNanAware(&x[low], &w[low], high-low+1, W0 - wLeft);

#if _DEBUG
    num_comparisons += 3;
#endif

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do {
        wSumLeft += w[ll];
        ll++;
#if _DEBUG
        num_comparisons++;
#endif
      } while (piv > x[ll]);
      do {
        hh--;
#if _DEBUG
        num_comparisons++;
#endif
      } while (x[hh] > piv);

#if _DEBUG
      num_comparisons += 2;
#endif

      if( hh < ll )
        break;

      // Swap those items and their corresponding weight
      swap(x[ll], x[hh]);
      swap(w[ll], w[hh]);
    }

    /* Swap middle item (in position low) back into correct position */
    swap(x[low], x[hh]);
    swap(w[low], w[hh]);
    //wSumLeft += w[hh]; Do not include the pivot

    // In case there were elements equal to the pivot, we need to adjust the sum of weights
    while( ll > hh+1 ) {
      ll--;
      wSumLeft -= w[ll];
    }

    // Note: wSumLeft is the sum of all weights left of the pivot excluding the
    // pivot itself.
    /* Re-set active partition */
    if( (wLeft + wSumLeft) <= W0 ) {
      if( (wLeft + wSumLeft + w[hh]) > W0 ) {
        // The pivot is the WM -> return
        return x[hh];
      }
      // The pivot was smaller than the WM -> discard elements to the left
      low = hh+1; // Set the left limit right to the pivot
      wLeft += w[hh] + wSumLeft;
    } else {
      // The pivot was larger than the WM -> discard elements to the right
      high = hh-1; // Set the right limit left to the pivot
      wRight = 2*W0 - wSumLeft - wLeft;
    }
    // Reset the counter
    wSumLeft = 0.0;
  }
}

// Returns the optimal M to a given N
// Where N is the size of the array and M is the size of the subset
int getM(int N)
{
  int M  = (int)(pow((double)N,(double)0.618764) - sqrt((double)N) + 5.0);
  M |= 0x01; // Make M odd
  M = std::max(3, M); // Make M at least 3
  return M;
}

// Returns the optimal new M to a given old M and a factor
int getScaledM(int M, double scale)
{
  return (int)((double)M*pow(scale, 0.618764));
}

// Returns the optimal OS (k) for a given alpha & M
int getK(int M_in, double alph_in)
{
  double alph = alph_in;
  double Md = (double)M_in;
  
  // Make alph always 0..1/2
  if( alph > 0.5 )
    alph = 1 - alph_in;

  // Fitting was done with data a in 1/2..1
  double a = 1-alph;
  // This is the function from Eureqa (MSE ~ 500) (fits from M = 50..12000)
  double xold = 0.998202*a*M_in + 11.421*M_in/pow(92.9553*M_in,a) - 1.86877*pow(a*M_in - a*a*M_in,0.545273);
  xold = M_in - xold; // Invert since a is in 0..1/2

  // Ad hoc adjustment for a starting point left of the optimal point
  xold = 0.5*alph*M_in + 0.5*xold;

  if( M_in*alph < NORMAL_POISSON_LIMIT ) {
    // For poisson use lower starting point to ensure convergence
    xold = 0.6*alph*M_in + 0.4*xold;
    // Make starting point never < 1 -> This helps extremly small values of alph
    if( xold < 1 )
      xold = 1;
  }

  // The final value will go here
  double xnew = 0.0;

  if( M_in < 12000 ) {
    // Fitting holds and Newton method will converge
    xnew = xold - fdivfp(xold, alph, M_in);
    xold = xnew;
    xnew = xold - fdivfp(xold, alph, M_in);
  } else {
    // M > 12000 -> Use bisection method
    double xleft  = alph*M_in;
    double xright = BISECTION_DELTA*xleft;
    double yleft  = cost_fct(xleft, alph, M_in);
    double yright = cost_fct(xright, alph, M_in);
    // First make sure that yleft & yright have different signs
    while( yright <= 0 ) {
      xleft = xright; yleft = yright; // Make the left point the right point
      xright = BISECTION_DELTA*xleft;
      yright = cost_fct(xright, alph, M_in);
    }
    // Now we ensured that yright is positive.
    if( xleft <= 0 ) {
      // This should never happen
      throw "xleft <= 0";
    }
    
    int j = 0;
    while( (xright - xleft)/M_in > BISECTION_TOL ) {
      // Take the midpoint
      double xmid = (xleft + xright)/2;
      double ymid = cost_fct(xmid, alph, M_in);
      if( ymid <= 0 ) {
        // It's on the left side of the cost function
        xleft = xmid; yleft = ymid;
      } else {
        xright = xmid; yright = ymid;
      }
      if( ++j > 8 ) {
        break; // Do at most 8 iterations
      }
    }
    // The final output is going to be the midpoint
    xnew = (xleft + xright)/2;
  }

  xnew += 1.0;

  // Check if first pivot was left/right of median.
  if( alph_in > 0.5 )
    xnew = M_in - xnew + 1.0;

  return (int)xnew;
}

double cost_fct(double x, double alph, int M_in)
{
  if( M_in*alph >= NORMAL_POISSON_LIMIT ) {
    // Use normal approx
    return (2*x-M_in+1)*normpdf(x, alph*M_in, alph*(1-alph)*M_in) + 1 - 2*betai(x+1, M_in-x, alph);
  } else {
    // Use poisson approx
    return (2*x-M_in+1)*poisspdf(x,alph*M_in) + 1;
  }
}


double fdivfp(double x, double alph, int M_in)
{
  double mu = M_in*alph;
  double num, den;

  if( mu >= NORMAL_POISSON_LIMIT ) {
    // Use normal approx
    //fdivfpN = @(k) (2*k-M+1 + (1-2*betainc(alph(Nidx),k+1,M-k))./normpdf(k,alph(Nidx)*M,sig(Nidx))) ...
    //./ (4+(2*k-M+1).*(alph(Nidx)*M-k)./sig(Nidx).^2);
    double sigma = sqrt(mu*(1-alph));
    num = 2*x-M_in+1 + (1-2*betai(x+1, M_in-x, alph)) / normpdf(x, mu, sigma);
    den = 4 + (2*x-M_in+1)*(mu-x)/(sigma*sigma);
  } else {
    // Use poisson approx
    //fdivfpP = @(k) (k - (M-1)/2 + 1/2./poisspdf(k,M*alph(Pidx))) ...
    //  ./ (1 + (k - (M-1)/2).*(log(M*alph(Pidx)./k) - 1./(2*k) + 1/12./k.^2));
    double mu = M_in*alph;
    double xn = x - 0.5*(M_in-1);
    num = xn + 0.5/poisspdf(x,mu);
    den = 1 + xn*log(mu/x) - 1/(2*x);
  }
  return num/den;
}


double approx_fct_1(double x)
{
  return x - x*x + x*x*x + 1.4*x*x*x*x;
}

/**
 * Interpolation function used for calculating the pivot when the method of
 * linear combination of max & min is used.
 */
double approx_fct(double a)
{
  if( a < .5 ) {
    return .5 - approx_fct_1( .5 - a );
  } else {
    return .5 + approx_fct_1(-.5 + a );
  }
}


template <typename Ty> static void inline swap(Ty &a, Ty &b)
{
  Ty temp = a;
  a = b;
  b = temp;
}

/**
 * Solving the WM problem by sorting. If W0 is unknown it can be passed a negative
 * value and the function will compute it.
 */
template <typename Ty> Ty wmSort(Ty* x, Ty* w, int N, Ty W0)
{
  ShellSortPair(x,w,N);

  if( W0 < 0.0 ) {
    W0 = 0.0;
    for( int i = 0; i < N; i++ ) {
      W0 += w[i];
    }
    W0 *= 0.5;
  }

  Ty wSum = 0.0;
  int i;

  for( i = 0; i < N; i++ )
  {
    wSum += w[i];
    if( wSum >= W0 ) {
      break;
    }
  }
  return x[i];
}

/**
 * Solving the WM problem by sorting. If W0 is unknown it can be passed a negative
 * value and the function will compute it.
 * This implementation can handle NaNs (Not A Number) in x.
 */
template <typename Ty> Ty wmSortNanAware(Ty* x, Ty* w, int N, Ty W0)
{
  ShellSortPairNanAware(x,w,N);

  if( W0 < 0.0 ) {
    W0 = 0.0;
    for( int i = 0; i < N; i++ ) {
      W0 += w[i];
    }
    W0 *= 0.5;
  }

  Ty wSum = 0.0;
  int i;

  for( i = 0; i < N; i++ )
  {
    wSum += w[i];
    if( wSum >= W0 ) {
      break;
    }
  }
  return x[i];
}

/**
 * A shell sorting routine which works on x and shuffles the pair (x,w) accordingly
 */
template<typename Ty>
void ShellSortPair(Ty* x, Ty* w, unsigned size)
{
    const unsigned hmax = size/9;
    unsigned h;
    for(h = 1; h <= hmax; h = 3*h+1);
    for(; h > 0; h /= 3)
    {
        for(unsigned i = h; i < size; ++i)
        {
            Ty v  = x[i];
            Ty vw = w[i];
            unsigned j = i;
            while(j >= h && v < x[j-h])
            {
                x[j] = x[j-h];
                w[j] = w[j-h];
                j -= h;
            }
            x[j] = v;
            w[j] = vw;
        }
    }
}

/**
 * A shell sorting routine which works on x and shuffles the pair (x,w) accordingly
 * This is the NaN (Not A Number) aware implementation
 */
template<typename Ty>
void ShellSortPairNanAware(Ty* x, Ty* w, unsigned size)
{
    const unsigned hmax = size/9;
    unsigned h;
    for(h = 1; h <= hmax; h = 3*h+1);
    for(; h > 0; h /= 3)
    {
        for(unsigned i = h; i < size; ++i)
        {
            Ty v  = x[i];
            Ty vw = w[i];
            if( isnan(v) )
            {
              unsigned j = i;
              while(j >= h)
              {
                  x[j] = x[j-h];
                  w[j] = w[j-h];
                  j -= h;
              }
              x[j] = v;
              w[j] = vw;
            } else {
              unsigned j = i;
              while(j >= h && v < x[j-h])
              {
                  x[j] = x[j-h];
                  w[j] = w[j-h];
                  j -= h;
              }
              x[j] = v;
              w[j] = vw;
            }
        }
    }
}

/**
 * A shell sorting routine which works on x
 */
template<typename ItemType>
void ShellSort(ItemType* x, unsigned size)
{
    const unsigned hmax = size/9;
    unsigned h;
    for(h = 1; h <= hmax; h = 3*h+1);
    for(; h > 0; h /= 3)
    {
        for(unsigned i = h; i < size; ++i)
        {
            ItemType v = x[i];
            unsigned j = i;
            while(j >= h && v < x[j-h])
            {
                x[j] = x[j-h];
                j -= h;
            }
            x[j] = v;
        }
    }
}


/**
 * Floyd and Rivest SELECT algorithm.
 */
#if _DEBUG
double select_floyd_rivest(double *a, int l, int r, int k, int& num_comparisons)
#else
double select_floyd_rivest(double *a, int l, int r, int k)
#endif
{
  int n, i, j, s, sd, ll, rr;
  double z, t;

  while( r > l )
  {
     if( (r-l) > 600 )
     {
       n = r-l+1;
       i = k-l+1;
       z = log((double)n);
       s = (int)(0.5 * exp(2*z/3));
       sd = (int)(0.5 * sqrt((double)z*s*(n-s)/n));
       if( i-n/2.0 < 0 ) {
         sd = -sd;
       }
       ll = (int)std::max((double)l, k-(double)i*s/n+sd);
       rr = (int)std::min((double)r, k+(double)(n-i)*s/n+sd);
#if _DEBUG
       select_floyd_rivest(a, ll, rr, k, num_comparisons);
#else
       select_floyd_rivest(a, ll, rr, k);
#endif
     }
     t = a[k];
     i = l;
     j = r;
     swap(a[l], a[k]);
     if (a[r] > t) {
       swap(a[r], a[l]);
     }
     while (i<j) {
       swap(a[i], a[j]);
       i++; j--;
#if _DEBUG
       num_comparisons += 2;
#endif
       while( a[i] < t ) {
         i++;
#if _DEBUG
         num_comparisons++;
#endif
       }
       while( t < a[j] ) {
         j--;
#if _DEBUG
         num_comparisons++;
#endif
       }
     }
     if( a[l] == t ) {
       swap(a[l], a[j]);
     } else {
       j++;
       swap(a[j], a[r]);
     }
     if( j <= k ) l = j+1;
     if( k <= j ) r = j-1;
#if _DEBUG
     num_comparisons += 2;
#endif
  }
  return a[k];
}

/**
 * Floyd and Rivest SELECT algorithm.
 * Takes an additional parameter delta which is an integer specifying the inter sample
 * distance
 */
#if _DEBUG
double select_floyd_rivest_delta(double *a, int l, int r, int k, int delta, int& num_comparisons)
#else
double select_floyd_rivest_delta(double *a, int l, int r, int k, int delta)
#endif
{
  int n, i, j, s, sd, ll, rr;
  double z, t;

  if( ((k-l) % delta != 0) || ((r-l)%delta != 0) )
    return -1.0; // assume k is l + some multiple of delta

  while( r > l )
  {
     if( (r-l) > 600*delta )
     {
       n = (r-l+1)/delta;
       i = (k-l+1)/delta;
       z = log((double)n);
       s = (int)(0.5 * exp(2*z/3));
       sd = (int)(0.5 * sqrt((double)z*s*(n-s)/n));
       if( i-n/2.0 < 0 ) {
         sd = -sd;
       }
       ll = std::max(l, k - ((int)((double)i*s/n+sd))*delta);
       rr = std::min(r, k + ((int)((double)(n-i)*s/n+sd))*delta);
#if _DEBUG
        select_floyd_rivest(a, ll, rr, k, num_comparisons);
#else
        select_floyd_rivest(a, ll, rr, k);
#endif
     }
     t = a[k];
     i = l;
     j = r;
     swap(a[l], a[k]);
     if (a[r] > t) {
       swap(a[r], a[l]);
     }
     while (i<j) {
       swap(a[i], a[j]);
       i += delta; j -= delta;
#if _DEBUG
       num_comparisons += 2;
#endif
       while( a[i] < t ) {
         i += delta;
#if _DEBUG
         num_comparisons++;
#endif
       }
       while( t < a[j] ) {
         j -= delta;
#if _DEBUG
         num_comparisons++;
#endif
       }
     }
     if( a[l] == t ) {
       swap(a[l], a[j]);
     } else {
       j += delta;
       swap(a[j], a[r]);
     }
     if( j <= k ) l = j+delta;
     if( k <= j ) r = j-delta;
#if _DEBUG
     num_comparisons += 2;
#endif
  }
  return a[k];
}

/**
 * A debug routine to check if the values Wl, W0, Wr etc. still are in agreement
 * with each other
 */
template <typename Ty> void checkConsistency(Ty* x, Ty* w, Ty W0_ref, Ty Wl_ref, Ty Wr_ref, int left, int right, int N)
{
#if _DEBUG
  Ty Wl = 0.0;
  Ty Wr = 0.0;
  Ty W0 = 0.0;

  int i = 0;
  for( ; i < left; i++ ) {
    Wl += w[i];
  }
  for( ; i <= right; i++ ) {
    W0 += w[i];
  }
  for( ; i < N; i++ ) {
    Wr += w[i];
  }
  W0 += Wr + Wl;
  W0 *= 0.5;

  if( abs(Wl - Wl_ref) > 1e-3 )
    throw "Wl not correct.\n";
  if( abs(Wr - Wr_ref) > 1e-3 )
    throw "Wr not correct.\n";
  if( abs(W0 - W0_ref) > 1e-3 )
    throw "W0 not correct.\n";
#endif
}



// Instantiate templates
#if _DEBUG
template double wmedianf<double>(double* x, double* w, double* pivot, int N, double W0, int& num_comparisons);
template double medianf<double>(double* x, double* pivot, int N, int& num_comparisons);
template double quickSelect<double>(double* x, int delta, int k, int n, int& num_comparisons);
template double wquickSelect(double* x, double* w, int N, double W0, int& num_comparisons, double Wl, double Wr);
#else
template double wmedianf<double>(double* x, double* w, double* pivot, int N, double W0);
template double medianf<double>(double* x, double* pivot, int N);
template double quickSelect<double>(double* x, int delta, int k, int n);
template double wquickSelect(double* x, double* w, int N, double W0, double Wl, double Wr);
#endif
template double wmSort(double* x, double* w, int N, double W0);
