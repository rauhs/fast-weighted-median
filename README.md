# Fast-weighted-median

This repository holds the source code for a C/C++ implementation for my
publication:

A. Rauh, G.R. Arce: **Optimal Pivot Selection in Fast Weighted Median Search**.
IEEE Transactions on Signal Processing (Vol. 60, Issue 8)

# Abstract
Weighted median filters are increasingly being used in signal processing
applications and thus fast implementations are of importance. This paper
introduces a fast algorithm to compute the weighted median (WM) of N samples
which has linear time and space complexity as opposed to O(N log N) which is
the time complexity of traditional sorting algorithms. A popular selection
algorithm often used to find the WM in large data sets is Quickselect whose
performance is highly dependent on how the pivots are chosen. We introduce an
optimization based pivot selection strategy which results in significantly
improved performance as well as a more consistent runtime compared to
traditional approaches. The selected pivots are order statistics of subsets. In
order to find the optimal order statistics as well as the optimal subset sizes,
a set of cost functions are derived, which when minimized lead to optimal
design parameters. We compare the complexity to Floyd and Rivest's algorithm
SELECT which to date has been the fastest median finding algorithm and we show
that the proposed algorithm compared with SELECT requires close to 30% fewer
comparisons. It is also shown that the proposed selection algorithm is
asymptotically optimal for large N.

# IEEE explore link

The journal article is available at:

http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6193457

# License
Copyright © 2014-2015 Andre Rauh.

Distributed under MIT License.

