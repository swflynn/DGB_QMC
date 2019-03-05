#!/usr/bin/env python
#=============================================================================80
#                       Gaussian Distributed Bins
#=============================================================================80
#       Discussion: 
# Python 2 implementation for computing 1D gaussian bins
#Computes the area under a Gaussian, and then generates bins such that the 
#area of each bin is constant.
#Code assumes symmetric gaussian grid centered at 0
#==============================================================================#
#       Modified:
# 13 September 2018
#       Author:
# Shane Flynn
#==============================================================================#
import numpy as np
from sys import argv
from scipy.special import erfinv
from scipy.special import erf
#==============================================================================#
# argv = script, lower_bound, upper_bound, Number_Bins
#==============================================================================#
script,a,b,Npoints = argv
a = float(a)
b = float(b)
Npoints = int(Npoints)
#==============================================================================#
#Compute the area under the gaussain
#Area := int dx e^{-x^2/2} = sqrt(pi/2)erf(x/sqrt(2))
#==============================================================================#
area = np.sqrt(np.pi/2)*(erf(b/np.sqrt(2.))-erf(a/np.sqrt(2.)))
#print('Area under Gaussian [a,b] ==> %r')% area
#==============================================================================#
# Compute constant area for each bin to have
#==============================================================================#
bins = area / Npoints
#print('Area/ Bin ==> %r')% bins
#==============================================================================#
x0 = 0
x = []
x.append(x0)
#==============================================================================#
for i in range(1,(Npoints+1)/2):
    xi = np.sqrt(2)* erfinv( bins / (np.sqrt(np.pi/2)) + erf(x[i-1] / np.sqrt(2))  )
    x.append(xi)
#==============================================================================#
# Gaussian is symmetic, if x0 := 0 then simply collect negative points
#==============================================================================#
y = []
y.append(-x0)
for i in range(1,(Npoints+1)/2):
    yi = -x[i]
    y.append(yi)
#==============================================================================#
# Remove initial point from y to avoid redundancy (x[0] = y[0] := 0)
#==============================================================================#
y.remove(y[0])
#==============================================================================#
# merge lists for convenience
#==============================================================================#
data = x + y
#==============================================================================#
#                       Generate Grid as Permutations
#==============================================================================#
k=0
for i in range(0,(Npoints)):
    for j in range(0,(Npoints)):
        k +=1
        print('%r   %r')%(data[i],data[j])
