# Fix the script to work for eigenvalues as a function of iteration.
# 1-8-18
#=============================================================================80
#               Distributed Gaussian Basis (Eigen Error Analysis)
#==============================================================================#
#       Discussion:
#Python 2 script for computing simple error analysis for eigenvalue calculations 
#'Theory' values are computed in the Fortran code 
#Error is taken as the deviation from theory value
#==============================================================================#
#       Modified:
#   19 November 2018
#       Author:
#   Shane Flynn 
#==============================================================================#
import pandas as pd
import numpy as np
from sys import argv
#==============================================================================#
#               Discussion:
#script             ==> Name of this script (error.py)
#eigenvalues        ==> datafile containing computed eigenvalues
#theory_values      ==> datafile containing theoretical eigenvalues
#df_eigen           ==> python dataframe for eigenvalues
#df_theory          ==> python dataframe for theory_values
#df_sort            ==> python dataframe with sorted theory_values
#errors             ==> difference between computed/theory eigenvalues
#index              ==> eigenvalue index (for plotting convenience)
#df_analysis        ==> store all data to write to file (for convenience)
#==============================================================================!
#                    Read In Eigenvalues and Theory Values
#==============================================================================!
script,eigenvalues,theory_values=argv
df_eigen=pd.read_csv(eigenvalues,delim_whitespace=True,header=None,dtype=np.float64)
df_theory=pd.read_csv(theory_values, delim_whitespace=True,header=None,dtype=np.float64)
print 'eigen input'
print(df_eigen)
print 'theory input'
print(df_theory)
# From the eigen matrix, determine how many iterations were computed
num_iter = int(df_eigen.shape[0])  #number of rows
# From the eigen matrix, determine how many eigenvalues were computed
num_eig_val = int(df_eigen.shape[1])     #number of columns
print'iterations'
print num_iter
print num_eig_val
print 'eigenvalues'
#==============================================================================!
#                         Sort Theory Eigenvalues
#LLAPACK Eigensolver sorts eigenvalues smallest to largest
#==============================================================================!
df_sort=df_theory.sort_values([0])
df_sort=df_sort.reset_index(drop=True)
size_theory = int(df_sort.shape[0])
#==============================================================================!
#                   Remove extra sorted theory values
#Fortran code computes all combinatorics for i,j index
#Only compute i eigenvalues (true=i^2, but eigenvalues=i)
#So remove tail (theory - number of eigenvalues computed)
#==============================================================================!
df_sort.drop(df_sort.tail((size_theory-num_eig_val)).index,inplace=True)
print'sorted and truncated data frame'
print df_sort
#==============================================================================!
#                    Compute Error ===> abs(computed-theory)
#want to compute error for each eigenvalue and each iteration
# compute error (loop over iterations, always compare to same theory value)
#==============================================================================!
errors = []
for i in range(0,num_iter):
    for j in range(0,num_eig_val):
        errors.append(abs(df_eigen[j][i] - df_sort[0][j]))
        print'error test'
        print 'i = ', i
        print 'j = ', j
        print 'df_eigen[j,i]', df_eigen[j][i]
        print 'df_sort[0,j]', df_sort[0][j]
        print abs(df_eigen[j][i] - df_sort[0][j])
#==============================================================================!
#               Reshape the error to be consistent with input data
#==============================================================================!
df_error=pd.DataFrame(np.array(errors).reshape(num_iter,num_eig_val))
print 'error boyyy'
print df_error
#==============================================================================!
#                       Write Data to CSV File
# I am here, need to make a plot with xmgrace and figure out the most convenient
#way to plot and store the data
#==============================================================================!
#df_analysis=pd.concat([df_eigen, df_sort, df_error],axis=1)
#df_analysis=pd.concat([df_index, df_eigen, df_sort, df_error],axis=1)
#df_analysis.to_csv('error.dat',sep='\t',index=False)
