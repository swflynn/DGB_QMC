#=============================================================================80
#               Distributed Gaussian Basis (Eigen Error Analysis)
#==============================================================================#
#       Discussion:
#Python 2 script for computing simple error analysis for eigenvalue calculations 
#Error is simply the absolute difference between computed and exact result
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
#errors             ==> absolute difference between computed/theory eigenvalues
#index              ==> eigenvalue index (for plotting convenience)
#df_analysis        ==> store all data to write to file (for convenience)
#==============================================================================!
#                    Read In Eigenvalues and Theory Values
#==============================================================================!
script,eigenvalues,theory_values=argv
df_eigen=pd.read_csv(eigenvalues,delim_whitespace=True,header=None,dtype=np.float64)
df_theory=pd.read_csv(theory_values,delim_whitespace=True,header=None,dtype=np.float64)
#print 'Input Eigenvalues'
#print df_eigen.head()
#print 'theory input'
#print df_theory.head()
#==============================================================================!
#           Determine the Number of Iterations and Eigenvalues Computed 
#==============================================================================!
num_iter=int(df_eigen.shape[0])  #number of rows
num_eig_val=int(df_eigen.shape[1])     #number of columns
#print'Number of Iterations => ', num_iter
#print'Number of Eigenvalues => ', num_eig_val
#==============================================================================!
#                       Sort Theory Eigenvalues
#LLAPACK Eigensolver sorts eigenvalues smallest to largest
#==============================================================================!
df_sort=df_theory.sort_values([0])
df_sort=df_sort.reset_index(drop=True)
size_theory = int(df_sort.shape[0])
#print('Sorted Theory Eigenvalues')
#print df_sort.head()
#print df_sort.tail()
#==============================================================================!
#                   Remove extra sorted theory values
#Fortran code computes all combinatorics for i,j index
#Only compute i eigenvalues (true=i^2, but eigenvalues=i)
#Remove tail values (theory - number of eigenvalues computed)
#==============================================================================!
df_sort.drop(df_sort.tail((size_theory-num_eig_val)).index,inplace=True)
#print'sorted and truncated data frame'
#print df_sort.head()
#==============================================================================!
#                    Compute Error ===> abs(computed-theory)
#want to compute error for each eigenvalue and each iteration
#==============================================================================!
errors=[]
index=[]
for i in range(0,num_iter):
    index.append(i+1)
    for j in range(0,num_eig_val):
        errors.append(abs(df_eigen[j][i] - df_sort[0][j]))
#==============================================================================!
#                       Write Data to CSV File
#==============================================================================!
df_index=pd.DataFrame(np.array(index))
df_error=pd.DataFrame(np.array(errors).reshape(num_iter,num_eig_val))
df_eigen.to_csv('eigens.dat',sep='\t',index=False)
df_sort.to_csv('sort_theory.dat',sep='\t',index=False)
df_error.to_csv('errors.dat',sep='\t',index=False)
df_analysis=pd.concat([df_index, df_eigen, df_error],axis=1)
df_analysis.to_csv('all.dat',sep='\t',index=False)
