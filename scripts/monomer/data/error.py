#=============================================================================80
#               Distributed Gaussian Basis (Eigen Error Analysis)
#==============================================================================#
#       Discussion:
#Python 2 script for computing simple error analysis for eigenvalue calculations 
#'Theory' values are computed in the Fortran code 
#Error is simply the deviation from the theory value
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
#==============================================================================!
#                         Sort Theory Eigenvalues
#LLAPACK Eigensolver sorts eigenvalues smallest to largest
#==============================================================================!
df_sort=df_theory.sort_values([0])
df_sort=df_sort.reset_index(drop=True)
#==============================================================================!
#                   Remove extra sorted theory values
#Fortran code computes all combinatorics for i,j index
#Only compute i eigenvalues (true=i^2, but eigenvalues=i)
#==============================================================================!
df_sort.drop(df_sort.tail((df_sort.shape[0]-df_eigen.shape[0])).index,inplace=True)
#==============================================================================!
#                    Compute Error ===> abs(computed-theory)
#==============================================================================!
errors = []
index = []
for i in range(0,df_eigen.shape[0]):
    errors.append(abs(df_eigen[0][i] - df_sort[0][i]))
    index.append(i+1)
#==============================================================================!
#                       Write Data to CSV File
#==============================================================================!
df_index=pd.DataFrame(np.array(index).reshape(df_eigen.shape[0]))
df_error=pd.DataFrame(np.array(errors).reshape(df_eigen.shape[0]))
df_analysis=pd.concat([df_index, df_eigen, df_sort, df_error],axis=1)
df_analysis.to_csv('error.dat',sep='\t',index=False)
