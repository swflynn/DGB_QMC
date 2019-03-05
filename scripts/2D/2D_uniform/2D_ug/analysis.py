# Python script to compute error for 2D seperable Eigenvalues

import pandas as pd
import numpy as np
from sys import argv

# eigenvalues and true_values are generated in the fortran code
script, eigenvalues, true_values = argv
df_eigen = pd.read_csv(eigenvalues, delim_whitespace=True,header=None,dtype=np.float64)
df_true = pd.read_csv(true_values, delim_whitespace=True,header=None,dtype=np.float64)
#  Need to sort true values (they are not in chronological order, eigen energies are)
df_sort=df_true.sort_values([0])
df_sort=df_sort.reset_index(drop=True)
# Only keep same number of true values as eigenvalues
# The fortran code computes all combinatorics for i,j Index, but only computes
# i eigenvalues, therefore true = i^2, but eigenvlaues = i
df_sort.drop(df_sort.tail((df_sort.shape[0]-df_eigen.shape[0])).index,inplace=True)
# Loop over eigenvalues and compute the abs difference between true/compute
errors = []
for i in range(0,df_eigen.shape[0]):
    errors.append(abs(df_eigen[0][i] - df_sort[0][i]))
# Generate index for eigenvalue being computed, for plotting convenience
index = []
for i in range(0,df_eigen.shape[0]):
    index.append(i+1)
# make the errors/index into dataframes and write them to a file
df_index = pd.DataFrame(np.array(index).reshape(df_eigen.shape[0]))
df_error = pd.DataFrame(np.array(errors).reshape(df_eigen.shape[0]))
df_analysis = pd.concat([df_index, df_eigen, df_sort, df_error],axis=1)
df_analysis.to_csv('all.dat',sep='\t',index=False)
