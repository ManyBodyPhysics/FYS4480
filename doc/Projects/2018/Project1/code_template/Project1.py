import numpy as np
import sys
from project1helper import *

Z = 4				 #Atomic number

h  = OBME(Z) 		 #One-body matrix elements
v_spatial     = spatial_TBME(Z)  #Spatial integral (pq|v|rs) given in project text

n_spatial = v_spatial.shape[0] #Number of spatial basis functions, i.e 1s,2s,3s orbitals
L         = 2*n_spatial        #Total number of basis function with spin included

vAS = antiSymmetrize(v_spatial) #vAS = <pq|v|rs>-<pq|v|sr> with proper spin antisymmetry and has shape (L,L,L,L)

part = slice(0,Z) #particle indices
hole = slice(Z,L) #hole indices

#Use the einstein summation convention to evaluate the reference energy
Eref  = np.einsum('ii->',h[part,part])+0.5*np.einsum('ijij->',vAS[part,part,part,part])
print("Eref: %g" % Eref)