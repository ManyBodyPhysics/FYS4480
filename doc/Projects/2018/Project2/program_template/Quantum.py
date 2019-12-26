import numpy as np
import math

class Quantum:
	def __init__(self,n_fermi,TBME):
		"""
		n_fermi - Number of electrons (number of states below fermi level)
		TBME - function which returns the two-body matrix elements of the Hamiltonian.
		"""
		self.n_fermi = n_fermi
		self.TBME = TBME
	def set_Z(self,Z):
		"""
		Z - Atomic number / coloumb interaction
		"""
		self.Z = Z
		self.u = self.TBME(Z)
	def set_states(self,states):
		"""
		states - list containing lists [n,m_s] (basis states), where n refers to quantum number n of the state
				 and m_s is the secondary spin number of the state.
		"""
		self.n_basis = len(states)
		self.states = states
	def ihj(self,i,j):
		if self.states[i,1] == self.states[j,1] and self.states[i,0] == self.states[j,0]:
			return(-self.Z**2/(2*self.states[i,0]**2))
		else:
			return(0)
	def ijvkl(self,i,j,k,l):
		if self.states[i,1] == self.states[k,1] and self.states[j,1] == self.states[l,1]:
			if self.states[i,0] == self.states[j,0] and self.states[i,1] == self.states[j,1]:
				return(0)
			elif self.states[k,0] == self.states[l,0] and self.states[k,1] == self.states[l,1]:
				return(0)
			else:
				return(self.u[self.states[i,0]-1,self.states[j,0]-1,\
						self.states[k,0]-1,self.states[l,0]-1])
		else:
			return(0)
	

	def ijvklAS(self,i,j,k,l):
		return(self.ijvkl(i,j,k,l) - self.ijvkl(i,j,l,k))

	


		


						




