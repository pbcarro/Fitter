'''
Python numpy program to solve the rigid rotor asymmetric top eigenvalue problem

Todo: 
	-Hard code a set of test functions for the analyiticly solvable levels

'''

import numpy as np
from numpy import linalg as LA

def F(kappa):	
	return 0.5*(kappa-1.0)	#F function Table 7.3 G&C, solved in Ir representation

def H(kappa):
	return -0.5*(kappa+1.0)	#H function Table 7.3 G&C, solved in Ir representation

def f1(J,K):
	return 0.25*(J*(J+1.0)-K*(K+1.0))*(J*(J+1.0)-(K+1.0)*(K+2.0))	#Equation 7.13 G&C

def EKK(J,K,kappa):
	return F(kappa)*(J*(J+1.0)-K**2.0)+K**2.0	#Equation 7.11 G&C

def EKK2(J,K,kappa):
	return H(kappa)*np.sqrt(f1(J,K))			#Equation 7.12 G&C




delta = 1.0E-3									#The resolution/step size of kappa used to solve the eigenvalue equation
Delta = np.arange(-1.0,1.0+delta,delta)			#Build an array of the kappa values from -1 to 1
JMax = 25										#Max j value to be solved
JRange = range(JMax+1)							#Array of the J values, the +1 is to account for the 0 indexing so this solves to JMax and not JMax-1
for TargetJ	in JRange:							#Iterate through values of J
	MatSize = 2*TargetJ+1						#Each J has 2J+1 eigenvalues
	CurrentMat = np.zeros((MatSize,MatSize))	#Make a matrix of 0s to be filled with our Hamiltonian
	E_tau = np.zeros((len(Delta),MatSize))		#Make a matrix to store the eigenvalues
	print ("Building Matrix %d" % TargetJ)		
	for Index,Kap in enumerate(Delta):			#Iterate through values of Kappa
		for i in range(MatSize):				#Iterate through the first axis of the Hamiltonian
			for j in range(MatSize):			#Iterate through the 2nd axis of the Hamiltonian
				if (i == j):
					CurrentMat[i][j] = EKK(TargetJ,i-TargetJ,Kap)	#Diagonal elements of the Hamiltonian
				if (j == i+2):
					CurrentMat[i][j] = EKK2(TargetJ,i-TargetJ,Kap)	#Off diagonal elements on the first axis
				if (i == j+2):
					CurrentMat[i][j] = EKK2(TargetJ,j-TargetJ,Kap)	#Off diagonal elements on the second axis
		E_tau[Index] = np.sort(LA.eigvals(CurrentMat))	#Compute the eigenvalues for this value of Kappa then sort them 
	E_tau = np.swapaxes(E_tau,0,1)	#Once the matrix is computed flip the axes for my sanity in saving it
	#Save the matrix to file
	FileStr = ""
	for i in range(MatSize):
		for j in range(len(Delta)):
			FileStr += "%.16f\t" % E_tau[i][j]	#Currently saving way more precision than is needed, but it also costs basically nothing
		FileStr += "\n"
	FileHandle = open ("EigenVals/%d_dk_3.dat" % (TargetJ),'w')	#Change file name as needed
	FileHandle.write(FileStr)
	FileHandle.close()
