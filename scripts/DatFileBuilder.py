#!/usr/bin/env python

'''
Program for merging individual eigenvalue files into a full file
-Currently it's own program so I can build catalogs only as large as I want to
'''

import os
import json
import numpy as np

if os.path.exists("EigenVals") is False:
	os.mkdir("EigenVals")

with open("parameters.json") as read_file:
	params = json.load(read_file)

delta = params["delta"]                         #The resolution/step size of kappa used to solve the eigenvalue equation
JStop = params["JMax"]                           #Max j value to be solved			
JStart = 0		#Starting J
Delta = np.arange(-1.0,1.0+delta,delta)	
StateCount = 0
for i in range(JStart,JStop+1):	#Count up how many states there are in total
	StateCount += i*2+1

DatTable = np.zeros((StateCount,len(Delta)))	#Build a matrix to hold all the data
Count = 0
total = 0

for i in range(JStart,JStop+1):		#As with the other programs this is to get the last J value instead of JStop-1
	CurrentState = np.loadtxt("EigenVals/%d_dk_3.dat" % i)
	if (i == 0):		#The J=0 eigenvalue table is a 1D array, not 2D, so it needs slightly different code
		DatTable[Count] = CurrentState
		Count +=1
	else:				#For J>0 we pull the values iteratively
		total += len(CurrentState)
		print(total)
		for j in range(len(CurrentState)):
			DatTable[Count] = CurrentState[j]
			Count +=1

#Write the merged values to a file
FileStr = ""
for i in range(len(DatTable)):
	for j in range(len(Delta)):
		FileStr += "%.16f\t" % DatTable[i][j]
	FileStr += "\n"
FileHandle = open ("J%d_%d_dk3.dat" % (JStart,JStop),'w')
FileHandle.write(FileStr)
FileHandle.close()
