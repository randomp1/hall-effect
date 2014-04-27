#!/usr/bin/python2.7
# coding=UTF-8

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=10.0)
plt.rc('figure', autolayout=True)
plt.rc('text.latex', preamble = ','.join('''
 \usepackage{siunitx}
 '''.split()))
# Pts to inches conversion factor
pointsToInches=1/72.27
widthPts = 450.0
widthInches = widthPts*pointsToInches
heightInches = widthInches/1.61803398875
plt.rc('figure', figsize=(widthInches,heightInches))
outputFolder = r'./graphs_raw/'

# Finds bounds on a quadratic given its coefficients and standard errors
def quadraticCurveBounds(coefficients,coeffErrors,xVals):
	yVals = [ coefficients[0]*x**2 + coefficients[1]*x + coefficients[2] for x in xVals ]
	yMaxVals = list(yVals)
	yMinVals = list(yVals)
	for i in range(len(xVals)):
		if(xVals[i] >= 0.0):
			yMaxVals[i] += coeffErrors[0]*xVals[i]**2 + coeffErrors[1]*xVals[i] + coeffErrors[2]
			yMinVals[i] -= coeffErrors[0]*xVals[i]**2 + coeffErrors[1]*xVals[i] + coeffErrors[2]
		else:
			yMaxVals[i] += coeffErrors[0]*xVals[i]**2 - coeffErrors[1]*xVals[i] + coeffErrors[2]
			yMinVals[i] -= coeffErrors[0]*xVals[i]**2 - coeffErrors[1]*xVals[i] + coeffErrors[2]
	
	return yVals,yMinVals,yMaxVals

# Finds xs in a list of x values with negative and positive signs
def separateSigns(X):
	positiveStart = None
	positiveEnd = None
	negativeStart = None
	negativeEnd = None
	
	if(X[0] < 0.0):
		negativeStart = 0
	else:
		positiveStart = 0
	
	if(X[-1] < 0.0) and (negativeStart == 0):
		negativeEnd = len(X)
	elif(X[-1] >= 0.0) and (positiveStart == 0):
		positiveEnd = len(X)
	else:
		for i in range(len(X)):
			if(X[i] >= 0.0 and X[0] < 0.0):
				negativeEnd = i
				if(i != len(X)-1):
					positiveStart = i
					positiveEnd = len(X)
				break
			if(X[i] < 0.0 and X[0] >= 0.0):
				positiveEnd = i
				if(i != len(X)-1):
					negativeStart = i
					negativeEnd = len(X)
				break
	
	negatives = (negativeStart,negativeEnd)
	positives = (positiveStart,positiveEnd)

	return negatives,positives

def findRuns(Y,Yerr):
	# Examining parameter space of x-y, with y held fixed at values from y_min to y_max
	# while x is varied for a single 'run'. We consider a 'run' to have y-value fixed
	# to within y±(y_err*2), where y_err is the maximum y error found
	
	runIndices = [] # List of form [(run start index,run end index+1),(etc.)]
	startIndex = 0
	tolerance = 10.0
	
	maxError = max(Yerr)
	lowerYBound = Y[0]-2.0*maxError
	upperYBound = Y[0]+2.0*maxError
	
	for i in range(len(Y)):
		if(Y[i] < lowerYBound) or (Y[i] > upperYBound):
			# Not part of the current run
			runIndices.append((startIndex,i))
			#print 'Y[' + str(runIndices[-1][0]) + '] = ' + str(Y[runIndices[-1][0]])
			#print 'Y[' + str(runIndices[-1][1]-1) + '] = ' + str(Y[runIndices[-1][1]-1])
			startIndex = i
			lowerYBound = Y[i]-tolerance*maxError
			upperYBound = Y[i]+tolerance*maxError
		elif (Y[i] - tolerance*maxError < lowerYBound):
			lowerYBound = Y[i]-tolerance*maxError
		elif (Y[i] + tolerance*maxError > upperYBound):
			upperYBound = Y[i]+tolerance*maxError
	
	runIndices.append((startIndex,len(Y)))
	
	return runIndices

# Wrapper to polyfit
def getCoeffsAndCovMatrix(X,Y,errY,degree):
	# For numpy.polyfit the weights are 1/sigma, where sigma is standard deviation
	# Feed logged values into the fitting function to get coefficients and weirdly scaled covariance matrix
	fitResults = np.polyfit(x=X,y=Y,deg=degree,full=False,cov=True)
	coefficients = fitResults[0]
	covarianceMatrix = fitResults[-1]

	# Do this again to get residuals - for some reason function won't return residuals AND covariance matrix
	residuals = np.polyfit(x=X,y=Y,deg=degree,full=True,cov=False)[1][0]

	# the covariance matrix we have currently is scaled by residuals/(len(xs)-(deg+1)-2.0)
	# for some strange reason. Unscale it:
	covarianceMatrix = np.multiply(residuals/(len(X)-(degree+1.0)-2.0),covarianceMatrix)
	
	return coefficients,covarianceMatrix

# results in format: [ [B /T], [B error], [I /A], [I error], [T /C], [T error /C], [V /V], [V error /V] ]
results = np.load('results.npz')

pTypeGe = results['pTypeGe']
nTypeGe = results['nTypeGe']
tungsten = results['tungsten']
silver = results['silver']

print '\n'+\
      '********************************************************************************'+\
      '                             p-Type Ge analysis                                 '+\
      '********************************************************************************'
# Runs done in constant I
pTypeGeRunIndices = findRuns(pTypeGe[2],pTypeGe[3])
print pTypeGeRunIndices
pTypeGeRuns = [ zip(*zip(*pTypeGe)[a:b]) for a,b in pTypeGeRunIndices ]

for run in pTypeGeRuns:
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.errorbar(run[0],run[6],yerr=run[7],fmt='bo')
	ax.set_xlim(-0.3,0.3)
	
	negatives,positives = separateSigns(run[0])
	if(positives[0] != None):
		positiveCoefficients,posCovMatrix = getCoeffsAndCovMatrix(run[0][positives[0]:positives[1]],run[6][positives[0]:positives[1]],run[7][positives[0]:positives[1]],2)
		print 'Positive coefficients: ' + str(positiveCoefficients)
		print 'Errors: ' + str([ math.sqrt(posCovMatrix[i][i]) for i in range(len(posCovMatrix)) ])
		positiveSpace = np.linspace(0.0,ax.get_xlim()[1])
		positiveValues,posLowerBound,posUpperBound = quadraticCurveBounds(positiveCoefficients,[ math.sqrt(posCovMatrix[i][i]) for i in range(len(posCovMatrix)) ],positiveSpace)
		ax.plot(positiveSpace,positiveValues,'g-')
		#ax.fill_between(positiveSpace, posLowerBound, posUpperBound, facecolor='yellow',alpha=0.5,linestyle='--')
	if(negatives[0] != None):
		negativeCoefficients,negCovMatrix = getCoeffsAndCovMatrix(run[0][negatives[0]:negatives[1]],run[6][negatives[0]:negatives[1]],run[7][negatives[0]:negatives[1]],2)
		print 'Negative coefficients: ' + str(negativeCoefficients)
		print 'Errors: ' + str([ math.sqrt(negCovMatrix[i][i]) for i in range(len(negCovMatrix)) ])
		negativeSpace = np.linspace(ax.get_xlim()[0],0.0)
		negativeValues,negLowerBound,negUpperBound = quadraticCurveBounds(negativeCoefficients,[ math.sqrt(negCovMatrix[i][i]) for i in range(len(negCovMatrix)) ],negativeSpace)
		ax.plot(negativeSpace,negativeValues,'r-')
		#ax.fill_between(negativeSpace, negLowerBound, negUpperBound, facecolor='yellow',alpha=0.5,linestyle='--')
	
	ax.set_title('p-type Ge, I=' + str(run[2][0])[:6] + 'A')
	ax.set_xlabel('B /T')
	ax.set_ylabel('Hall voltage /V')
	print 'Max error: ' + str(max(run[7]))
	ax.grid()
	fig.savefig(str(outputFolder) + 'pTypeGe_' + str(run[2][0]) + '.pdf')

runCoefficientsList = []
runCoefficientErrorsList = []
transformedDataList = []
for run in pTypeGeRuns:
	negatives,positives = separateSigns(run[0])
	if(positives[0] != None):
		fitFields = [ b for b in run[0][positives[0]:positives[1]] ]
		# Positive B, check sign of I
		if(np.mean(run[2]) < 0.0):
			fitVoltages = [ -v for v in run[6][positives[0]:positives[1]] ]
			print '(-I;+B): B\'=B, V\'=-V, I\'= ' + str(abs(np.mean(run[2]))) + ', error= ' + str(max(run[3]))
		else:
			fitVoltages = [ v for v in run[6][positives[0]:positives[1]] ]
			print '(+I;+B): B\'=B, V\'=V, I\'= ' + str(abs(np.mean(run[2]))) + '±' + str(max(run[3]))
		# For numpy.polyfit the weights are 1/sigma, where sigma is standard deviation
		# Feed logged values into the fitting function to get coefficients and weirdly scaled covariance matrix
		fitResults = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=False,cov=True)
		coefficients = fitResults[0]
		covarianceMatrix = fitResults[-1]

		# Do this again to get residuals - for some reason function won't return residuals AND covariance matrix
		residuals = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=True,cov=False)[1][0]

		# the covariance matrix we have currently is scaled by residuals/(len(xs)-(deg+1)-2.0)
		# for some strange reason. Unscale it:
		covarianceMatrix = np.multiply(residuals/(len(fitFields)-2.0-2.0),covarianceMatrix)
		runCoefficientErrorsList.append([ math.sqrt(covarianceMatrix[i][i]) for i in range(len(covarianceMatrix)) ])
		runCoefficientsList.append(coefficients)
		print 'Positive coefficients: ' + str(coefficients[0]) + ', ' + str(coefficients[1]) + ', ' + str(coefficients[2]) + '\nerrors: ' + str(runCoefficientErrorsList[-1][0]) + ', ' + str(runCoefficientErrorsList[-1][1]) + ', ' + str(runCoefficientErrorsList[-1][2])
	if(negatives[0] != None):
		fitFields = [ -b for b in run[0][negatives[0]:negatives[1]] ]
		# Negative B, check sign of I
		if(np.mean(run[2]) < 0.0):
			fitVoltages = [ v for v in run[6][negatives[0]:negatives[1]] ]
			print '(-I,-B): B\'=-B, V\'=V, I\'= ' + str(abs(np.mean(run[2]))) + '±' + str(max(run[3]))
		else:
			fitVoltages = [ -v for v in run[6][negatives[0]:negatives[1]] ]
			print '(+I,-B): B\'=-B, V\'=-V, I\'= ' + str(abs(np.mean(run[2]))) + '±' + str(max(run[3]))
		# For numpy.polyfit the weights are 1/sigma, where sigma is standard deviation
		# Feed logged values into the fitting function to get coefficients and weirdly scaled covariance matrix
		fitResults = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=False,cov=True)
		coefficients = fitResults[0]
		covarianceMatrix = fitResults[-1]

		# Do this again to get residuals - for some reason function won't return residuals AND covariance matrix
		residuals = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=True,cov=False)[1][0]

		# the covariance matrix we have currently is scaled by residuals/(len(xs)-(deg+1)-2.0)
		# for some strange reason. Unscale it:
		covarianceMatrix = np.multiply(residuals/(len(fitFields)-2.0-2.0),covarianceMatrix)
		runCoefficientErrorsList.append([ math.sqrt(covarianceMatrix[i][i]) for i in range(len(covarianceMatrix)) ])
		runCoefficientsList.append(coefficients)
		print 'Negative coefficients: ' + str(coefficients[0]) + ', ' + str(coefficients[1]) + ', ' + str(coefficients[2]) + '\nerrors: ' + str(runCoefficientErrorsList[-1][0]) + ', ' + str(runCoefficientErrorsList[-1][1]) + ', ' + str(runCoefficientErrorsList[-1][2])

combinedCoefficients = [ np.mean(x) for x in zip(*runCoefficientsList) ]
fig = plt.figure()
ax = fig.add_subplot(111)
combinedSpace = np.linspace(0.0,0.3)
combinedValues = [ combinedCoefficients[0]*(x**2) + combinedCoefficients[1]*x + combinedCoefficients[2] for x in combinedSpace ]
for data in transformedDataList:
	ax.plot(data[0],data[1],'o')
ax.plot(combinedSpace,combinedValues)

ax.set_title('p-type Ge, combined')
ax.set_xlabel('B /T')
ax.set_ylabel('Hall Voltage /V')
ax.grid()
fig.show()


print 'Combined coefficients: ' + str(combinedCoefficients)

hallCoefficient = -1e-4
thickness = 5e-5
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.set_title('p-type Ge')
ax.set_xlabel(r'B /T')
ax.set_ylabel(r'I /A')
ax.set_zlabel(r'Hall Voltage /V')
z = pTypeGe[6]
x = pTypeGe[0]
y = pTypeGe[2]
ax.scatter(x, y, z)
#minB = ax.get_xlim()[0]
#maxB = ax.get_xlim()[1]
#minI = ax.get_ylim()[0]
#maxI = ax.get_ylim()[1]
#x = np.linspace(*ax.get_xlim(),num=3)
#y = np.linspace(*ax.get_ylim(),num=3)
#x,y = np.meshgrid(x,y)
#x = np.reshape(x,3*3)
#y = np.reshape(y,3*3)
#z = [ hallCoefficient*I*B/thickness for I,B in zip(x,y) ]
#print x
#print y
#print z
#ax.plot_surface(np.reshape(x,(3,3)), np.reshape(y,(3,3)), np.reshape(z,(3,3)), rstride=1, cstride=1, linewidth=0,shade=True)
fig.show()

print '\n'+\
      '********************************************************************************'+\
      '                             n-Type Ge analysis                                 '+\
      '********************************************************************************'

nTypeGeRunIndices = findRuns(nTypeGe[2],nTypeGe[3])
print nTypeGeRunIndices
nTypeGeRuns = [ zip(*zip(*nTypeGe)[a:b]) for a,b in nTypeGeRunIndices ]

for run in nTypeGeRuns:
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.errorbar(run[0],run[6],yerr=run[7],fmt='bo')
	ax.set_xlim(-0.3,0.3)
	
	negatives,positives = separateSigns(run[0])
	if(positives[0] != None):
		positiveCoefficients = np.polyfit(run[0][positives[0]:positives[1]],run[6][positives[0]:positives[1]],2)
		print 'Positive coefficients: ' + str(positiveCoefficients)
		positiveSpace = np.linspace(0.0,ax.get_xlim()[1])
		positiveValues = [ positiveCoefficients[0]*(x**2) + positiveCoefficients[1]*x + positiveCoefficients[2] for x in positiveSpace ]
		ax.plot(positiveSpace,positiveValues,'g-')
	if(negatives[0] != None):
		negativeCoefficients = np.polyfit(run[0][negatives[0]:negatives[1]],run[6][negatives[0]:negatives[1]],2)
		print 'Negative coefficients: ' + str(negativeCoefficients)
		negativeSpace = np.linspace(ax.get_xlim()[0],0.0)
		negativeValues = [ negativeCoefficients[0]*(x**2) + negativeCoefficients[1]*x + negativeCoefficients[2] for x in negativeSpace ]
		ax.plot(negativeSpace,negativeValues,'r-')
	
	ax.set_title('n-type Ge, I=' + str(run[2][0])[:6] + 'A')
	ax.set_xlabel('B /T')
	ax.set_ylabel('Hall Voltage /V')
	print 'Max error: ' + str(max(run[7]))
	ax.grid()
	fig.savefig(str(outputFolder) + 'nTypeGe_' + str(run[2][0]) + '.pdf')


runCoefficientsList = []
runCoefficientErrorsList = []
transformedDataList = []
for run in nTypeGeRuns:
	negatives,positives = separateSigns(run[0])
	if(positives[0] != None):
		fitFields = [ b for b in run[0][positives[0]:positives[1]] ]
		# Positive B, check sign of I
		if(np.mean(run[2]) < 0.0):
			fitVoltages = [ -v for v in run[6][positives[0]:positives[1]] ]
			print '(-I,+B): B\'=B, V\'=-V, I\'= ' + str(abs(np.mean(run[2]))) + '±' + str(max(run[3]))
		else:
			fitVoltages = [ v for v in run[6][positives[0]:positives[1]] ]
			print '(+I,+B): B\'=B, V\'=V, I\'= ' + str(abs(np.mean(run[2]))) + '±' + str(max(run[3]))
		# For numpy.polyfit the weights are 1/sigma, where sigma is standard deviation
		# Feed logged values into the fitting function to get coefficients and weirdly scaled covariance matrix
		fitResults = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=False,cov=True)
		coefficients = fitResults[0]
		covarianceMatrix = fitResults[-1]

		# Do this again to get residuals - for some reason function won't return residuals AND covariance matrix
		residuals = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=True,cov=False)[1][0]

		# the covariance matrix we have currently is scaled by residuals/(len(xs)-(deg+1)-2.0)
		# for some strange reason. Unscale it:
		covarianceMatrix = np.multiply(residuals/(len(fitFields)-2.0-2.0),covarianceMatrix)
		runCoefficientErrorsList.append([ math.sqrt(covarianceMatrix[i][i]) for i in range(len(covarianceMatrix)) ])
		runCoefficientsList.append(coefficients)
		print 'Positive coefficients: ' + str(coefficients[0]) + ', ' + str(coefficients[1]) + ', ' + str(coefficients[2]) + '\nerrors: ' + str(runCoefficientErrorsList[-1][0]) + ', ' + str(runCoefficientErrorsList[-1][1]) + ', ' + str(runCoefficientErrorsList[-1][2])
	if(negatives[0] != None):
		fitFields = [ -b for b in run[0][negatives[0]:negatives[1]] ]
		# Negative B, check sign of I
		if(np.mean(run[2]) < 0.0):
			fitVoltages = [ v for v in run[6][negatives[0]:negatives[1]] ]
			print '(-I,-B): B\'=-B, V\'=V, I\'= ' + str(abs(np.mean(run[2]))) + '±' + str(max(run[3]))
		else:
			fitVoltages = [ -v for v in run[6][negatives[0]:negatives[1]] ]
			print '(+I,-B): B\'=-B, V\'=-V, I\'= ' + str(abs(np.mean(run[2]))) + '±' + str(max(run[3]))
		# For numpy.polyfit the weights are 1/sigma, where sigma is standard deviation
		# Feed logged values into the fitting function to get coefficients and weirdly scaled covariance matrix
		fitResults = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=False,cov=True)
		coefficients = fitResults[0]
		covarianceMatrix = fitResults[-1]

		# Do this again to get residuals - for some reason function won't return residuals AND covariance matrix
		residuals = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=True,cov=False)[1][0]

		# the covariance matrix we have currently is scaled by residuals/(len(xs)-(deg+1)-2.0)
		# for some strange reason. Unscale it:
		covarianceMatrix = np.multiply(residuals/(len(fitFields)-2.0-2.0),covarianceMatrix)
		runCoefficientErrorsList.append([ math.sqrt(covarianceMatrix[i][i]) for i in range(len(covarianceMatrix)) ])
		runCoefficientsList.append(coefficients)
		print 'Negative coefficients: ' + str(coefficients[0]) + ', ' + str(coefficients[1]) + ', ' + str(coefficients[2]) + '\nerrors: ' + str(runCoefficientErrorsList[-1][0]) + ', ' + str(runCoefficientErrorsList[-1][1]) + ', ' + str(runCoefficientErrorsList[-1][2])



fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.set_title('n-type Ge')
ax.set_xlabel(r'B /T')
ax.set_ylabel(r'I /A')
ax.set_zlabel(r'Hall Voltage /V')
z = nTypeGe[6]
x = nTypeGe[0]
y = nTypeGe[2]
ax.scatter(x, y, z)
fig.show()

print '\n'+\
      '********************************************************************************'+\
      '                              tungsten analysis                                 '+\
      '********************************************************************************'

tungstenRunIndices = findRuns(tungsten[0],tungsten[1])
print tungstenRunIndices
tungstenRuns = [ zip(*zip(*tungsten)[a:b]) for a,b in tungstenRunIndices ]

for run in tungstenRuns:
	sortedRun = [ list(i) for i in zip(*sorted(zip(*run),key=lambda tup: tup[2])) ]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.errorbar(sortedRun[2],[ s*1e6 for s in sortedRun[6] ],yerr=[ e*1e6 for e in sortedRun[7] ],fmt='bo')
	ax.set_xlim(-12,12)
	
	negatives,positives = separateSigns(sortedRun[2])
	if(positives[0] != None):
		positiveCoefficients = np.polyfit(sortedRun[2][positives[0]:positives[1]],sortedRun[6][positives[0]:positives[1]],2)
		print 'Positive coefficients: ' + str(positiveCoefficients)
		positiveSpace = np.linspace(0.0,ax.get_xlim()[1])
		positiveValues = [ y*1e6 for y in [ positiveCoefficients[0]*(x**2) + positiveCoefficients[1]*x + positiveCoefficients[2] for x in positiveSpace ] ]
		ax.plot(positiveSpace,positiveValues,'g-')
	if(negatives[0] != None):
		negativeCoefficients = np.polyfit(sortedRun[2][negatives[0]:negatives[1]],sortedRun[6][negatives[0]:negatives[1]],2)
		print 'Negative coefficients: ' + str(negativeCoefficients)
		negativeSpace = np.linspace(ax.get_xlim()[0],0.0)
		negativeValues = [ y*1e6 for y in [ negativeCoefficients[0]*(x**2) + negativeCoefficients[1]*x + negativeCoefficients[2] for x in negativeSpace ] ]
		ax.plot(negativeSpace,negativeValues,'r-')
	
	ax.set_title('Tungsten, B=' + str(run[0][0])[:6] + 'T')
	ax.set_xlabel('I /A')
	ax.set_ylabel(r'Hall Voltage /\SI{}{\micro\volt}')
	print 'Max error: ' + str(max(run[7]))
	ax.grid()
	fig.savefig(str(outputFolder) + 'tungsten_' + str(run[0][0]) + '.pdf')

runCoefficientsList = []
runCoefficientErrorsList = []
transformedDataList = []
for run in tungstenRuns:
	sortedRun = [ list(i) for i in zip(*sorted(zip(*run),key=lambda tup: tup[2])) ]
	negatives,positives = separateSigns(sortedRun[2])
	print negatives
	print positives
	if(positives[0] != None):
		fitFields = [ i for i in sortedRun[2][positives[0]:positives[1]] ]
		# Positive I, check sign of B
		if(np.mean(sortedRun[0]) < 0.0):
			fitVoltages = [ -v for v in sortedRun[6][positives[0]:positives[1]] ]
			print '(-B,+I): I\'=I, V\'=-V, B\'= ' + str(abs(np.mean(sortedRun[0]))) + '±' + str(max(sortedRun[1]))
		else:
			fitVoltages = [ v for v in sortedRun[6][positives[0]:positives[1]] ]
			print '(+B,+I): I\'=I, V\'=V, B\'= ' + str(abs(np.mean(sortedRun[0]))) + '±' + str(max(sortedRun[1]))
		# For numpy.polyfit the weights are 1/sigma, where sigma is standard deviation
		# Feed logged values into the fitting function to get coefficients and weirdly scaled covariance matrix
		fitResults = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=False,cov=True)
		coefficients = fitResults[0]
		covarianceMatrix = fitResults[-1]

		# Do this again to get residuals - for some reason function won't return residuals AND covariance matrix
		residuals = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=True,cov=False)[1][0]

		# the covariance matrix we have currently is scaled by residuals/(len(xs)-(deg+1)-2.0)
		# for some strange reason. Unscale it:
		covarianceMatrix = np.multiply(residuals/(len(fitFields)-2.0-2.0),covarianceMatrix)
		runCoefficientErrorsList.append([ math.sqrt(covarianceMatrix[i][i]) for i in range(len(covarianceMatrix)) ])
		runCoefficientsList.append(coefficients)
		print 'Positive coefficients: ' + str(coefficients[0]) + ', ' + str(coefficients[1]) + ', ' + str(coefficients[2]) + '\nerrors: ' + str(runCoefficientErrorsList[-1][0]) + ', ' + str(runCoefficientErrorsList[-1][1]) + ', ' + str(runCoefficientErrorsList[-1][2])
	if(negatives[0] != None):
		fitFields = [ -i for i in sortedRun[2][negatives[0]:negatives[1]] ]
		# Negative I, check sign of B
		if(np.mean(sortedRun[0]) < 0.0):
			fitVoltages = [ v for v in sortedRun[6][negatives[0]:negatives[1]] ]
			print '(-B,-I): I\'=-I, V\'=V, B\'= ' + str(abs(np.mean(sortedRun[0]))) + '±' + str(max(sortedRun[1]))
		else:
			fitVoltages = [ -v for v in sortedRun[6][negatives[0]:negatives[1]] ]
			print '(+B,-I): I\'=-I, V\'=-V, B\'= ' + str(abs(np.mean(sortedRun[0]))) + '±' + str(max(sortedRun[1]))
		# For numpy.polyfit the weights are 1/sigma, where sigma is standard deviation
		# Feed logged values into the fitting function to get coefficients and weirdly scaled covariance matrix
		fitResults = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=False,cov=True)
		coefficients = fitResults[0]
		covarianceMatrix = fitResults[-1]

		# Do this again to get residuals - for some reason function won't return residuals AND covariance matrix
		residuals = np.polyfit(x=fitFields,y=fitVoltages,deg=2,full=True,cov=False)[1][0]

		# the covariance matrix we have currently is scaled by residuals/(len(xs)-(deg+1)-2.0)
		# for some strange reason. Unscale it:
		covarianceMatrix = np.multiply(residuals/(len(fitFields)-2.0-2.0),covarianceMatrix)
		runCoefficientErrorsList.append([ math.sqrt(covarianceMatrix[i][i]) for i in range(len(covarianceMatrix)) ])
		runCoefficientsList.append(coefficients)
		print 'Negative coefficients: ' + str(coefficients[0]) + ', ' + str(coefficients[1]) + ', ' + str(coefficients[2]) + '\nerrors: ' + str(runCoefficientErrorsList[-1][0]) + ', ' + str(runCoefficientErrorsList[-1][1]) + ', ' + str(runCoefficientErrorsList[-1][2])

print '\n'+\
      '********************************************************************************'+\
      '                                Silver analysis                                 '+\
      '********************************************************************************'
silverRunIndices = findRuns(silver[0],silver[1])
print silverRunIndices
silverRuns = [ zip(*zip(*silver)[a:b]) for a,b in silverRunIndices ]

for run in silverRuns:
	sortedRun = [ list(i) for i in zip(*sorted(zip(*run),key=lambda tup: tup[2])) ]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.errorbar(sortedRun[2],[ s*1e6 for s in sortedRun[6] ],yerr=[ e*1e6 for s in sortedRun[7] ],fmt='bo')
	ax.set_xlim(-12,12)
	
	negatives,positives = separateSigns(sortedRun[2])
	if(positives[0] != None):
		positiveCoefficients = np.polyfit(sortedRun[2][positives[0]:positives[1]],sortedRun[6][positives[0]:positives[1]],2)
		print 'Positive coefficients: ' + str(positiveCoefficients)
		positiveSpace = np.linspace(0.0,ax.get_xlim()[1])
		positiveValues = [ s*1e6 for s in [ positiveCoefficients[0]*(x**2) + positiveCoefficients[1]*x + positiveCoefficients[2] for x in positiveSpace ] ]
		ax.plot(positiveSpace,positiveValues,'g-')
	if(negatives[0] != None):
		negativeCoefficients = np.polyfit(sortedRun[2][negatives[0]:negatives[1]],sortedRun[6][negatives[0]:negatives[1]],2)
		print 'Negative coefficients: ' + str(negativeCoefficients)
		negativeSpace = np.linspace(ax.get_xlim()[0],0.0)
		negativeValues = [ s*1e6 for s in [ negativeCoefficients[0]*(x**2) + negativeCoefficients[1]*x + negativeCoefficients[2] for x in negativeSpace ] ]
		ax.plot(negativeSpace,negativeValues,'r-')
	
	ax.set_title('Silver, B=' + str(run[0][0])[:6] + 'T')
	ax.set_xlabel('I /A')
	ax.set_ylabel(r'Hall Voltage /\SI{}{\micro\volt}')
	print 'Max error: ' + str(max(run[7]))
	ax.grid()
	fig.savefig(str(outputFolder) + 'silver_' + str(run[0][0]) + '.pdf')

'''
nTypeGeRuns = findRuns(nTypeGe[2],nTypeGe[3])
tungstenRuns = findRuns(tungsten[0],tungsten[1])
silverRuns = findRuns(silver[0],silver[1])

print tungstenRuns
print silverRuns

# Process p-type Germanium sample results
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_title('p-type Ge')
ax.set_xlabel(r'B /T')
ax.set_ylabel(r'I /A')
ax.set_zlabel(r'Hall Voltage /V')
z = pTypeGe[6]
x = pTypeGe[0]
y = pTypeGe[2]
ax.scatter(x, y, z)
plt.show()


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_title('n-type Ge')
ax.set_xlabel(r'B /T')
ax.set_ylabel(r'I /A')
ax.set_zlabel(r'Hall Voltage /V')
z = nTypeGe[6]
x = nTypeGe[0]
y = nTypeGe[2]
ax.scatter(x, y, z)
plt.show()
'''
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.set_title('Tungsten')
ax.set_xlabel(r'B /T')
ax.set_ylabel(r'I /A')
ax.set_zlabel(r'Hall Voltage /V')
z = tungsten[6]
x = tungsten[0]
y = tungsten[2]
ax.scatter(x, y, z)
fig.show()

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.set_title('Silver')
ax.set_xlabel(r'B /T')
ax.set_ylabel(r'I /A')
ax.set_zlabel(r'Hall Voltage /V')
z = silver[6]
x = silver[0]
y = silver[2]
ax.scatter(x, y, z)
fig.show()

