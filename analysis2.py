#!/usr/bin/python2.7
# coding=UTF-8

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
from scipy import integrate
from scipy.optimize import minimize

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

correctedDirectory = r'./correctedGraphs/'
fitnessDirectory = r'./fitnessGraphs/'

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
	lowerYBound = Y[0]-tolerance*maxError
	upperYBound = Y[0]+tolerance*maxError
	
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

def fitQuadraticsToSurface(A,groupDataList):
	totalError = 0.0
	for group in groupDataList:
		#totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - A*x*group[2])**2,0.0,group[3])[0]
		#totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - A*x*group[2])**2/(group[1][0]**2*x**4+group[1][1]**2*x**2 + group[1][2]**2),0.0,group[3])[0]
		totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - A*x*group[2])**2/group[4]**2,0.0,group[3])[0]
	
	totalError *= 1.0/len(groupDataList)
	return totalError

# Runs done in constant Y
def getMaterialQuadratics(X,Y,Z,Xerr,Yerr,Zerr,materialName,XQuantity,XUnit,YQuantity,YUnit,rtol=100.0):
	runIndices = findRuns(Y,Yerr)
	print 'Found ' + str(len(runIndices)) + ' runs ignoring sign, indices:\n' + str(runIndices)
	runs = [ zip(*zip(*[X,Y,Z,Xerr,Yerr,Zerr])[a:b]) for a,b in runIndices ]
	positiveNegativeRuns = []
	for run in runs:
		# Sort run by X values
		sortedRun = [ list(i) for i in zip(*sorted(zip(*run),key=lambda tup: tup[0])) ]
		
		# Separate run into positive and negative part
		negatives,positives = separateSigns(sortedRun[0]) # Separate signs based on X values
		if(negatives[0] != None):
			positiveNegativeRuns.append(zip(*zip(*sortedRun)[negatives[0]:negatives[1]]))
		if(positives[0] != None):
			positiveNegativeRuns.append(zip(*zip(*sortedRun)[positives[0]:positives[1]]))
		#print 'Negative indices: ' + str(negatives) + ', positive indices: ' + str(positives)
	
	print 'Found ' + str(len(positiveNegativeRuns)) + ' signed runs.'
	
	# Group runs by abs(Y) value
	groupedRuns = []
	for run in positiveNegativeRuns:
		absAverageY = abs(np.mean(run[1]))
		maxErrorY = max(run[4])
		#print 'Looking for groups with |Y| = ' + str(absAverageY) + '±' + str(rtol*maxErrorY)
		foundMatch = False
		for g in groupedRuns:
			groupMeanY = np.mean([ abs(np.mean(r[1])) for r in g ])
			if(abs(groupMeanY - absAverageY) < rtol*maxErrorY):
				#print ' -> Group with |Y| = ' + str(groupMeanY) + ' matches!'
				g.append(run)
				foundMatch = True
				break
			#else:
				#print ' -> Group with |Y| = ' + str(groupMeanY) + ' doesn\'t match :('
		
		if(foundMatch == False):
			#print ' -> No matching group, creating new group.'
			groupedRuns.append([run])
	
	# Discard groups which don't have the right set of runs
	groupsToDiscard = []
	for i in range(len(groupedRuns)):
		fullSet = set([(1,1),(1,-1),(-1,1),(-1,-1)])
		currentSet = set()
		for run in groupedRuns[i]:
			if(np.mean(run[0]) > 0):
				if(np.mean(run[1]) > 0):
					currentSet.add((1,1))
				else:
					currentSet.add((1,-1))
			else:
				if(np.mean(run[1]) > 0):
					currentSet.add((-1,1))
				else:
					currentSet.add((-1,-1))
		
		if(fullSet == currentSet):
			print 'Group has correct set of runs'
		else:
			print 'Group doesn\'t have correct set of runs, deleting it'
			groupsToDiscard.append(i)
	
	for g in groupsToDiscard:
		groupedRuns.pop(g)
	
	print '\nFound the following groups:'
	for group in groupedRuns:
		groupMeanY = np.mean([ abs(np.mean(r[1])) for r in group ])
		print 'Group with |Y| = ' + str(groupMeanY) + ' contains ' + str(len(group)) + ' runs.'
		for i in range(len(group)):
			if(np.mean(group[i][0]) < 0.0):
				print ' - Run ' + str(i+1) + ', Z(-X,' + str(np.mean(group[i][1])) + '), ' + str(len(group[i][0])) + ' readings.'
			else:
				print ' - Run ' + str(i+1) + ', Z(+X,' + str(np.mean(group[i][1])) + '), ' + str(len(group[i][0])) + ' readings.'
	
	# Generate correct quadratic for each group
	groupCoefficientsList = []
	groupCoefficientErrorsList = []
	groupYValueList = []
	groupMaxXList = []
	groupErrorZList = []
	for group in groupedRuns:
		groupCoefficients = []
		groupCoefficientErrors = []
		groupYValues = []
		groupErrorZ = 0.0
		maxX = 0
		
		for run in group:
			transformedRun = copy.deepcopy(run)
			if(np.mean(run[0]) < 0):
				transformedRun[0] = [ -1.0*x for x in transformedRun[0] ]
			if(np.mean(run[1]) < 0):
				transformedRun[1] = [ -1.0*y for y in transformedRun[1] ]
			if((np.mean(run[0]) < 0) and (np.mean(run[1]) > 0)) or ((np.mean(run[0]) > 0) and (np.mean(run[1]) < 0)):
				transformedRun[2] = [ -1.0*z for z in transformedRun[2] ]
			
			if(max(transformedRun[0]) > maxX):
				maxX = max(transformedRun[0])
			
			groupYValues.append(list(transformedRun[1]))
			groupErrorZ = max(list(transformedRun[5]) + [groupErrorZ])
			
			# For numpy.polyfit the weights are 1/sigma, where sigma is standard deviation
			# Feed logged values into the fitting function to get coefficients and weirdly scaled covariance matrix
			fitResults = np.polyfit(x=transformedRun[0],y=transformedRun[2],deg=2,full=False,cov=True)
			coefficients = fitResults[0]
			covarianceMatrix = fitResults[-1]

			# Do this again to get residuals - for some reason function won't return residuals AND covariance matrix
			residuals = np.polyfit(x=transformedRun[0],y=transformedRun[2],deg=2,full=True,cov=False)[1][0]

			# the covariance matrix we have currently is scaled by residuals/(len(xs)-(deg+1)-2.0)
			# for some strange reason. Unscale it:
			covarianceMatrix = np.multiply(residuals/(len(transformedRun[0])-2.0-2.0),covarianceMatrix)
			#print covarianceMatrix
			
			groupCoefficients.append(coefficients)
			groupCoefficientErrors.append([ math.sqrt(covarianceMatrix[i][i]) for i in range(len(covarianceMatrix)) ])
		
		groupCoefficientsList.append([ np.mean(c) for c in zip(*groupCoefficients)])
		groupCoefficientErrorsList.append([ math.sqrt(np.sum(np.power(e,2)))/len(e) for e in zip(*groupCoefficientErrors) ])
		groupYValueList.append(np.mean([item for sublist in groupYValues for item in sublist]))
		groupMaxXList.append(maxX)
		groupErrorZList.append(groupErrorZ/math.sqrt(len(group)))
	
	for g in zip(groupCoefficientsList,groupCoefficientErrorsList,groupYValueList,groupMaxXList,groupErrorZList):
		print '\nFor group with |Y| = ' + str(g[2]) + ':'
		print ' - Coefficients are:\n' + str(g[0])
		print ' - Coefficient errors are:\n' + str(g[1])
		print ' - Maximum X values are: ' + str(g[3])
		print ' - Standard error is: ' + str(g[4])
	
	groupDataList = zip(groupCoefficientsList,groupCoefficientErrorsList,groupYValueList,groupMaxXList,groupErrorZList)
	positiveResult = minimize(lambda x: math.log10(fitQuadraticsToSurface(math.pow(10.0,x),groupDataList)), -1.0, method='nelder-mead')
	negativeResult = minimize(lambda x: math.log10(fitQuadraticsToSurface(-math.pow(10.0,x),groupDataList)),-1.0, method='nelder-mead')
	print 'Positive fit result: ' + str(positiveResult)
	print 'Negative fit result: ' + str(negativeResult)
	print '+R_H/t = ' + str(math.pow(10.0,positiveResult['x'])) + ', min error: ' + str(math.pow(10.0,positiveResult['fun']))
	print '-R_H/t = ' + str(-math.pow(10.0,negativeResult['x'])) + ', min error: ' + str(math.pow(10.0,negativeResult['fun']))
	
	rValues = np.logspace(-8,2,1000)
	totalErrors = []
	for r in rValues:
		totalError = 0.0
		for group in zip(groupCoefficientsList,groupCoefficientErrorsList,groupYValueList,groupMaxXList,groupErrorZList):
			totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - r*x*group[2])**2/group[4]**2,0.0,group[3])[0]
			#totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - r*x*group[2])**2/(group[1][0]**2*x**4+group[1][1]**2*x**2 + group[1][2]**2),0.0,group[3])[0]
			#totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - r*x*group[2])**2,0.0,group[3])[0]
		
		totalError *= 1.0/len(zip(groupCoefficientsList,groupCoefficientErrorsList,groupYValueList,groupMaxXList,groupErrorZList))
		totalErrors.append(totalError)
	
	# Plot positive coefficient fit result
	for group in zip(groupCoefficientsList,groupCoefficientErrorsList,groupYValueList,groupMaxXList,groupErrorZList):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		xSpace = np.linspace(0.0,group[3])
		yValsEmpirical,yMinVals,yMaxVals = quadraticCurveBounds(group[0],np.multiply(1e4,group[1]),xSpace)
		yMinVals = [ y-group[4] for y in yValsEmpirical ]
		yMaxVals = [ y+group[4] for y in yValsEmpirical ]
		yValsTheoretical = [ math.pow(10.0,positiveResult['x'][0])*x*group[2] for x in xSpace ]
		if(max([ abs(y) for y in yValsEmpirical]) < 1e-3):
			yValsEmpirical = np.multiply(yValsEmpirical,1e6)
			yMinVals = np.multiply(yMinVals,1e6)
			yMaxVals = np.multiply(yMaxVals,1e6)
			yValsTheoretical = np.multiply(yValsTheoretical,1e6)
			ax.set_ylabel(r'Hall Voltage /\SI{}{\micro\volt}')
		elif(max([ abs(y) for y in yValsEmpirical]) < 1.0):
			yValsEmpirical = np.multiply(yValsEmpirical,1e3)
			yMinVals = np.multiply(yMinVals,1e3)
			yMaxVals = np.multiply(yMaxVals,1e3)
			yValsTheoretical = np.multiply(yValsTheoretical,1e3)
			ax.set_ylabel(r'Hall Voltage /\SI{}{\milli\volt}')
		else:
			ax.set_ylabel(r'Hall Voltage /\SI{}{\volt}')
		ax.plot(xSpace,yValsTheoretical,'g-')
		ax.plot(xSpace,yValsEmpirical,'b-')
		ax.fill_between(xSpace, yMinVals, yMaxVals, facecolor='yellow',alpha=0.5,linestyle='--')
		ax.set_title(str(materialName) + ', ' + str(YQuantity) + '=' + str(group[2])[:6] + ' ' + str(YUnit))
		ax.set_xlabel(str(XQuantity) + ' /' + str(XUnit))
		ax.set_xlim(0.0,group[3])
		ax.grid()
		fig.savefig(str(correctedDirectory) + str(materialName) + '_corrected_' + str(group[2])[:6] + '_positive.pdf')
		
	# Plot negative coefficient fit result
	for group in zip(groupCoefficientsList,groupCoefficientErrorsList,groupYValueList,groupMaxXList,groupErrorZList):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		xSpace = np.linspace(0.0,group[3])
		yValsEmpirical,yMinVals,yMaxVals = quadraticCurveBounds(group[0],np.multiply(1e4,group[1]),xSpace)
		yMinVals = [ y-group[4] for y in yValsEmpirical ]
		yMaxVals = [ y+group[4] for y in yValsEmpirical ]
		yValsTheoretical = [ -math.pow(10.0,negativeResult['x'][0])*x*group[2] for x in xSpace ]
		if(max([ abs(y) for y in yValsEmpirical]) < 1e-3):
			yValsEmpirical = np.multiply(yValsEmpirical,1e6)
			yMinVals = np.multiply(yMinVals,1e6)
			yMaxVals = np.multiply(yMaxVals,1e6)
			yValsTheoretical = np.multiply(yValsTheoretical,1e6)
			ax.set_ylabel(r'Hall Voltage /\SI{}{\micro\volt}')
		elif(max([ abs(y) for y in yValsEmpirical]) < 1.0):
			yValsEmpirical = np.multiply(yValsEmpirical,1e3)
			yMinVals = np.multiply(yMinVals,1e3)
			yMaxVals = np.multiply(yMaxVals,1e3)
			yValsTheoretical = np.multiply(yValsTheoretical,1e3)
			ax.set_ylabel(r'Hall Voltage /\SI{}{\milli\volt}')
		else:
			ax.set_ylabel(r'Hall Voltage /\SI{}{\volt}')
		ax.plot(xSpace,yValsTheoretical,'g-')
		ax.plot(xSpace,yValsEmpirical,'b-')
		ax.fill_between(xSpace, yMinVals, yMaxVals, facecolor='yellow',alpha=0.5,linestyle='--')
		ax.set_title(str(materialName) + ', ' + str(YQuantity) + '=' + str(group[2])[:6] + ' ' + str(YUnit))
		ax.set_xlabel(str(XQuantity) + ' /' + str(XUnit))
		ax.set_xlim(0.0,group[3])
		ax.grid()
		fig.savefig(str(correctedDirectory) + str(materialName) + '_corrected_' + str(group[2])[:6] + '_negative.pdf')
	
	# Plot positive coefficient errors vs. R_H
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(rValues,totalErrors)
	ax.plot(math.pow(10.0,positiveResult['x'][0]),math.pow(10.0,positiveResult['fun']),'ro')
	ax = plt.axes()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel(r'$R_H/t$ /\SI{}{\square\metre\per\coulomb}')
	ax.set_ylabel('Integral of squared errors',rotation=90)
	ax.grid()
	fig.savefig(str(fitnessDirectory) + str(materialName) + '_errors_positive.pdf')
	
	rValuesNegative = np.dot(-1.0,np.logspace(-8,2,1000))
	totalErrors = []
	for r in rValuesNegative:
		totalError = 0.0
		for group in zip(groupCoefficientsList,groupCoefficientErrorsList,groupYValueList,groupMaxXList,groupErrorZList):
			totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - r*x*group[2])**2/group[4]**2,0.0,group[3])[0]
			#totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - r*x*group[2])**2/(group[1][0]**2*x**4+group[1][1]**2*x**2 + group[1][2]**2),0.0,group[3])[0]
			#totalError += (1.0/group[3])*integrate.quad(lambda x: (group[0][0]*x**2 + group[0][1]*x + group[0][2] - r*x*group[2])**2,0.0,group[3])[0]
		
		totalError *= 1.0/len(zip(groupCoefficientsList,groupCoefficientErrorsList,groupYValueList,groupMaxXList,groupErrorZList))
		totalErrors.append(totalError)
	
	# Plot negative coefficient errors vs. R_H
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(rValues,totalErrors)
	ax.plot(math.pow(10.0,negativeResult['x'][0]),math.pow(10.0,negativeResult['fun']),'ro')
	ax = plt.axes()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel(r'$-R_H/t$ /\SI{}{\square\metre\per\coulomb}')
	ax.set_ylabel('Integral of squared errors',rotation=90)
	ax.grid()
	fig.savefig(str(fitnessDirectory) + str(materialName) + '_errors_negative.pdf')
	
	print '+R_H/t = ' + str(math.pow(10.0,positiveResult['x'])) + ', min error: ' + str(math.pow(10.0,positiveResult['fun']))
	print '-R_H/t = ' + str(-math.pow(10.0,negativeResult['x'])) + ', min error: ' + str(math.pow(10.0,negativeResult['fun']))
	
	possibleHallConstant = [math.pow(10.0,positiveResult['x']),-math.pow(10.0,negativeResult['x'])]
	functionValues = [math.pow(10.0,positiveResult['fun']),math.pow(10.0,negativeResult['fun'])]
	return possibleHallConstant[functionValues.index(min(functionValues))],min(functionValues)
	
# results in format: [ [B /T], [B error], [I /A], [I error], [T /C], [T error /C], [V /V], [V error /V] ]
results = np.load('results.npz')

pTypeGe = results['pTypeGe'] # I fixed
nTypeGe = results['nTypeGe'] # I fixed
tungsten = results['tungsten'] # B fixed
silver = results['silver'] # B fixed

pTypeGeThickness = 1e-3
nTypeGeThickness = 1e-3
tungstenThickness = 5e-5
silverThickness = 5e-5

print '\n'+\
      '********************************************************************************'+\
      '                             p-Type Ge analysis                                 '+\
      '********************************************************************************'
possibleHallConstant,error = getMaterialQuadratics(pTypeGe[0],pTypeGe[2],pTypeGe[6],pTypeGe[1],pTypeGe[3],pTypeGe[7],'p-type Ge','B','T','I','A')
print 'Hall Constant (p-type Ge): ' + str(pTypeGeThickness*possibleHallConstant) + ', literature value: 6.6*10^-3, min error: ' + str(error)

print '\n'+\
      '********************************************************************************'+\
      '                             n-Type Ge analysis                                 '+\
      '********************************************************************************'
possibleHallConstant,error = getMaterialQuadratics(nTypeGe[0],nTypeGe[2],nTypeGe[6],nTypeGe[1],nTypeGe[3],nTypeGe[7],'n-type Ge','B','T','I','A')
print 'Hall Constant (n-type Ge): ' + str(nTypeGeThickness*possibleHallConstant) + ', literature value: 5.6*10^-3, min error: ' + str(error)

print '\n'+\
      '********************************************************************************'+\
      '                              Tungsten analysis                                 '+\
      '********************************************************************************'
possibleHallConstant,error = getMaterialQuadratics(tungsten[2],tungsten[0],tungsten[6],tungsten[3],tungsten[1],tungsten[7],'Tungsten','I','A','B','T',rtol=10.0)
print 'Hall Constant (tungsten): ' + str(tungstenThickness*possibleHallConstant) + ', literature value: 1.18*10^-10, min error: ' + str(error)

print '\n'+\
      '********************************************************************************'+\
      '                               Silver analysis                                  '+\
      '********************************************************************************'
possibleHallConstant,error = getMaterialQuadratics(silver[2],silver[0],silver[6],silver[3],silver[1],silver[7],'Silver','I','A','B','T',rtol=10.0)
print 'Hall Constant (silver): ' + str(silverThickness*possibleHallConstant) + ', literature value: 8.9*10^-11, min error: ' + str(error)

