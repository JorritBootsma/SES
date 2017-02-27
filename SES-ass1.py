######################################################################
#   Structure and Evolution of Stars - SES-ass1.py
#   Date: 		20/02/2017
#   Name: 		Jorrit Bootsma
#	Studentid: 	10251499
######################################################################

from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spInt
from sympy import *
print ''

def ode(y, ksi, n, a):
	'''
	The set of differential equations that needs to be solved.
	@param y: list of functions that one takes the derivative of.
	@param ksi: variable with respect to which one takes the derivative to.
	@n: polytropic index
	@a: dummy variable required by odeint, set to 1.
	'''

	theta, phi = y
	dydksi = [phi/ksi**2., -ksi**2. * a * theta**n]

	return dydksi

def n0_anal(x):
	'''
	Relation between theta and xi (analytical solution, Benacquista p. 113)
	'''
	return (1 - (1. / 6) * x**2)

def n1_anal(x):
	'''
	Relation between theta and xi (analytical solution)
	'''
	return (np.sin(x)/x)


# Define constants
Msun = 1.9*10**33 						# gram
Mns = 1.4 * Msun						# gram
radius = 1000000						# cm
c = 2.99792458*10**10 					# cm
hbar = 1.0546*10**-27 					# erg*s
mp = 1.6726*10**-24 					# g
G = 6.6743*10**8 						#cm3 g-1 s-2
mueWD = 2

##Ks = ......

Kwd = 3**(1./3) * np.pi**(2./3) * hbar*c / (2**(4./3) * 4 * mp**(4./3))
##a_book = (((n + 1) * Kns ) / (4 * np.pi * G))**(1./2)


# Input for ode solver
y0 = [1., 0.]
ksi = np.arange(0.001, 20, 0.001)
nvalues = [0, 1, 3./2, 2, 3, 4]
a = 1.

# Some initializations
index = 0
sol = []
theta = []
k = []
rhoPerRhoC = []

# Generate solutions for all values of n
for n in nvalues:
	sol.append([])
	sol[index] = spInt.odeint(ode, y0, ksi, args=(n, a))
	index += 1


# Purpose of the next piece of code: Only save theta values of zero and higher
# Loop over all solutions (one for every value of n)
for j in range(0, index):
	theta.append([])
	k.append([])
	rhoPerRhoC.append([])
	
	# Loop over all datapoints per value of n
	for i in range(0, len(sol[j])):
		if sol[j][i][0] < 0:
			break

		# Only save positive values of theta (and keep track of index)
		theta[j].append(sol[j][i][0])

		# Also save theta**n, since this is rho/rho_c
		rhoPerRhoC[j].append(sol[j][i][0]**nvalues[j])

		k[j] = i

n1x = np.arange(0.0001, 5, 0.0001)



# Calculate analytical value for n=0, n=1
n0 = n0_anal(ksi)
n1_dev = n1_anal(ksi)
n1 = n1_anal(n1x)


maxDev1, index1 = max(theta[0] - n0[:k[0]+1]), np.argmax(theta[0] - n0[:k[0]+1])
maxDev2, index2 = max(theta[1] - n1[:k[1]+1]), np.argmax(theta[1] - n1[:k[1]+1])



# print maxDev1, ', index = ', index1, ', theta value = ', theta[0][index1]
# print maxDev2, ', index = ', index2, ', theta value = ', theta[1][index2]
# print ''

# print 'Relative deviation from analytical value for n = 0: %f' % (maxDev1 / theta[0][index1])
# print 'Relative deviation from analytical value for n = 1: %f' % (maxDev2 / theta[1][index2])



#M_WD = 



fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(20,12))

# Plot all solutions (one for every value of n)
for num in range(index): 
	ax1.plot(ksi[:k[num]+1], theta[num], lw=2, label='n = %.2g' % nvalues[num])
	ax2.plot(ksi[:k[num]+1], rhoPerRhoC[num], lw=2, label='n = %.2g' % nvalues[num])

ax1.plot(ksi, n0, linestyle = '--', lw=2, c='black', label='Analytical solution (n=0)')
ax1.plot(n1x, n1, linestyle = '--', lw=2, c='black', label='Analytical solution (n=1)')

ax3.plot(ksi[:k[0]+1], theta[0] - n0[:k[0]+1], lw=2, label='deviation n = %.2g' % nvalues[0])
ax3.plot(ksi[:k[1]+1], theta[1] - n1_dev[:k[1]+1], lw=2, label='deviation n = %.2g' % nvalues[1])

# Plot lay-out
ax1.legend()
ax1.set_xlabel(r'$\xi$', fontsize=20)
ax1.set_ylabel(r'$\theta$($\xi$)', fontsize=20)
ax1.axis([0, 20, 0, 1])

ax2.legend()
ax2.set_xlabel(r'$\xi$', fontsize=20)
ax2.set_ylabel(r'$\rho$/$\rho_c$', fontsize=20)
ax2.axis([0, 10, 0, 1])

ax3.legend()
ax3.set_xlabel(r'$\xi$', fontsize=20)
ax3.set_ylabel('Numerical - analytical value', fontsize=20)


plt.show()