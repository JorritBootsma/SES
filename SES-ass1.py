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


# Input for ode solver
y0 = [1., 0.]
ksi = np.arange(0.01, 20, 0.0001)
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


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15,10))

# Plot all solutions (one for every value of n)
for num in range(index): 
	ax1.plot(ksi[:k[num]+1], theta[num], label='n = %.2g' % nvalues[num])
	ax2.plot(ksi[:k[num]+1], rhoPerRhoC[num], label='n = %.2g' % nvalues[num])

# Plot lay-out
ax1.legend()
ax1.set_xlabel(r'$\xi$', fontsize=20)
ax1.set_ylabel(r'$\theta$($\xi$)', fontsize=20)
ax1.axis([0, 20, 0, 1])

ax2.legend()
ax2.set_xlabel(r'$\xi$', fontsize=20)
ax2.set_ylabel(r'$\rho$/$\rho_c$', fontsize=20)
ax2.axis([0, 10, 0, 1])


plt.show()