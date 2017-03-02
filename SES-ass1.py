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
import scipy.stats as spStat
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
Msun = 1.99*10**33 						# gram
Mns = 1.4 * Msun						# gram
radiusNS = 1000000						# cm
c = 2.99792458*10**10 					# cm
hbar = 1.0546*10**-27 					# erg*s
mp = 1.6726*10**-24 					# g
G = 6.6743*10**-8 						#cm3 g-1 s-2
mueWD = 2

##Ks = ......

Kwd = 3**(1./3) * np.pi**(2./3) * hbar*c / (2**(4./3) * 4 * mp**(4./3))
##a_book = (((n + 1) * Kns ) / (4 * np.pi * G))**(1./2)


# Input for ode solver
stepsize = 0.0001
y0 = [1., 0.]
ksi = np.arange(stepsize, 20, stepsize)
nvalues = [0, 1, 3./2, 2, 3, 4]
a = 1.

# Some initializations
index = 0
ksi1 = []
sol = []
theta = []
phi = []
k = []
rhoPerRhoC = []

# Generate solutions for all values of n
for n in nvalues:
	sol.append([])
	sol[index] = spInt.odeint(ode, y0, ksi, args=(n, a))
	index += 1
bla = 0

# Purpose of the next piece of code: Only save theta values of zero and higher
# Loop over all solutions (one for every value of n)

for j in range(0, index):
	theta.append([])
	phi.append([])
	k.append([])
	rhoPerRhoC.append([])
	
	# Loop over all datapoints per value of n
	for i in range(0, len(sol[j])):
		if sol[j][i][0] < 0:
			ksi1.append(i-1)
			bla = sol[j][i][1]
			break

		# Only save positive values of theta (and keep track of index)
		theta[j].append(sol[j][i][0])
		phi[j].append(sol[j][i][1])

		# Also save theta**n, since this is rho/rho_c
		rhoPerRhoC[j].append(sol[j][i][0]**nvalues[j])

		k[j] = i

n1x = np.arange(stepsize, 5, stepsize)



# Calculate analytical value for n=0, n=1
n0 = n0_anal(ksi)
n1_dev = n1_anal(ksi)
n1 = n1_anal(n1x)


maxDev1, index1 = max(theta[0] - n0[:k[0]+1]), np.argmax(theta[0] - n0[:k[0]+1])
maxDev2, index2 = max(theta[1] - n1[:k[1]+1]), np.argmax(theta[1] - n1[:k[1]+1])



print maxDev1, ', index = ', index1, ', theta value = ', theta[0][index1]
print maxDev2, ', index = ', index2, ', theta value = ', theta[1][index2]
print ''

print 'Relative deviation from analytical value for n = 0: %f' % (maxDev1 / theta[0][index1])
print 'Relative deviation from analytical value for n = 1: %f' % (maxDev2 / theta[1][index2])


ttestn0 = spStat.ttest_ind(n0[:k[0]+1], theta[0])
ttestn1 = spStat.ttest_ind(n1[:k[1]+1], theta[1])

ttestn0_rel = spStat.ttest_rel(n0[:k[0]+1], theta[0])
ttestn1_rel = spStat.ttest_rel(n1[:k[1]+1], theta[1])

# print ttestn0
# print ttestn1

# print 'rel0 = ', ttestn0_rel
# print 'rel1 = ', ttestn1_rel

	
meanDensNS = Mns / ((4. / 3) * np.pi * radiusNS**3)
rhoCentNS = meanDensNS / ((-3. * phi[1][ksi1[1]] / ksi[ksi1[1]]**3))

print '\nAverage density NS: %.4g g/cm3 \nRho center NS: %.4g g/cm3' % (meanDensNS, rhoCentNS)


mass_WD = -1 * ( 4 * phi[4][ksi1[3]] / np.sqrt(np.pi) ) * (Kwd / G)**1.5

print '\nMass WD = %.4g Mdot\n' % float(mass_WD/Msun)



# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18,12))

# # Plot all solutions (one for every value of n)
# for num in range(index): 
# 	ax1.plot(ksi[:k[num]+1], theta[num], lw=2, label='n = %.2g' % nvalues[num])
# 	ax2.plot(ksi[:k[num]+1], rhoPerRhoC[num], lw=2, label='n = %.2g' % nvalues[num])

# ax1.plot(ksi, n0, linestyle = '--', lw=2, c='black', label='Analytical solution (n=0)')
# ax1.plot(n1x, n1, linestyle = '--', lw=2, c='black', label='Analytical solution (n=1)')


# # Plot lay-out
# ax1.legend()
# ax1.set_title(r'$\theta(\xi)$ for different n', fontsize=20)
# ax1.set_xlabel(r'$\xi$', fontsize=20)
# ax1.set_ylabel(r'$\theta$($\xi$)', fontsize=20)
# ax1.axis([0, 20, 0, 1])

# ax2.legend()
# ax2.set_title(r'$\frac{\rho}{\rho_c}(\xi)$ for different n', fontsize=20)
# ax2.set_xlabel(r'$\xi$', fontsize=20)
# ax2.set_ylabel(r'$\frac{\rho}{\rho_c}$    ', fontsize=25, rotation=0)
# ax2.axis([0, 10, 0, 1])

# plt.tight_layout()
# # plt.show()



# fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(18,12))

# ax1.axhline(0, color='r', linestyle='dashed', lw=2, label='Zero deviation')
# ax1.vlines(ksi[:k[0]+1], 0, theta[0] - n0[:k[0]+1], color='blue', lw=2, alpha=0.5, label='deviation n = %.2g' % nvalues[0])
# ax1.plot(ksi[:k[0]+1], theta[0] - n0[:k[0]+1], 'bo', color='blue', lw=2, label='deviation n = %.2g' % nvalues[0])

# ax2.axhline(0, color='r', linestyle='dashed', lw=2, label='Zero deviation')
# ax2.vlines(ksi[:k[1]+1], 0, theta[1] - n1_dev[:k[1]+1], color='green', lw=2, alpha=0.5, label='deviation n = %.2g' % nvalues[1])
# ax2.plot(ksi[:k[1]+1], theta[1] - n1_dev[:k[1]+1], 'bo', color='green', lw=2, label='deviation n = %.2g' % nvalues[1])

# ax1.legend(loc='lower right')
# ax2.legend(loc='lower right')

# ax1.set_title('Deviation from analytical value for n = 0', fontsize=20)
# ax1.set_xlabel(r'$\xi$', fontsize=20)
# ax1.set_ylabel('Numerical - analytical value', fontsize=20)

# ax2.set_title('Deviation from analytical value for n = 1', fontsize=20)
# ax2.set_xlabel(r'$\xi$', fontsize=20)
# ax2.set_ylabel('Numerical - analytical value', fontsize=20)

# # if stepsize == 0.1:
# # 	ax1.axis([0, 3.5, 0, 0.005])
# # 	ax2.axis([0, 3.5, 0, 0.005])
# # elif stepsize == 0.01:
# # 		ax1.axis([0, 3.5, 0, 0.00006])	
# # 		ax2.axis([0, 3.5, 0, 0.00006])
# # elif stepsize == 0.001:
# # 		ax1.axis([0, 3.5, -0.0000018, 0.0000005])	
# # 		ax2.axis([0, 3.5, -0.0000018, 0.0000005])
# # elif stepsize == 0.001:
# # 		ax1.axis([0, 3.5, -0.0000035, 0.0000005])	
# # 		ax2.axis([0, 3.5, -0.0000035, 0.0000005])

# plt.tight_layout()
# #plt.show()


fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(18,12))

ax1.vlines(ksi[:k[0]+1], 0, (theta[0] - n0[:k[0]+1])/n0[:k[0]+1], color='blue', lw=2, alpha=0.5)
ax1.plot(ksi[:k[0]+1], (theta[0] - n0[:k[0]+1])/n0[:k[0]+1], 'bo', color='blue', lw=2, label='deviation n = %.2g' % nvalues[0])
# if stepsize == 0.1:
# 	ax1.axhline(0.01, color='r', linestyle='dashed', lw=2, label=r'$1.0\%$ limit')	
# elif stepsize == 0.01:
# 	ax1.axhline(0.0001, color='r', linestyle='dashed', lw=2, label=r'$0.01\%$ limit')
# 	ax1.axhline(0.0005, color='g', linestyle='dashed', lw=2, label=r'$0.05\%$ limit')

ax2.vlines(ksi[:k[1]+1], 0, (theta[1] - n1_dev[:k[1]+1])/n1_dev[:k[1]+1], color='green', lw=2, alpha=0.5)
ax2.plot(ksi[:k[1]+1], (theta[1] - n1_dev[:k[1]+1])/n1_dev[:k[1]+1], 'bo', color='green', lw=2, label='deviation n = %.2g' % nvalues[1])
# if stepsize == 0.1:
# 	ax2.axhline(0.004, color='r', linestyle='dashed', lw=2, label=r'$0.4\%$ limit')	
# elif stepsize == 0.01:
# 	ax2.axhline(0.0001, color='r', linestyle='dashed', lw=2, label=r'$0.01\%$ limit')
# elif stepsize == 0.001:
# 	ax2.axhline(0.0001, color='r', linestyle='dashed', lw=2, label=r'$0.01\%$ limit')
# elif stepsize == 0.0001:
# 	ax2.axhline(0.0001, color='r', linestyle='dashed', lw=2, label=r'$0.01\%$ limit')

ax1.legend(loc='upper left')
ax2.legend(loc='upper left')
ax1.grid(lw=1)
ax2.grid(lw=1)

ax1.set_title('Relative deviation from analytical value for n = 0 (stepsize = %s)' % stepsize, fontsize=20)
ax1.set_xlabel(r'$\xi$', fontsize=24)
ax1.set_ylabel(r'$\%$ deviation', fontsize=20)
ax1.tick_params(labelsize=15)

ax2.set_title('Relative deviation from analytical value for n = 1 (stepsize = %s)' % stepsize, fontsize=20)
ax2.set_xlabel(r'$\xi$', fontsize=24)
ax2.set_ylabel(r'$\%$ deviation', fontsize=20)
ax2.tick_params(labelsize=15)

plt.tight_layout()
plt.show()




#### COMPARISON WITH ANAL. SOLUTIONS
# figCOMP, (ax1) = plt.subplots(1, 1, figsize=(10,10))

# # Plot all solutions (one for every value of n)
# index = 2
# for num in range(index): 
# 	ax1.plot(ksi[:k[num]+1], theta[num], lw=3, alpha=0.4, label='n = %.2g' % nvalues[num])
# 	# ax2.plot(ksi[:k[num]+1], rhoPerRhoC[num], lw=2, label='n = %.2g' % nvalues[num])

# ax1.plot(ksi, n0, linestyle = '--', lw=2, c='darkblue', label='Analytical solution (n=0)')
# ax1.plot(n1x, n1, linestyle = '--', lw=2, c='darkgreen', label='Analytical solution (n=1)')


# # Plot lay-out
# ax1.legend()
# ax1.set_title('Comparison with analytical solutions (stepsize = %s)' % (stepsize), fontsize=20)
# ax1.set_xlabel(r'$\xi$', fontsize=20)
# ax1.set_ylabel(r'$\theta$($\xi$)', fontsize=20)
# ax1.axis([0, 3.5, 0, 1])
# ax1.tick_params(labelsize=15)

# # ax2.legend()
# # ax2.set_title(r'$\frac{\rho}{\rho_c}(\xi)$ for different n', fontsize=20)
# # ax2.set_xlabel(r'$\xi$', fontsize=20)
# # ax2.set_ylabel(r'$\frac{\rho}{\rho_c}$    ', fontsize=25, rotation=0)
# # ax2.axis([0, 10, 0, 1])

# plt.tight_layout()
# plt.show()