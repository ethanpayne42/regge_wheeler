#!python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import lambertw

def BH_tort(rs,M):
    # Calculates the values of r given r*

    return 2*M*(1+lambertw(np.exp(rs/(2*M)-1)))

def BH_potential(rs,M, l=2):
    r = BH_tort(rs,M)

    return (1-(2*M)/r)*(l*(l+1)/r**2-6*M/r**3)

def WH_r0(r, M):
    # Calculate the value at which the throats are connected
    return r+ 2*M*np.log(r/(2*M)-1)

def WH_tort(rs,M):

    shift = -WH_r0(2.001,M)

    return np.piecewise(rs, [rs < 0, rs >= 0],
    [lambda rs: BH_tort(-rs-shift,M),
    lambda rs: BH_tort(rs-shift,M)])


def WH_potential(rs, M, l=2):

    r = WH_tort(rs, M)

    return (1-(2*M)/r)*(l*(l+1)/r**2-6*M/r**3)

def ECO_mass(rs, r, M):
    # Piecewise mass function for the ECO case, where we have some radial mass function for rs < 0
    r0 = -2.001*M

    return np.piecewise(rs, [rs < 0, rs >= 0],
    [M * (r/r0)**3,M])

def ECO_inner_tort(rs,M):
    # Calculates value of r given r* for the inner ECO region (centrifugal barrier)
    r0 = 2.001 * M
    print(r0)

    return (r0**(3/2) * np.tanh( np.sqrt(2) * np.sqrt(M/(r0**3)) * rs) / np.sqrt(2 * M))

def ECO_tort(rs, M):
    # Piecewise tortoise coordinate for ECO

    shift = -WH_r0(2.001,M)

    return np.piecewise(rs, [rs < -shift, ((-shift < rs) & (rs < 0)), rs >= 0],
    [lambda rs: 0,
    lambda rs: ECO_inner_tort(-rs+shift,M),
    lambda rs: BH_tort(rs-shift,M)])

def ECO_potential(rs, M, l=2):
    # Potential for ECO case

    r = ECO_tort(rs, M)
    r0 = 2.001 * M
    m = ECO_mass(rs,r[0:500],M)

    return (1-(2*M)/r)*(l*(l+1)/r**2-6*m/r**3)

    # if rs < 0:
    #     return (((l * (l + 1) / r**2) - (6 * M / r0**3)) * (1 - 2 * M * r**2 / r0**3))
    # elif rs >=0:
    #     return (1-(2*M)/r)*(l*(l+1)/r**2-6*M/r**3)


M = 1

xs = np.linspace(-50,50,1000)*M

plt.plot(xs, WH_tort(xs,M))
plt.plot(xs, BH_tort(xs,M))
plt.plot(xs, ECO_inner_tort(xs,M))
plt.show()

BH_pot = ECO_potential(xs,M)

print(float(ECO_potential(-10.,M)))

plt.plot(xs/M,BH_pot*M**2)
plt.plot(xs/M,0.15*np.exp(-(xs/M-10/M)**2/(6**2)))
plt.xlim(-50,50)
# plt.ylim(0,0.16)
plt.show()

potPlot = ECO_potential(xs,M)

plt.plot(xs/M, potPlot*M**2)
plt.xlim(-50,50)
plt.ylim(0,0.16)
plt.show()