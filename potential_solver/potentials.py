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

M = 1

xs = np.linspace(-50,50,1000)*M

plt.plot(xs, WH_tort(xs,M))
plt.plot(xs, BH_tort(xs,M))
plt.show()

BH_pot = WH_potential(xs,M)

print(float(WH_potential(10.,M)))

plt.plot(xs/M,BH_pot*M**2)
plt.plot(xs/M,0.15*np.exp(-(xs/M-10/M)**2/(6**2)))
plt.xlim(-50,50)
plt.show()
