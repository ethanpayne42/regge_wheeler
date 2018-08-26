#!python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import lambertw

def shifter(r, M):
    # Calculate the value at which the throats are connected
    return r+ 2*M*np.log(r/(2*M)-1)

def BH_tort(rs,M):
    # Calculates the values of r given r*

    return 2*M*(1+lambertw(np.exp(rs/(2*M)-1)))

def BH_potential(rs,M, l=2):
    r = BH_tort(rs+shifter(2.001,M),M)

    return (1-(2*M)/r)*(l*(l+1)/r**2-6*M/r**3)

def WH_tort(rs,M):

    shift = -shifter(2.001*M,M)

    return np.piecewise(rs, [rs < 0, rs >= 0],
    [lambda rs: BH_tort(-rs-shift,M),
    lambda rs: BH_tort(rs-shift,M)])


def WH_potential(rs, M, l=2):

    r = WH_tort(rs, M)

    return (1-(2*M)/r)*(l*(l+1)/r**2-6*M/r**3)

def ECO_mass(r, M):
    # Piecewise mass function for the ECO case, where we have some radial mass function for rs < 0
    r0 = 2.001*M

    return np.piecewise(r, [r < r0, r >= r0],
    [lambda r: M * (r/r0)**3,
    lambda r: M])

def ECO_inner(rs,M):
    r0 = 2.001*M
    numer = r0**(1.5)*np.tanh(np.sqrt(2)*\
            np.sqrt(M)/r0**(1.5)*(rs+20)) # TODO weird shift
    denom = np.sqrt(2*M)

    r_inner = numer/denom

    return r_inner

def ECO_tort(rs, M):
    # Function to convert r* to r in ECO coordinates
    r0 = 2.001*M

    shift = -shifter(r0,M)

    rstar = np.piecewise(rs, [rs < 0, rs >= 0],
                        [lambda rs: ECO_inner(rs,M),
                        lambda rs: BH_tort(rs-shift,M)])

    return rstar

def ECO_potential(rs, M, l=2):

    r = ECO_tort(rs,M)

    r0 = 2.001*M

    return np.piecewise(r, [rs < 0, rs>= 0,r < ECO_tort(-20,M)], #TODO clean up
    [lambda r: (1-2*M*r**2/r0**3)*(l*(l+1)/r**2-6*M/r0**3),
    lambda r: (1-(2*M)/r)*(l*(l+1)/r**2-6*M/r**3),
    lambda r: np.inf])

if __name__ == '__main__':
    M = 1

    xs = np.linspace(-50,50,1000)*M

    plt.plot(xs, WH_tort(xs,M))
    plt.plot(xs, BH_tort(xs,M))
    plt.plot(xs, ECO_tort(xs,M))
    plt.show()

    BH_pot = BH_potential(xs,M)
    WH_pot = WH_potential(xs,M)
    ECO_pot = ECO_potential(xs,M)

    plt.plot(xs/M,BH_pot*M**2,label='BH')
    plt.plot(xs/M,WH_pot*M**2,linestyle='-.',label='WH')
    plt.plot(xs/M,ECO_pot*M**2,linestyle='--',label='ECO')
    plt.xlim(-50,50)
    plt.ylim(0,0.16)
    plt.legend()
    plt.show()
