#!python3
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import lambertw

def shifter(r, M):
    # Calculate the value at which the throats are connected
    return r+ 2*M*np.log(r/(2*M)-1)

def BH_tort(rs,M):
    # Calculates the values of r given r*

    return 2*M*(1+np.real(lambertw(np.exp(rs/(2*M)-1))))

def BH_potential(rs,M, l=2):
    r0 = 2.000001*M
    shift = 0.+shifter(r0,M)
    r = BH_tort(rs+shift,M)

    print(r.min())

    return (1-2*M/r)*(l*(l+1)/r**2-6*M/r**3)

def WH_tort(rs,M):

    r0 = 2.000001*M

    shift = -shifter(r0,M)

    return np.piecewise(rs, [rs < 0, rs >= 0],
    [lambda rs: BH_tort(-rs-shift,M),
    lambda rs: BH_tort(rs-shift,M)])


def WH_potential(rs, M, l=2):

    r = WH_tort(rs, M)

    return (1-(2*M)/r)*(l*(l+1)/r**2-6*M/r**3)

def ECO_mass(r, M):
    # Piecewise mass function for the ECO case, where we have some radial mass function for rs < 0
    r0 = 2.000001*M

    return np.piecewise(r, [r < r0, r >= r0],
    [lambda r: M * (r/r0)**3,
    lambda r: M])

def ECO_inner(rs,M):
    r0 = 2.000001*M
    numer = r0**(1.5)*np.tanh(np.sqrt(2)*\
            np.sqrt(M)/r0**(1.5)*(rs+20)) # TODO weird shift
    denom = np.sqrt(2*M)

    r_inner = numer/denom

    # return np.piecewise(rs, [rs > -18, rs <= -18.],
    #                     [lambda rs: r0**(1.5)*np.tanh(np.sqrt(2)*\
    #                             np.sqrt(M)/r0**(1.5)*(rs+20))/np.sqrt(2*M),
    #                     lambda rs: 2*M)

    return r_inner



def ECO_tort(rs, M):
    # Function to convert r* to r in ECO coordinates
    r0 = 2.000001*M

    shift = -shifter(r0,M)

    rstar = np.piecewise(rs, [rs < 0, rs >= 0],
                        [lambda rs: ECO_inner(rs,M),
                        lambda rs: BH_tort(rs-shift,M)])

    return rstar

def ECO_potential(rs, M, l=2):

    r = ECO_tort(rs,M)

    r0 = 2.000001*M

    return np.piecewise(rs, [rs < 0, rs>= 0], #TODO clean up
    [lambda rs: 10*np.exp(-(rs+20)**2/4),#(1-2*M*r**2/r0**3)*(l*(l+1)/r**2-6*M/r0**3),
    lambda rs: (1-(2*M)/ECO_tort(rs,M))*(l*(l+1)/ECO_tort(rs,M)**2-6*M/ECO_tort(rs,M)**3)])

def TS_tort(rs, M):
    # Function to convert r* to r in TS coordinates
    r0 = 2.000001*M

    shift = -shifter(r0,M)

    rstar = np.piecewise(rs, [rs < 0, rs >= 0],
                        [lambda rs: rs,
                        lambda rs: BH_tort(rs-shift,M)])

    return rstar

def TS_potential(rs, M, l=2):

    r = TS_tort(rs,M)

    r0 = 2.000001*M

    return np.piecewise(r, [rs < 0, rs>= 0,r < TS_tort(-0.1,M)], #TODO clean up
    [lambda r: (l*(l+1)/r**2),
    lambda r: (1-(2*M)/r)*(l*(l+1)/r**2-6*M/r**3),
    lambda r: 10])

if __name__ == '__main__':
    M = 1
    L = 2

    xs = np.arange(-200,200.01,1e-1)*M

    print(len(xs))

    eco_tort = ECO_tort(xs,M)

    BH_pot = BH_potential(xs,M,l=L)
    WH_pot = WH_potential(xs,M,l=L)
    ECO_pot = ECO_potential(xs,M,l=L)
    TS_pot = TS_potential(xs,M,l=L)

    plt.plot(xs,ECO_pot)
    plt.show()

    potentials = np.array([xs, BH_pot, WH_pot,
                           ECO_pot, TS_pot]).T

    np.savetxt('potentials_data.dat',potentials)

    plt.plot(xs/M,BH_pot*M**2,label='BH')
    plt.plot(xs/M,WH_pot*M**2,linestyle='-.',label='WH')
    plt.plot(xs/M,ECO_pot*M**2,linestyle='--',label='ECO')
    plt.plot(xs/M,TS_pot*M**2,label='TS')
    plt.xlim(min(xs),max(xs))
    plt.ylim(-0.16,0.16)
    plt.legend()
    plt.show()
