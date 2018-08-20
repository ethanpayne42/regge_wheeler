import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl

import sys

f=open('output.dat')
lines=f.readlines()

data = []
for i, line in enumerate(lines):

    data_l = []

    line = line.replace(" ",'').strip('(').strip(')\n')

    line = line.split(')(')

    for entry in line:
        entry = entry.split(',')
        data_l.append(float(entry[0])+1j*float(entry[1]))

    if i%1000 == 0:
        sys.stdout.write('\r Up to entry: {} of {}'.format(i,len(lines)-1))
        sys.stdout.flush()

    data.append(data_l)

# Now data contains all the wave information
f.close()

f2 = open('coord.dat')
lines = f2.readlines()

xs = []
ts = []

for i,line in enumerate(lines):
    cols = line.split()
    if i == 0:
        xs = np.array(cols,dtype=float)
    elif i == 1:
        ts = np.array(cols,dtype=float)

wave = np.array(data)

#plt.plot(data.T[:,0]-data.T[:,10])
plt.plot(xs, np.real(wave[0,:]))
plt.plot(xs, np.real(wave[3900,:]))
plt.ylim(-1.1,1.1)
plt.xlim(-6,6)
plt.show()
