#! /usr/bin/env python
#
# Python dynamic status plotter
#
import sys
import glob
import os

from mpltools import style

style.use('style1')

# Defaults:
Linestyle1 = '-r'
Linestyle2 = '-b'
xLabel = r'$\sigma_{xx}$'
yLabel = r'$\sigma_{yy}$'

import numpy as np
import matplotlib
#matplotlib.use('GTKAgg')

import matplotlib.pyplot as plt
#plt.rc('font',family='serif')
#plt.rc('font',size=16)


#fig = plt.figure(1,(6,6))
fig = plt.figure(1,(3.2,3.2))
ax = fig.add_axes([0.17,0.1,0.78,0.8])
#ax.set_ylim(0, 1)
#plot 0 lines
ax.plot([-1000,1000],[0,0],'-k')
ax.plot([0,0],[-1000,1000],'-k')

ax.set_ylabel(yLabel)
ax.set_xlabel(xLabel)
ax.set_xlim((-250,450))
ax.set_ylim((-250,450))
ax.grid(True)
ax.set_aspect('equal')

fig.savefig('grid.png', transparent=True)
plt.show()


