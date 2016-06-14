#!/usr/bin/env python
import numpy as np
import pylab as P
import numpy
from matplotlib import pyplot
import random

import fileinput

#P.figure()
 


#n,bins, patches = P.hist(x,bins, histtype='step' , cumulative=True)
#
#y= P.normpdf( bins, 13, 4).cumsum()
#y/= y[-1]
#l= P.plot(bins,y, 'k--', linewidth=1.5)
#
# 
#
#n, bins , patches  = P.hist(y, bins=bins, normed=1, histtype= 'step', cumulative= True)
#
#n,bins, patches = P.hist(x, bins=bins, normed=1, histtype='step', cumulative=-1)
#
#P.show()
np.set_printoptions(precision=4)
bins = np.arange(1, 22)
Outfile = "fakedat.txt"

def gaussian(x, mu=3.0, sigma=1.0):
    tss = 2 * np.power(sigma, 2)
    return np.exp( -np.power(x-mu, 2) / tss ) / np.sqrt( tss * np.pi )

l = np.random.randint(3,high=8)
q = np.random.randint(1,high=5) 
with open(Outfile, "w") as fout:
    for i in bins: 
        col1 = i  
        col2 = gaussian(i,q,1)
        col3 = gaussian(i,l,1)
        fout.write('{0:d}\t\t{1:.5f}\t\t{2:.5f} \n'.format(col1,col2,col3))
        


                










