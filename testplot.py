#!/usr/bin/env python
import numpy as np
import pylab as P

mu, sigma = 200, 25
x = mu + sigma*P.randn(10000)
bins = [100,125,150,160,170,180,190,200,210,220,230,240,250,275,300]


P.figure()
n,bins, patches = P.hist(x,50,normed=1, histtype='step' , cumulative=True)

y= P.normpdf( bins, mu, sigma).cumsum()
y/= y[-1]
l= P.plot(bins,y, 'k--', linewidth=1.5)

sigma2= 15. 
x = mu + sigma2*P.randn(1000)

n, bins , patches  = P.hist(x, bins=bins, normed=1, histtype= 'step', cumulative= True)

n,bins, patches = P.hist(x, bins=bins, normed=1, histtype='step', cumulative=-1)

P.show()



















