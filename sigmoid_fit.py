from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.mlab as mlab
import pylab as P
import pdb
import time






def sigmoid(x,b0,b1):
    return 1/( 1 + np.exp(-( b0 + b1*x ) ))
