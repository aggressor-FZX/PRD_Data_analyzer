#!/usr/bin/env python

"""prd_chanparse

Usage:
    prd_chanparse.py B-round_fileINPUT SIM-TIME [OUTPUT]

Options:
    INPUT   Name of input file
    TIME    Simulated time
    OUTPUT  Result of background_analysis, default: INPUT_out
"""
from docopt import docopt
import prd_auxfunctions
import logging
import textwrap  
import fileinput
import numpy as np  
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pylab as P
import pdb
from datetime import datetime
# POST MAY20 updates#
# a rad is .01 joules per kilogram
# Polimaster with several channels
# STE has no upper bound energy
# Polimaster 33-3000 kev
# STE 45 - ? kev
#-----   HARD CODED STUFF FOR NOW -------#

_log = logging.getLogger(__name__)
class Sortdat(object):
    
    detector1 = 'ste'
    detector2 = 'Polimaster'
    detector3 = 'minirad'
    det1_cut = 45
    det3_cut = 30
    PM_chan = [33,120,240,745]#from plimaster .inc email
    high_cut=3000 
    def __init__( self, filename,time): 
        self.time = time
        self.d1counter = 0
        self.d3counter = 0
        self.chan1_counter = 0
        self.chan2_counter = 0
        self.chan3_counter = 0
        self.chan4_counter = 0
        self.name_er = 0
        self.all_val = []
        self.det1_val = []
        self.det2_val= []
        self.det3_val= []
        self.det1_esum = 0 
        self.det2_esum = 0
        self.det3_esum = 0
        data = self.genspecdat(filename)        
        self.fill_class(data) 
        self.bins = 3001
        self.gethist(self.bins,self.high_cut)  
        

    def fill_class( self, data ):
        for n, col in enumerate(data):
            #col[7] is energy
            self.all_val.append(col[7])
           
            if (col[2]== self.detector1.encode('utf-8')):
                self.det1_val.append(col[7]) #for plot only
                if (col[7] >= self.det1_cut and col[7]<=self.high_cut):
                    self.det1_esum += col[7]
                    self.d1counter += 1
                    
            elif (col[2] == self.detector2.encode('utf-8')):
                self.det2_val.append(col[7]) #for plot
                #energy sum
                if (col[7] >= self.PM_chan[0]):
                     self.det2_esum += col[7]
                     self.chan4_counter+= 1
                #count channels
                if (col[7] >= self.PM_chan[0] and col[7] <= self.PM_chan[1]):
                    self.chan1_counter+= 1
                elif(col[7] <= self.PM_chan[2]):
                    self.chan2_counter+= 1
                elif(col[7] <= self.PM_chan[3]):
                    self.chan3_counter+= 1

            if (col[2]== self.detector3.encode('utf-8')):
                self.det3_val.append(col[7]) #for plot only
                if (col[7] >= self.det3_cut and col[7]<=self.high_cut):
                    self.det3_esum += col[7]
                    self.d3counter += 1

        self.all_hits = len(self.all_val)

    #find the uncertianty of the energy square(sum [e^2.n]) 
    def gethist(self, bins, high_cut):
        nbins =  np.linspace( 0, high_cut, bins ) 
        self.det1_n, self.det1_e = np.histogram(self.det1_val, bins=nbins)
        self.det2_n, self.det2_e = np.histogram(self.det2_val, bins=nbins)
        self.det3_n, self.det3_e = np.histogram(self.det3_val, bins=nbins)
        
        self.sig_e1 = np.sqrt(np.dot( self.det1_e[1:]**2, self.det1_n))
        self.sig_e2 = np.sqrt(np.dot( self.det2_e[1:]**2, self.det2_n))
        self.sig_e3 = np.sqrt(np.dot( self.det3_e[1:]**2, self.det1_n))

        self.e1_total =  np.dot( self.det1_e[1:], self.det1_n)
        self.e2_total =  np.dot( self.det2_e[1:], self.det2_n)
        self.e3_total =  np.dot( self.det3_e[1:], self.det2_n)
    
    def formattxt(self, filename):
        with open(filename, 'r') as fin:
           
            combine_hit = []
            
            for line in fin:
                if line.startswith('1' ):
                    combine_hit.append(line.encode('utf-8' ) )
                else:
        
                    columns = line.split()
                    row1 = columns[:9]
                    row1 = '\t'.join(row1 )
                    combine_hit.append(row1.encode('utf-8' ) )

                    row2 = ['1','1' ] + columns[9:18]
                    row2 = '\t'.join(row2 )
                    combine_hit.append(row2.encode('utf-8' ) )

                    if line.startswith('3'):
                        row3 = ['1','1' ] + columns[18:]
                        row3 = '\t'.join(row2 )
                        combine_hit.append(row3.encode('utf-8' ) )

                    if int(line.plit()[0]) > 3:
                       print("holly crap! went through all detectors and then some!!!!!!!!!!!!!!")
                       print("check input .dat file. If its good, remove line 146 from this code")
                       sys.exit(0 )

        return combine_hit
    
    def genspecdat(self,filename):
        fromtxt = self.formattxt(filename )
        data = np.genfromtxt(fromtxt, dtype=None )
        return(data)

#------------------Main Funcion----------------#
def main(filename, sim_time, output=None): 
    start_time = datetime.now() 
    ratedat = Sortdat(filename,sim_time) 
    print('class object creation time %s sec', (datetime.now()-start_time))
    plotter(ratedat,filename) 
    statprint(ratedat,sim_time,filename)

#makes specs 
def plotter(fromdat,filename):
    
    plt.figure() 
    bins = fromdat.bins
    plt.hist(fromdat.all_val, bins=bins, color=(0, 0, 0, 1 ),
                 histtype='step',label = 'All Hits' )
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title('All Detectors Spectrum\n'+ filename )
    plt.legend(loc='upper right' ) 
    plt.show() 

    plt.figure() 
    his_det1 = plt.hist(fromdat.det1_val, bins=bins, color=(0, 0, 0, 0.7),
                 histtype='step', label = fromdat.detector1 )
   
    his_det2 = plt.hist(fromdat.det2_val, bins=bins, color=(0, 1, 0, 0.7 ),
                 histtype='step', label = fromdat.detector2 )
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title('Overlay Plot Both Spectrum \n ' + filename)
    plt.legend(loc='upper right' ) 
    plt.show()

    his_det3 = plt.hist(fromdat.det3_val, bins=bins, color=(0, 0, 0, 0.5 ),
             histtype='step',label = fromdat.detector3 )
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title( fromdat.detector3)
    plt.legend(loc='upper right' ) 
    plt.show() 


#Calculates does rates in micro-Rad/hour
def statprint(ratedat,time,filename): 

    so_ste = Ste(ratedat) 
    bk_ste = Bkr_ste()
    stat_ste = detstat(so_ste, bk_ste)
    
    so_minirad = Minirad(ratedat)
    bk_minirad = Bkr_minirad()
    stat_mini = detstat(so_minirad, bk_minirad)
    
    so_polimaster = Polimaster(ratedat) 
    bk_poli = Bkr_Poli()
    ps1 ps2 = polistat(so_polimaster, bk_poli)
    
    _log.info('\n\n---------------------- {0} ----------------------\n'.format(stat_ste[Name])
    for key, value in stat_ste.iteritems(): 
        if key == 'Name':
            pass
        _log.info('\n %s = %s'.format(str(key),str(val)))

    _log.info('\n\n---------------------- {0} ----------------------\n'.format(stat_min[Name])
    for key, value in stat_mini.iteritems(): 
        if key == 'Name':
            pass
        _log.info('\n %s = %s'.format(str(key),str(val)))

    _log.info('\n\n---------------------- {0} ----------------------\n'.format(ps1[Name])
    _log.info('\n\n---------------------- {0} ----------------------\n'.format(ps1[Title])
    for key, value in ps1.iteritems(): 
        if key == 'Name' or key == 'Title':
            pass
        _log.info('\n %s = %s'.format(str(key),str(val)))

    _log.info('\n\n---------------------- {0} ----------------------\n'.format(ps2[Name])
    for key, value in ps2.iteritems(): 
        if key == 'Name' or 'Title':
            pass
        _log.info('\n %s = %s'.format(str(key),str(val)))

#ooooOOOOoooooOOOOooooOOOOOoooOOOOOOoooooOOOOOOooooo

if __name__=='__main__':

    args = docopt(__doc__)
    if args['OUTPUT'] is None:
        args['OUTPUT'] = args['INPUT'] + '_out'

    sim_time = float(args['TIME'])

    logging.basicConfig(level=logging.INFO,
                        format='%(message)s',
                        filename = args['OUTPUT'])
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter('%(message)s'))
    logging.getLogger('').addHandler(console)

    main(args['INPUT'], sim_time, output=args['OUTPUT'])

    sys.exit(0 )
