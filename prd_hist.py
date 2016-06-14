#!/usr/bin/env python
"""prd_chanparse

Usage:
    prd_chanparse.py INPUT TIME [OUTPUT]

Options:
    INPUT   Name of input file
    TIME    Simulated time
    OUTPUT  Result of analysis, default: INPUT_out
"""
from docopt import docopt

import logging
import textwrap  
import fileinput
import numpy as np  
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pylab as P
import pdb
import time
# POST MAY20 updates#
# a rad is .01 joules per kilogram

# Polimaster with several channels
# STE has no upper bound energy
# Polimaster 33-3000 kev
# STE 45 - ? kev
#-----   HARD CODED STUFF FOR NOW -------#
mass_det1 = 22.748 #mass of STE crystal
mass_det2 = 34.076 #mass of Polimaster crystal
mass_det3 = 14.28  #mass of minirad crystal 
mass_d1kg = mass_det1/1000
mass_d2kg = mass_det2/1000
mass_d3kg = mass_det3/1000

_log = logging.getLogger(__name__)
class Sortdat(object):
    
    detector1 = 'ste'
    detector2 = 'Polimaster'
    detector3 = 'minirad'
    det1_cut = 45
    det3_cut = 30
    PM_chan = [33,120,240,745]
    high_cut=3000 
    def __init__( self, filename): 

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
                        
        return combine_hit
    
    def genspecdat(self,filename):
        fromtxt = self.formattxt(filename )
        data = np.genfromtxt(fromtxt, dtype=None )
        return(data)

#------------------Main Funcion----------------#
def main(filename, sim_time, output=None): 
    start_time = time.time() 
    ratedat = Sortdat(filename) 
    print('class object creation time %f sec', (time.time()-start_time))
    plotter(ratedat,filename) 
    doserate(ratedat,sim_time,filename)

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
def doserate(ratedat,time,filename): 
    _log.info('\n---------------------- %s ----------------------', filename)
    _log.info('\tTime = \t{0}\t seconds\n'.format(time))  
    conv = 1.60218*10**-16#conversion factor kev to 
    
    #STE doserate
    micro_rad_per_hr = ( (100*ratedat.det1_esum*conv) / (mass_d1kg*(time/3600)) )*(10**6)    
    sig_m_rad_1 = ( (100*ratedat.sig_e1*conv) / (mass_d1kg*(time/3600)) )*(10**6)    

    #STE results
    _log.info(textwrap.dedent("""\
    {0}\t cut =\t {1}\tkev 
    \thit rate:\t {2}\t counts per second
    \tmicro rad per hour:\t {3}
    \tuncertainty in (u)rads/hour :\t{4}\tmicro rads per hour
    """)
    .format(ratedat.detector1,ratedat.det1_cut,ratedat.d1counter/time, micro_rad_per_hr, sig_m_rad_1))
    print("the value total energy one %f" % ratedat.det1_esum)
    print("the value of uncertianty for energy one %f" % ratedat.sig_e1)

    #Polimaster doserate
    micro_rad_per_hr = ( (100*ratedat.det2_esum*conv) / (mass_d2kg*(time/3600)) )*(10**6)    
    sig_m_rad_2 = ( (100*ratedat.sig_e2*conv) / (mass_d2kg*(time/3600)) )*(10**6)    
    
    #Polimaster results   
    _log.info(textwrap.dedent("""\
    {0}\t cut =\t {1}\tkev 
    \tmicro rad per hour:\t {2}
    \tuncertainty in (u)rads/hour :\t{3}\tmicro rads per hour
    """)
    .format(ratedat.detector2,ratedat.PM_chan[0], micro_rad_per_hr, sig_m_rad_2))
    print("the value of uncertianty for energy two %f" % ratedat.sig_e2)
    print("the value total energy two %f" % ratedat.det2_esum)

    _log.info('\tchannel 1 rate {0} to {1} Kev:\t {2}\t counts per second'.
        format(ratedat.PM_chan[0],ratedat.PM_chan[1],ratedat.chan1_counter/time))

    _log.info('\tchannel 2 rate {0} to {1} Kev:\t {2}\t counts per second'.
        format(ratedat.PM_chan[1],ratedat.PM_chan[2],ratedat.chan2_counter/time))

    _log.info('\tchannel 3 rate {0} to {1} Kev:\t {2}\t counts per second'.
        format(ratedat.PM_chan[2],ratedat.PM_chan[3],ratedat.chan3_counter/time))

    _log.info('\tchannel 4 all hits in PM range {0} to {1} Kev:\t {2}\t counts per second'.
        format(ratedat.PM_chan[0],ratedat.PM_chan[3],ratedat.chan4_counter/time))


    #minirad doserate
    micro_rad_per_hr = ( (100*ratedat.det3_esum*conv) / (mass_d3kg*(time/3600)) )*(10**6)    
    sig_m_rad_3 = ( (100*ratedat.sig_e1*conv) / (mass_d3kg*(time/3600)) )*(10**6)    

    #minirad results
    _log.info(textwrap.dedent("""\
    {0}\t cut =\t {1}\tkev 
    \thit rate:\t {2}\t counts per second
    \tmicro rad per hour:\t {3}
    \tuncertainty in (u)rads/hour :\t{4}\tmicro rads per hour
    """)
    .format(ratedat.detector3,ratedat.det3_cut,ratedat.d3counter/time, micro_rad_per_hr, sig_m_rad_3))
    print("the value total energy one %f" % ratedat.det3_esum)
    print("the value of uncertianty for energy one %f" % ratedat.sig_e3)

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
