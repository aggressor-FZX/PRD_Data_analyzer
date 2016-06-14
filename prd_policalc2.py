#!/usr/bin/env python

"""prd_chanparse

Usage:
    prd_policalc.py INPUT TIME 

Options:
    INPUT   Name of input file
    TIME    Simulated time
"""
from docopt import docopt
import time
import prd_auxfunctions2 as aux
import logging
import textwrap  
import fileinput
import numpy as np  
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pylab as P
import pdb
from os import path
from datetime import datetime
import prd_allplots as pap
# POST MAY20 updates#
# a rad is .01 joules per kilogram
# Polimaster with several channels
# STE has no upper bound energy
# Polimaster 33-3000 kev
# STE 45 - ? kev
#-----   HARD CODED STUFF FOR NOW -------#
config = 'myconfig.pyth'

class Sortdat(object):   

    detector1 = 'ste'
    detector2 = 'Polimaster'
    detector3 = 'minirad'
    det1_cut = 45
    det3_cut = 30
    PM_chan = [33,120,240,745]#from plimaster .inc email
    high_cut = 3000 

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
        data , remainders = self.genspecdat(filename)        
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
            weirdos = []
            for x, line in enumerate(fin):
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
                    
                    if int(line.split()[0]) > 3:
                        weirdos.append(line)
                        print("holly crap! went through all detectors and then some!!!!!!!!!!!!!!")
                        print('line number {0}'.format(x))
                        print(line)
                    

                        #sys.exit(0 )

        return combine_hit, weirdos
    
    def genspecdat(self,filename):
        fromtxt , remiander = self.formattxt(filename )
        data = np.genfromtxt(fromtxt, dtype=None )
        return data, remiander


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
def statprint(ratedat,sim_time,filename,sys_time,date,dist,allopt, all_poli, all_mini): 

    so_ste = aux.Ste(ratedat) 
    bk_ste = aux.Bkr_ste(ratedat)
    stat_ste = aux.stestat(so_ste, bk_ste)
    
    so_minirad = aux.Minirad(ratedat)
    bk_minirad = aux.Bkr_minirad(ratedat)
    stat_mini = aux.ministat(so_minirad, bk_minirad)
    
    so_polimaster = aux.Polimaster(ratedat) 
    bk_poli = aux.Bkr_poli(ratedat)
    ps1, ps2 = aux.polistat(so_polimaster, bk_poli)
    
    if allopt:    
        file_ste= open('subreport-ste','a')
        file_minirad = open('subreport-minirad','a')
        file_poli = open('subreport-polimaster','a')
        files = [file_ste,file_minirad,file_poli] 
        #import pdb; pdb.set_trace() 
        all_poli.fill(ps1,ps2,dist,sim_time)
        all_mini.fill(stat_mini,dist,sim_time) 
    else:
        file_ste= open('report-ste','a')
        file_minirad = open('report-minirad','a')
        file_poli = open('report-polimaster','a')
        files = [file_ste,file_minirad,file_poli] 
#Writes out to each prd file
    for f in files: 
        f.write('\n \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ')
        f.write('\n\n\---------------------- {0} ----------------------'.format(filename))
        f.write('\n Date:\t{0}/{1}/{2}\n time:\t{3}'.format(date.day,date.month,date.year,sys_time))
        f.write('\n\--------------Source Sim Time\t{0} ----------'.format(sim_time))
    
    for key, value in stat_ste.items(): 
        file_ste.write('\n{0}\t{1}'.format(str(key),str(value)))
    file_ste.close()

    for key, value in stat_mini.items(): 
        file_minirad.write('\n{0}\t{1}'.format(str(key),str(value)))
    file_minirad.close()

    for key, value in ps1.items(): 
        file_poli.write('\n{0}\t{1}'.format(str(key),str(value)))
    for key, value in ps2.items(): 
        file_poli.write('\n{0}\t{1}'.format(str(key),str(value)))
    file_poli.close()

    return files

def allout(config):
    allopt = True
    all_poli = aux.Polistore()
    all_minirad = aux.Ministore()
    with open(config, 'r') as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith('#'):
               continue 
            line = line.split()
            path = line[0]
            sim_time = float(line[1])
            distance = line[2]
            main( path, sim_time, distance , allopt, all_poli, all_minirad )
    return all_poli, all_minirad 

#------------------Main Funcion----------------#

def main(filename, sim_time, dist = None, allopt = False, all_poli = None, all_minirad = None): 
    start_time = datetime.now() 
    date = datetime.now()
    sys_time = repr(time.strftime("%H:%M:%S"))
    ratedat = Sortdat(filename,sim_time) 
    # plotter(ratedat,filename) 
    files = statprint(ratedat,sim_time,filename,sys_time,date,dist,allopt,all_poli,all_minirad)
    return files 


#ooooOOOOoooooOOOOooooOOOOOoooOOOOOOoooooOOOOOOooooo
#[0] script name, [input_file], [Time]

if __name__=='__main__':
    print('use argv "runconfig" to load myconfig.pyth\n')
    if sys.argv[1] == 'runconfig' : 
        all_poli, all_minirad = allout(config)
        pap.plotpoli(all_poli)
        #pap.plotmini(all_mini)
    else:       
        filename = sys.argv[1]
        sim_time = float(sys.argv[2])
        files = main(filename, sim_time)
        for name in files:
            print('writing to file: '.format(name))
    sys.exit(0)



