#!/usr/bin/env python

import fileinput
import numpy as np  
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pylab as P
import pdb
import time

#   HARD CODED STUFF FOR NOW
#   mass in grams
mass_det1 = 22.748 
mass_det2 = 34.076



def plotter(fromdat):
    
    print(('\n\n  hits from unknown detector(s): {0} \n'  '   hits from {1}: {2}\n'  '  hits from {3}: {4}\n\n').format(
        fromdat.name_er, fromdat.detector1, fromdat.d1counter, fromdat.detector2, fromdat.d2counter))
    
    plt.figure() 
    bins =  np.linspace(0, 3000, 0.001*len(fromdat.all_val ) ) 
    plt.hist(fromdat.all_val, bins=bins, color=(0, 0, 0, 1 ),
                 histtype='step',label = 'All Hits' )
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title('All Detectors\n',filename )
    plt.legend(loc='upper right' ) 
    plt.show() 

    plt.figure() 
    his_det1 = plt.hist(fromdat.det1_val, bins=bins, color=(0, 0, 0, 0.7),
                 histtype='step', label = fromdat.detector1 )
   
    his_det2 = plt.hist(fromdat.det2_val, bins=bins, color=(0, 1, 0, 0.7 ),
                 histtype='step', label = fromdat.detector2 )
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title(filename)
    plt.legend(loc='upper right' ) 
    plt.show()
   
    if(fromdat.name_er != 0 ):
        his_err = plt.hist(fromdat.err_val, bins=bins, color=(0, 0, 0, 0.5 ),
                 histtype='step',label = 'Other Detectors' )
        plt.ylabel('Counts' )
        plt.xlabel('Energy kev' )
        plt.title('Other Detectors',filename )
        plt.legend(loc='upper right' ) 
        plt.show() 
        

class Sortdat(object):
    
    def __init__( self, filename, det1 , det2, thresh1 =50, thresh2 = 50 ): 
        self.d1counter = 0
        self.d2counter = 0
        self.name_er = 0
        self.all_val = []
        self.detector1 = det1
        self.detector2 = det2
        self.det1_val = []
        self.det2_val= []
        self.err_val = []
        self.thresh1 = 0 # thesholds in kev
        self.thresh2 = 0
        data = self.genspecdat(filename)        
        self.fill_class(data) 
 
    def fill_class( self, data ):
        for n, col in enumerate(data):#x is line section
             
            self.all_val.append(col[7])
           
            if (col[2]== self.detector1.encode('utf-8')):
               
                self.det1_val.append(col[7])
                self.d1counter += 1
                
            elif (col[2] == self.detector2.encode('utf-8')):
                self.det2_val.append(col[7])
                self.d2counter += 1
                
            else: 
                self.name_er += 1
                self.err_val.append(col[7])

    def gethist(self, bins ):
        plt.ioff() 
        m_fig = plt.figure() 
       
        self.nbins =  np.linspace(0, 3000, bins  ) 
        self.det1_n, self.det1_e, _ = plt.hist(self.det1_val, bins=self.nbins )
        self.det2_n, self.det2_e, _ = plt.hist(self.det2_val, bins=self.nbins )
        self.all_n, self.all_e, _  = plt.hist(self.all_val, bins=self.nbins )
        self.err_n, self.err_e, _ = plt.hist(self.err_val, bins=self.nbins )
        
        plt.close(m_fig)
        
    
    def formattxt(self, filename):
        with open(filename, 'r') as fin:
           
            combine_hit = []
            
            for line in fin:
                if line.startswith('1' ):
                    combine_hit.append(line.encode('utf-8' ) )
                else:
        
                    columns = line.split()
                    row1 = columns[:9]
                    row2 = ['1','1' ] + columns[9:]
                    row1 = '\t'.join(row1 )
                    row2 = '\t'.join(row2 )
                    combine_hit.append(row1.encode('utf-8' ) )
                    combine_hit.append(row2.encode('utf-8' ) )
        
        return combine_hit
    
    def genspecdat(self,filename):
        fromtxt = self.formattxt(filename )
        data = np.genfromtxt(fromtxt, dtype=None )
        pdb.set_trace()
        return(data)


def myrate(ratedat,filename, sim_time, cut):
    bins = 0.01*len(ratedat.all_val )
    ratedat.gethist(bins)
    
    bin_w = ratedat.nbins[2] #width for 3k bins    
    Ndet1 = []
    Ndet2 = []
    xaxis =  np.arange(cut,1500,10 ) 
    
    while (cut < 1500):

        
        bin_num = np.searchsorted(ratedat.det1_e, cut )
        partial_spec = ratedat.det1_n[bin_num : ]
        new_count = np.sum(partial_spec )
        Ndet1.append(new_count)

        bin_num = np.searchsorted( ratedat.det2_e, cut)
        partial_spec = ratedat.det2_n[bin_num : ]
        new_count = np.sum(partial_spec )
        Ndet2.append(new_count )

        cut += 10
   
    plt.figure() 
    det1plot = plt.plot(xaxis, Ndet1,label=str(sys.argv[2] ))
    det2plot = plt.plot(xaxis, Ndet2,label=str(sys.argv[3] ))
    plt.legend()
    plt.ylabel('Counts' )
    plt.title(filename)
    plt.xlabel('Low Threshold in KeV' )
    plt.draw()
    plt.show()


def main(filename):         
    sim_time = input('time of background simulation in seconds: ')
    if(not sim_time):  #These are the simulation times
        sim_time = 20 
    else:
        sim_time = int(sim_time)
    cut = input('where to start the cut (Kev): ')
    if(not cut):  #These are the threshold value
        cut = 33 
    else:
        cut = int(cut)
    
    start_time = time.time() 
    ratedat = Sortdat(filename, str(sys.argv[2] ), str(sys.argv[3] ))
    print('class creaton time' ,(time.time()-start_time),'sec')
    plotter(ratedat) 
    myrate(ratedat,filename,sim_time,cut)
#
if __name__=='__main__':

    filename = str(sys.argv[1] )
    main(filename )
    sys.exit(0 )

