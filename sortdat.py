#!/usr/bin/env python

import fileinput
import numpy as np  
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pylab as P
import pdb
import time
def plotter(fromdat):
    
    print(('\n{0} hits from unknown detector\n\n'
           '{1} hits from  detector1 {1} \n\n'
           '{2} hits from detector2 {1} \n\n')
           .format(fromdat.name_er,
                   fromdat.d1counter,
                   fromdat.detector1,
                   fromdat.d2counter,
                   fromdat.detector2)
    )
    
    plt.figure() 
    bins =  np.linspace(0, 3000, 0.001*len(fromdat.all_val ) ) 
    plt.hist(fromdat.all_val, bins=bins, color=(1, 0, 0, 1 ),
                 histtype='step',label = 'All Hits' )
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title('All Detectors' )
    plt.legend(loc='upper right' ) 
    plt.show() 
    plt.figure() 
    plt.hist(fromdat.det1_val, bins=bins, color=(1, 0, 0, 0.7),
                 histtype='step', label = fromdat.detector1 )
   
    plt.hist(fromdat.det2_val, bins=bins, color=(0, 1, 0, 0.7 ),
                 histtype='step', label = fromdat.detector2 )
    
    if(fromdat.name_er != 0 ): 
        plt.hist(fromdat.err_val, bins=bins, color=(1, 1, 0, 0.5 ),
                 histtype='step',label = 'Other Detectors' )
    
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title(filename)
    plt.legend(loc='upper right' ) 
   
    plt.show()


class Sortdat(object):

    def formattxt(self, filename):
        with open(filename, 'r') as fin:
            #combine_hit = []
            for line in fin:
                if line.startswith('1' ):
                    first_row = line.encode('utf-8' )
                    yield first_row
                else:
        
                    columns = line.split()
                    row1 = columns[:9]
                    row2 = ['1','1' ] + columns[9:]
                    row1 = '\t'.join(row1 )
                    row2 = '\t'.join(row2 )
                    first_row = row1.encode('utf-8' )
                    second_row = row2.encode('utf-8' )
                    yield first_row 
                    yield second_row
    
    def buildspec( self, col ):
             
            self.all_val.append(col.item(0)[7])
           
            if (col.item(0)[2]== self.detector1.encode('utf-8')):
               
                self.det1_val.append(col.item(0)[7])
                self.d1counter += 1
                
            elif (col.item(0)[2] == self.detector2.encode('utf-8')):
                self.det2_val.append(col.item(0)[7])
                self.d2counter += 1
                
            else: 
                self.name_er += 1
                self.err_val.append(col.item(0)[7])

    def histdat(self, bins ):
        plt.ioff() 
        m_fig = plt.figure() 
       
        self.nbins =  np.linspace(0, 3000, bins  ) 
        self.det1_n, self.det1_e, _ = plt.hist(self.det1_val, bins=self.nbins )
        self.det2_n, self.det2_e, _ = plt.hist(self.det2_val, bins=self.nbins )
        self.all_n, self.all_e, _  = plt.hist(self.all_val, bins=self.nbins )
        self.err_n, self.err_e, _ = plt.hist(self.err_val, bins=self.nbins )
        pdb.set_trace() 
        
        plt.close(m_fig)

        #return combine_hit
    def __init__(self, filename, det1 , det2 ): 
        start_time = time.time() 
        self.d1counter = 0
        self.d2counter = 0
        self.name_er = 0
        self.all_val = []
        self.detector1 = det1
        self.detector2 = det2
        self.det1_val = []
        self.det2_val= []
        self.err_val = []
        self.row_gen = self.formattxt(filename)
        try:
            while True:
                self.row = [next(self.row_gen)]
                self.rowdat = np.genfromtxt(self.row, dtype=None,delimiter="\t" )
                self.buildspec(self.rowdat )
        except StopIteration:
            pass
        print( 'class create time ', ( time.time()-start_time),' sec') 

def myrate(filename, time, cut):
    
    ratedat = Sortdat(filename, str(sys.argv[2] ), str(sys.argv[3] ))
    bins = 0.01*len(ratedat.all_val )
    ratedat.histdat(bins )
    
    bin_w = ratedat.nbins[2]    
    Ndet1 = []
    Ndet2 = []
    xaxis = np.arange(cut,1500,10 ) 
   
    while (cut < 1500 ):

        bin_num = np.searchsorted(ratedat.det1_e, cut )
        partial_spec = ratedat.det1_n[bin_num : ]
        new_count = np.sum(partial_spec )
        Ndet1.append(new_count )

        bin_num = np.searchsorted(ratedat.det2_e, cut )
        partial_spec = ratedat.det2_n[bin_num : ]
        new_count = np.sum(partial_spec )
        Ndet2.append(new_count )

        cut += 10
   
    plt.figure() 
    det1plot = plt.plot(xaxis, Ndet1,label=str(sys.argv[2] ))
    det2plot = plt.plot(xaxis, Ndet2,label=str(sys.argv[3] ))
    plt.legend()
    plt.ylabel('Counts' )
    plt.xlabel('Low Threshold in KeV' )
    plt.show()


def Eplots(filename):    

    data = genspecdat(filename)
    fromdat = Sortdat(str(sys.argv[2] ), str(sys.argv[3] ))

    fromdat.buildspec(data )
    plotter(fromdat ) 



def main(filename):         
    timeb = input('time of background simulation in seconds: ')
    if(not timeb):  #These are the simulation times
        timeb = 3600 
    cut = input('where to start the cut: ')
    if(not cut):  #These are the threshold value
        cut = 50 
    answer = input( 'plot? y/n ')
    if(answer == 'y'):
        Eplots(filename)

    myrate(filename, timeb, cut )

if __name__=='__main__':

    filename = str(sys.argv[1] )
    main(filename )
    sys.exit(0 )

