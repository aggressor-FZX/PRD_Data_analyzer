
#!/usr/bin/env python
# a rad is .01 joules per kilogram
# Polimaster with several channels
# STE has no upper bound energy

# POST MAY20 updates
#----------------------Polimaster [India]Notes----------------------#
# Polimaster 33-3000 kev
#The analysis is being done within two time intervals 1s and 2s . 
#During calibration, the average count rate for 1 and 2 seconds are being calculated for every energy window 
#and for a combination of several windows combined. 
#
#The sigmas are equal to square root of a count rate.
#8 quarter second intervals are being stored for every energy window at any time. 
#For 1s analysis slot, the most recent 4 quarter-second intervals are being analyzed. 
#----------------------------------Minirad [Tango]----------------------#
# The background is calculated for 30 seconds 
# this is divided by 30 to get the rate
# Limited to 30 CPSec
# backgroun is mapped to a 4sigma table
#A Hex Switch setting mappes to a hex swtich adjust value table
#Formula : base alarm-count = background +4sigma(background)+(Hexswtich adjustvalue)
# unit alarms when 2 sec av goes above count rate base count rate: alarm level is determined by a table
# uses calibration to establish a threshold
# recalibrates itself somehow 
# Gamma range goes from 30 Kev 3MeV 
#----------------------STE [Sierra] Notes----------------------#
#--------From the Manual & design spec sheet:
# Low energy cut is 45KeV
from scipy.stats import norm
from scipy.stats import poisson 
from scipy import interpolate
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np  
from collections import OrderedDict
from scipy.optimize import curve_fit
import pylab

conv = 1.60218*10**-16 #conversion factor kev to 
dist = [451,488,525,563] #distances
class Bkr_ste(object):


    rate = 15.14
    itime = 1

    def __init__(self,sortdat):
        self.counts =  self.rate*self.itime
        self.sig_counts = np.sqrt(self.counts)  
        self.engeryrate = 4.92
        self.sig_energyrate = .05
        self.source_time = sortdat.time

class Bkr_minirad(object):
    #values are hard coded based on background run
    rate = 8.915
    alarm = 5.4
    itime = 2
    hex_set = 0
    #from the spec sheet
    four_sigma = {x: np.around(np.sqrt(x)*4) for x in np.arange(1,31)}
    hex_adjust = {x: x*2 for x in np.arange(6)}

    def __init__(self,sortdat):
        self.counts = self.rate*self.itime
        self.sig_counts = np.sqrt(self.counts)  
        self.foursig = self.four_sigma[np.around(self.rate)]
        self.hexval = self.hex_adjust[np.around(self.hex_set)]
        self.base_alarm = self.rate + self.hexval + self.foursig
        self.source_time = sortdat.time
        
# make detector class containers
class Polimaster(object):

    codename = 'India'
    alarm = 9.8   # sigma above bkrd alarm setting
    mass = 34.076 # mass of Polimaster crystal
    itime1 = 1
    itime2 = 2
    mass_kg = mass/1000
    actualprob =    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.8, .63, .03, 0.0, 0.0, 0.0, 0.0, 0.0]
    actualdist =    [282, 316, 349, 372, 406, 451, 488, 525, 563, 601, 638, 676, 789, 902]
    actualsig_pos = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .10, .13, .11, .09, .09, .09, .09, .09]
    actualsig_neg = [.09, .09, .09, .09, .09, .09, .15, .15, .02, .00, .00, .00, .00, .00]
        
    def __init__(self,sortdat):

        self.s_time = sortdat.time
        self.name = sortdat.detector2
        self.ch1_cut = sortdat.PM_chan[0]
        self.ch2_cut = sortdat.PM_chan[1]
        self.ch3_cut = sortdat.PM_chan[2]
        self.ch4_cut = sortdat.PM_chan[3]

        self.ch1_counts = sortdat.chan1_counter
        self.ch2_counts = sortdat.chan2_counter
        self.ch3_counts = sortdat.chan3_counter
        self.ch4_counts = sortdat.chan4_counter

        self.ch1_rate = self.ch1_counts/sortdat.time
        self.ch2_rate = self.ch2_counts/sortdat.time
        self.ch3_rate = self.ch3_counts/sortdat.time
        self.ch4_rate = self.ch4_counts/sortdat.time

        self.energy = sortdat.det2_val
        self.t_energy = sortdat.det2_esum
        self.sig_energy = sortdat.sig_e2
        self.high_cut = sortdat.PM_chan[3]
        self.source_time = sortdat.time

class Ste(object):
    codename = 'Sierra'
    alarm = 9.8    # sigma above bkrd alarm setting
    mass = 22.748  # mass of STE crystal
    itime = 1      # integration time guess
    high_cut = 3000 
    mass_kg = mass/1000
    actualdist =    [282, 316, 349, 372, 406, 451, 488, 525, 563, 601, 638, 676, 789, 902]
    actualprob =    [1.0, 1.0, .75, .53, .23, .08, .08, .05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    actualsig_pos = [0.0, 0.0, .11, .15, .15, .13, .13, .12, .09, .09, .09, .09, .09, .09]
    actualsig_neg = [.09, .09, .15, .15, .10, .05, .05, .04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    def __init__(self,sortdat):
        self.name = sortdat.detector1
        self.so_time = sortdat.time
        self.low_cut = sortdat.det1_cut
        self.counts = sortdat.d1counter
        self.rate = self.counts/sortdat.time
        self.energy = sortdat.det1_val
        self.t_energy = sortdat.det1_esum
        self.sig_energy = sortdat.sig_e1
        self.source_time = sortdat.time
     

class Minirad(object):
    codename = 'Tango'
    alarm = 0     # from spec sheet default
    mass = 14.28    # mass of minirad crystal 
    itime = 2    # from spec sheet doesnt really matters
    high_cut = 3000 # total guess
    mass_kg = mass/1000

    ex_values =  [15*pow(2,num) for num in range(8)] #array of values; number above base alarm rate
    levels = np.array(np.arange(2,10)) #alarm lvl at each excess_value
    alarm_con = dict(zip(ex_values,levels)) #dict mapping of ex and lvl
    alarm_con[0] = 1
    actualdist =    [282, 316, 349, 372, 406, 451, 488, 525, 563, 601, 638, 676, 789, 902]
    actualprob =    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, .95, .88, .73, .55, .10, 0.0]
    actualsig_pos = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .04, .07, .11, .14, .13, .09] 
    actualsig_neg = [.09, .09, .09, .09, .09, .09, .09, .09, .12, .14, .15, .15, .06, 0.0]

    def __init__(self,sortdat):
        self.name = sortdat.detector3
        self.source_time = sortdat.time
        self.low_cut = sortdat.det3_cut
        self.high_cut = 3000 
        self.counts = sortdat.d3counter
        self.rate = self.counts/sortdat.time
        self.energy = sortdat.det3_val
        self.t_energy = sortdat.det3_esum
        self.sig_energy = sortdat.sig_e3

    def alarm_lv(self,background):
        excess = self.rate + background.rate - background.base_alarm 
        if excess < 0:
            return 'No Alarm'
        else:
            alarm = self.alarm_con[0]
            rachet = 0
            for val in self.ex_values:
                if excess > val and rachet < val:
                    rachet = val
                    alarm = self.alarm_con[val]

        return alarm

class Bkr_poli(object):
    rate = [7.543,4.9167,2.0317,14.2]
    itime1 = 1
    itime2 = 2
    broundtime = 20*60 #20 min

    def __init__(self,sortdat):
        self.sig_ch1 = np.sqrt(self.rate[0]) #how poli computes sigma  
        self.sig_ch2 = np.sqrt(self.rate[1]) 
        self.sig_ch3 = np.sqrt(self.rate[2]) 
        self.sig_ch4 = np.sqrt(self.rate[3]) 
        self.real_sig = np.sqrt(np.array(self.rate))*self.broundtime #how statistics computes sigma
        self.source_time = sortdat.time

# calculate statistics and make dictionary
def stestat(detector,background):

    stat = OrderedDict() 
    stat['Name'] = detector.name
    stat['actual sigma +'] = detector.actualsig_pos
    stat['actual sigma -'] = detector.actualsig_neg
    stat['Alarm Setting'] = detector.alarm
    stat['Micro Rad Per Hr'] = ( (100*detector.t_energy*conv) / (detector.mass_kg*(detector.source_time/3600)) )*(10**6)    
    stat['Sigma Rad'] = ((100*detector.sig_energy*conv) / (detector.mass_kg*(detector.source_time/3600)) )*(10**6)    
    stat['Source Hits Per Sec'] = detector.rate
    stat['Source Sigma Hit Rate'] = np.sqrt(detector.counts)/detector.source_time 
    stat['Background Hit Rate'] = background.rate
    stat['Sigma Background Hit Rate'] = background.sig_counts
    stat['Combined Hit Rate'] = (detector.rate) + (background.rate)
    stat['Sigma Combined Hit Rate'] = np.sqrt( stat['Source Sigma Hit Rate']**2 + (background.sig_counts)**2 )
    stat['Sigma Above Background'] = detector.rate*detector.itime / background.sig_counts
    stat['Hits Required to alarm'] = stat['Sigma Background Hit Rate']*stat['Alarm Setting'] 
    mean = stat['Combined Hit Rate']
    sigma = stat['Sigma Combined Hit Rate']
    alarm =  background.sig_counts*detector.alarm # value needed to set off alarm
    stat['Probability to Alarm'] = poisson.sf(alarm,mean)

    return stat 

def ministat(detector, background):
    stat = OrderedDict()
    stat['Name'] = detector.name
    stat['Alarm Setting'] = detector.alarm
    stat['Actual Distances'] = detector.actualdist
    stat['Actual Probability'] = detector.actualprob
    stat['Actual Sigma +'] = detector.actualsig_pos
    stat['Actual Sigma -'] = detector.actualsig_neg
    stat['Micro Rad Per Hr'] = ( (100*detector.t_energy*conv) / (detector.mass_kg*(detector.source_time/3600)) )*(10**6)    
    stat['Sigma Rad'] = ((100*detector.sig_energy*conv) / (detector.mass_kg*(detector.source_time/3600)) )*(10**6)    
    stat['Source Hits Per Sec'] = detector.rate
    stat['Source Sigma Hit Rate'] = np.sqrt(detector.counts)/detector.source_time 
    stat['Background Hit Rate'] = background.rate
    stat['Sigma Background Hit Rate'] = background.sig_counts
    stat['Combined Hit Rate'] = (detector.rate) + (background.rate)
    stat['Sigma Combined Hit Rate'] = np.sqrt( stat['Source Sigma Hit Rate']**2 + (background.sig_counts)**2 )
    stat['Base Alarm Count'] = background.base_alarm
    stat['Rate Above Base Alarm Rate'] = stat['Combined Hit Rate'] - stat['Base Alarm Count'] 
    stat['alarm level'] = detector.alarm_lv(background)
    mean = stat['Combined Hit Rate']
    alarm = stat['Base Alarm Count']  # value needed to set off alarm
    stat['Probability to Alarm'] = poisson.sf(alarm,mean)

    return stat


def polistat(poli,background):
    # for one second integration

    stat_1sec = OrderedDict()
    stat_1sec['Actual Distances'] = poli.actualdist
    stat_1sec['Actual Probability'] = poli.actualprob
    stat_1sec['Actual Probability'] = poli.actualprob
    stat_1sec['Actual Sigma +'] = poli.actualsig_pos
    stat_1sec['Actual Sigma -'] = poli.actualsig_neg

    stat_1sec['Name'] = poli.name
    stat_1sec['Title'] = '1 Second Integration time'
    stat_1sec['Alarm Setting'] = poli.alarm

    stat_1sec['Channel 1 Bin'] = poli.ch1_cut
    stat_1sec['Channel 2 Bin'] = poli.ch2_cut
    stat_1sec['Channel 3 Bin'] = poli.ch3_cut
    stat_1sec['Channel 4 Bin'] = poli.ch4_cut

    stat_1sec['Channel 1 Source Hits Per Sec'] = poli.ch1_rate #source rate
    stat_1sec['Channel 2 Source Hits Per Sec'] = poli.ch2_rate 
    stat_1sec['Channel 3 Source Hits Per Sec'] = poli.ch3_rate 
    stat_1sec['Channel 4 Source Hits Per Sec'] = poli.ch4_rate 
    stat_1sec['Source Sigma Hit Rate'] = np.sqrt(poli.ch1_rate*poli.itime1)/poli.s_time  
    stat_1sec['Channel 2 Source Sigma Hit Rate'] = np.sqrt(poli.ch2_rate*poli.itime1)/poli.s_time  
    stat_1sec['Channel 3 Source Sigma Hit Rate'] = np.sqrt(poli.ch3_rate*poli.itime1)/poli.s_time 
    stat_1sec['Channel 4 Source Sigma Hit Rate'] = np.sqrt(poli.ch4_rate*poli.itime1)/poli.s_time 

    stat_1sec['Channel 1 Background Hit Rate'] = background.rate[0] 
    stat_1sec['Channel 2 Background Hit Rate'] = background.rate[1] 
    stat_1sec['Channel 3 Background Hit Rate'] = background.rate[2] 
    stat_1sec['Channel 4 Background Hit Rate'] = background.rate[3] 
    stat_1sec['Channel 1 Sigma Background Hit Rate'] = background.sig_ch1 
    stat_1sec['Channel 2 Sigma Background Hit Rate'] = background.sig_ch2 
    stat_1sec['Channel 3 Sigma Background Hit Rate'] = background.sig_ch3 
    stat_1sec['Channel 4 Sigma Background Hit Rate'] = background.sig_ch4 
    stat_1sec['Regular Background Sigma'] = background.real_sig



    stat_1sec['Combined Hit Rate Channel 1'] = poli.ch1_rate + background.rate[0]
    stat_1sec['Combined Hit Rate Channel 2'] = poli.ch2_rate + background.rate[1]
    stat_1sec['Combined Hit Rate Channel 3'] = poli.ch3_rate + background.rate[2]
    stat_1sec['Combined Hit Rate Channel 4'] = poli.ch4_rate + background.rate[3]
    stat_1sec['Sigma Combined Hit Rate Channel 1'] = np.sqrt( stat_1sec['Channel 1 Sigma Background Hit Rate']**2 + (background.sig_ch1)**2 )
    stat_1sec['Sigma Combined Hit Rate Channel 2'] = np.sqrt( stat_1sec['Channel 2 Sigma Background Hit Rate']**2 + (background.sig_ch2)**2 )
    stat_1sec['Sigma Combined Hit Rate Channel 3'] = np.sqrt( stat_1sec['Channel 3 Sigma Background Hit Rate']**2 + (background.sig_ch3)**2 )
    stat_1sec['Sigma Combined Hit Rate Channel 4'] = np.sqrt( stat_1sec['Channel 4 Sigma Background Hit Rate']**2 + (background.sig_ch4)**2 )


    stat_1sec['Channel 1 Sigma Above Background'] = stat_1sec['Channel 1 Source Hits Per Sec'] / background.sig_ch1
    stat_1sec['Channel 2 Sigma Above Background'] = stat_1sec['Channel 2 Source Hits Per Sec'] / background.sig_ch2
    stat_1sec['Channel 3 Sigma Above Background'] = stat_1sec['Channel 3 Source Hits Per Sec'] / background.sig_ch3
    stat_1sec['Channel 4 Sigma Above Background'] = stat_1sec['Channel 4 Source Hits Per Sec'] / background.sig_ch4
    stat_1sec['Channel 1 Source Rate Required to alarm'] = stat_1sec['Channel 1 Sigma Background Hit Rate']*stat_1sec['Alarm Setting'] 
    stat_1sec['Channel 2 Source Rate Required to alarm'] = stat_1sec['Channel 2 Sigma Background Hit Rate']*stat_1sec['Alarm Setting'] 
    stat_1sec['Channel 3 Source Rate Required to alarm'] = stat_1sec['Channel 3 Sigma Background Hit Rate']*stat_1sec['Alarm Setting'] 
    stat_1sec['Channel 4 Source Rate Required to alarm'] = stat_1sec['Channel 4 Sigma Background Hit Rate']*stat_1sec['Alarm Setting'] 

    ch1_mean = stat_1sec['Combined Hit Rate Channel 1']
    ch2_mean = stat_1sec['Combined Hit Rate Channel 2']
    ch3_mean = stat_1sec['Combined Hit Rate Channel 3']
    ch4_mean = stat_1sec['Combined Hit Rate Channel 4']
    ch1_sigma = stat_1sec['Sigma Combined Hit Rate Channel 1']
    ch2_sigma = stat_1sec['Sigma Combined Hit Rate Channel 2']
    ch3_sigma = stat_1sec['Sigma Combined Hit Rate Channel 3']
    ch4_sigma = stat_1sec['Sigma Combined Hit Rate Channel 4']

    ch1_alarm =  stat_1sec['Channel 1 Source Rate Required to alarm'] # value needed to set off alarm
    ch2_alarm =  stat_1sec['Channel 2 Source Rate Required to alarm'] 
    ch3_alarm =  stat_1sec['Channel 3 Source Rate Required to alarm'] 
    ch4_alarm =  stat_1sec['Channel 4 Source Rate Required to alarm'] 

    stat_1sec['Channel 1 Probability to Alarm'] = poisson.sf(ch1_alarm,ch1_mean)
    stat_1sec['Channel 2 Probability to Alarm'] = poisson.sf(ch2_alarm,ch2_mean)
    stat_1sec['Channel 3 Probability to Alarm'] = poisson.sf(ch3_alarm,ch3_mean)
    stat_1sec['Channel 4 Probability to Alarm'] = poisson.sf(ch4_alarm,ch4_mean)


    # for two second integration time
    stat_2sec = OrderedDict() 
    stat_2sec['Name'] = poli.name
    stat_2sec['Title'] = '2 Second Integration time'
    stat_2sec['Alarm Setting'] = poli.alarm

    stat_2sec['Channel 1 Bin'] = poli.ch1_cut
    stat_2sec['Channel 2 Bin'] = poli.ch2_cut
    stat_2sec['Channel 3 Bin'] = poli.ch3_cut
    stat_2sec['Channel 4 Bin'] = poli.ch4_cut

    stat_2sec['Channel 1 Source Hits Per Sec'] = poli.ch1_rate #source rate
    stat_2sec['Channel 2 Source Hits Per Sec'] = poli.ch2_rate 
    stat_2sec['Channel 3 Source Hits Per Sec'] = poli.ch3_rate 
    stat_2sec['Channel 4 Source Hits Per Sec'] = poli.ch4_rate 
    stat_2sec['Channel 1 Source Sigma Hit Rate'] = np.sqrt(poli.ch1_rate*poli.itime2) / poli.itime2  
    stat_2sec['Channel 2 Source Sigma Hit Rate'] = np.sqrt(poli.ch2_rate*poli.itime2) / poli.itime2  
    stat_2sec['Channel 3 Source Sigma Hit Rate'] = np.sqrt(poli.ch3_rate*poli.itime2) / poli.itime2 
    stat_2sec['Channel 4 Source Sigma Hit Rate'] = np.sqrt(poli.ch4_rate*poli.itime2) / poli.itime2
    stat_2sec['ch1 sigma stat all rate'] =  [ np.sqrt(poli.ch1_rate*poli.itime2) / poli.itime2, np.sqrt(poli.ch2_rate*poli.itime2) / poli.itime2, np.sqrt(poli.ch3_rate*poli.itime2) / poli.itime2, np.sqrt(poli.ch4_rate*poli.itime2) / poli.itime2]


    stat_2sec['Channel 1 Background Hit Rate'] = background.rate[0] 
    stat_2sec['Channel 2 Background Hit Rate'] = background.rate[1]
    stat_2sec['Channel 3 Background Hit Rate'] = background.rate[2]
    stat_2sec['Channel 4 Background Hit Rate'] = background.rate[3]
    stat_2sec['Channel 1 Sigma Background Hit Rate'] = background.sig_ch1 
    stat_2sec['Channel 2 Sigma Background Hit Rate'] = background.sig_ch2 
    stat_2sec['Channel 3 Sigma Background Hit Rate'] = background.sig_ch3 
    stat_2sec['Channel 4 Sigma Background Hit Rate'] = background.sig_ch4 

    stat_2sec['Combined Hit Rate Channel 1'] = poli.ch1_rate + background.rate[0]
    stat_2sec['Combined Hit Rate Channel 2'] = poli.ch2_rate + background.rate[1]
    stat_2sec['Combined Hit Rate Channel 3'] = poli.ch3_rate + background.rate[2]
    stat_2sec['Combined Hit Rate Channel 4'] = poli.ch4_rate + background.rate[3]
    stat_2sec['Sigma Combined Hit Rate Channel 1'] = np.sqrt( stat_2sec['Channel 1 Source Sigma Hit Rate']**2 + (background.sig_ch1)**2 )
    stat_2sec['Sigma Combined Hit Rate Channel 2'] = np.sqrt( stat_2sec['Channel 2 Source Sigma Hit Rate']**2 + (background.sig_ch2)**2 )
    stat_2sec['Sigma Combined Hit Rate Channel 3'] = np.sqrt( stat_2sec['Channel 3 Source Sigma Hit Rate']**2 + (background.sig_ch3)**2 )
    stat_2sec['Sigma Combined Hit Rate Channel 4'] = np.sqrt( stat_2sec['Channel 4 Source Sigma Hit Rate']**2 + (background.sig_ch4)**2 )

    stat_2sec['Channel 1 Sigma Above Background'] = stat_2sec['Channel 1 Source Hits Per Sec'] / background.sig_ch1
    stat_2sec['Channel 2 Sigma Above Background'] = stat_2sec['Channel 2 Source Hits Per Sec'] / background.sig_ch2
    stat_2sec['Channel 3 Sigma Above Background'] = stat_2sec['Channel 3 Source Hits Per Sec'] / background.sig_ch3
    stat_2sec['Channel 4 Sigma Above Background'] = stat_2sec['Channel 4 Source Hits Per Sec'] / background.sig_ch4

    stat_2sec['Channel 1 Source Rate Required to alarm'] = stat_2sec['Channel 1 Sigma Background Hit Rate']*stat_2sec['Alarm Setting'] 
    stat_2sec['Channel 2 Source Rate Required to alarm'] = stat_2sec['Channel 2 Sigma Background Hit Rate']*stat_2sec['Alarm Setting'] 
    stat_2sec['Channel 3 Source Rate Required to alarm'] = stat_2sec['Channel 3 Sigma Background Hit Rate']*stat_2sec['Alarm Setting'] 
    stat_2sec['Channel 4 Source Rate Required to alarm'] = stat_2sec['Channel 4 Sigma Background Hit Rate']*stat_2sec['Alarm Setting'] 

    ch1_mean = stat_2sec['Combined Hit Rate Channel 1']
    ch2_mean = stat_2sec['Combined Hit Rate Channel 2']
    ch3_mean = stat_2sec['Combined Hit Rate Channel 3']
    ch4_mean = stat_2sec['Combined Hit Rate Channel 4']
    ch1_sigma = stat_2sec['Sigma Combined Hit Rate Channel 1']
    ch2_sigma = stat_2sec['Sigma Combined Hit Rate Channel 2']
    ch3_sigma = stat_2sec['Sigma Combined Hit Rate Channel 3']
    ch4_sigma = stat_2sec['Sigma Combined Hit Rate Channel 4']

    ch1_alarm =  stat_2sec['Channel 1 Source Rate Required to alarm'] # value needed to set off alarm
    ch2_alarm =  stat_2sec['Channel 2 Source Rate Required to alarm'] # 
    ch3_alarm =  stat_2sec['Channel 3 Source Rate Required to alarm'] # 
    ch4_alarm =  stat_2sec['Channel 4 Source Rate Required to alarm'] # 
    stat_2sec['Channel 1 Probability to Alarm'] = poisson.sf(ch1_alarm,ch1_mean)
    stat_2sec['Channel 2 Probability to Alarm'] = poisson.sf(ch2_alarm,ch2_mean)
    stat_2sec['Channel 3 Probability to Alarm'] = poisson.sf(ch3_alarm,ch3_mean)
    stat_2sec['Channel 4 Probability to Alarm'] = poisson.sf(ch4_alarm,ch4_mean)

    return stat_1sec, stat_2sec


 #concatanate an array of files for storage
def concat(files, outname):
    with open(outname,'a') as outfile:
        for fname in files:
            with open(fname) as inflile:
                for line in infile:
                    outfile.write(line)
            
 #makes specs 
def plotter(fromdat,filename, show = False):
    
    plt.figure() 
    bins = fromdat.bins
    plt.hist(fromdat.all_val, bins=bins, color=(0, 0, 0, 1 ),
                 histtype='step',label = 'All Hits' )
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title('All Detectors Spectrum\n'+ filename )
    plt.legend(loc='upper right' ) 
    if show:
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
    if show:
        plt.show() 

    his_det3 = plt.hist(fromdat.det3_val, bins=bins, color=(0, 0, 0, 0.5 ),
             histtype='step',label = fromdat.detector3 )
    plt.ylabel('Counts' )
    plt.xlabel('Energy kev' )
    plt.title( fromdat.detector3)
    plt.legend(loc='upper right' ) 
    if show:
        plt.show() 

# fills with data from multiple runs with probs at various distances.
class Ministore(object):
    actual_dist = None
    actual_prob = None
    actualsig_pos = None
    actualsig_neg = None
    
    def __init__(self):
        self.rates = []
        self.prob = []
        self.dist_list = []

    def getactual(self,dict1):
        self.actual_dist = np.array(dict1['Actual Distances'])
        self.actual_prob = np.array(dict1['Actual Probability'])
        self.actualsig_pos = np.array(dict1['Actual Sigma +']) 
        self.actualsig_neg = np.array(dict1['Actual Sigma -'])
        
    def fill(self, lib, distance,sim_time):
        self.dist_list.append(distance)
        self.prob.append(lib['Probability to Alarm'])
        # only load actual dists once 
        if self.actual_dist is None:
            self.getactual(lib)
        # get rate for each call
        self.rates.append(lib['Combined Hit Rate'])
        
# fills with data from multiple runs with probs at various distances.
class Polistore(object):

    dist_list = []
    channels_1 = [33,120,240,745]#from plimaster .inc email
    actual_dist = None
    actual_prob = None
    actualsig_pos = None
    actualsig_neg = None

    def __init__(self):
        self.rates1 = None
        self.rates2 = None
        self.probs1 = None
        self.probs2 = None 
        self.channels = []        
        self.ratesig1 = []
        
    def getchannels(self,dict1,dict2):
        self.channels = [
            dict1['Channel 1 Bin'],
            dict1['Channel 2 Bin'],
            dict1['Channel 3 Bin'],
            dict1['Channel 4 Bin'],
            ]

    def getactual(self,dict1):
        self.actual_dist = np.array(dict1['Actual Distances'])
        self.actual_prob = np.array(dict1['Actual Probability'])
        self.actualsig_pos = np.array(dict1['Actual Sigma +']) 
        self.actualsig_neg = np.array(dict1['Actual Sigma -'])
    
    def fill(self, lib1, lib2, distance, simtime = None):
        self.dist_list.append(distance)
        self.getchannels(lib1,lib2)
        self.ratesig1.append(lib1['Regular Background Sigma']) 
        # only load actuals  once 
        if self.actual_dist is None:
            self.getactual(lib1)
            
        # get probs for each distance simulated 
        probability1 = np.array([
            lib1['Channel 1 Probability to Alarm'],
            lib1['Channel 2 Probability to Alarm'],
            lib1['Channel 3 Probability to Alarm'],
            lib1['Channel 4 Probability to Alarm'],
        ])
        probability2 = np.array([
            lib2['Channel 1 Probability to Alarm'],
            lib2['Channel 2 Probability to Alarm'],
            lib2['Channel 3 Probability to Alarm'],
            lib2['Channel 4 Probability to Alarm'],
        ])
        rate_1 = np.array([
                lib1['Combined Hit Rate Channel 1'],
                lib1['Combined Hit Rate Channel 2'],
                lib1['Combined Hit Rate Channel 3'],
                lib1['Combined Hit Rate Channel 4'],
        ])
        rate_2 = np.array([
                lib2['Combined Hit Rate Channel 1'],
                lib2['Combined Hit Rate Channel 2'],
                lib2['Combined Hit Rate Channel 3'],
                lib2['Combined Hit Rate Channel 4'],
                ])
        #appends new 4 elem array to matrix (unless it's first time)
        if self.rates1 is None:
            self.rates1 = rate_1
            self.rates2 = rate_2
            self.prob1 = probability1
            self.prob2 = probability2
        else:
            self.rates1 = np.vstack((self.rates1,rate_1))
            self.rates2 = np.vstack((self.rates2,rate_1))

            self.prob1 = np.vstack((self.prob1,probability1))
            self.prob2 = np.vstack((self.prob2,probability2))
            
    def __str__(self):
        printout = []
        for dist in self.rates1:
            printout.append(str(dist))
        return ''.join(printout)

        
def splinyplt(dist,prob): 
    chan1_interp = interpolate.interp1d(dist,prob,3)
    fill_xspace = np.linspace(np.amin(dist), np.amax(dist),100)
    prob_intrp = chan1_interp(fill_xspace)

    return fill_xspace, prob_intrp

def sigmoid(x,x0,k):
    y = 1 / ( 1 + np.exp(-1*( k + x0*x )) )
    return y










