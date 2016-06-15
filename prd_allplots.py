
# Plots ste,poli, and Minirad data from many runs
import prd_auxfunctions2 as aux
import numpy as np
import matplotlib.mlab as mlab 
import matplotlib.pyplot as plt 


# requres a Polistore instance poli
def plotpoli(poli, rate = False, prob = True):
    err_xpos = poli.actualsig_pos
    err_xneg = poli.actualsig_neg
    distances = np.array(poli.dist_list)
    exp_sec = 2

    #plot rate if plot = True
    if rate:
        plt.figure() 
        plt.title('1 Second Integration Polimaster')
        plt.ylabel('Count Rate' )
        plt.xlabel('Distance CM' )
        plt.plot(distances, poli.rates1[:,0],
                label = 'Channel 1 :  '+str(poli.channels[0]) + 'kev - '+str(poli.channels[1]) + ' kev')
        plt.plot(distances, poli.rates1[:,1],
                label = 'Channel 2 : '+ str(poli.channels[1])+ ' kev - ' + str(poli.channels[2]) + ' kev')
        plt.plot(distances, poli.rates1[:,2],
                label = 'Channel 3 : '+ str(poli.channels[2])+ ' kev - ' + str(poli.channels[3])+ ' kev')
        plt.plot(distances, poli.rates1[:,3],
                label = 'Channel 4 : '+ str(poli.channels[0])+ ' kev - ' + str(poli.channels[3])+ ' kev')
        plt.legend(loc=1, borderaxespad=0.,frameon = False)
        plt.axis([300,600,0,1.1])
        plt.show()

        plt.figure() 
        plt.title('2 Second Integration Polimaster')
        plt.ylabel('Count Rate' )
        plt.xlabel('Distance CM' )
        plt.plot(distances, poli.rates2[:,0],
                label = 'Channel 1 : ' +str(poli.channels[0]) + ' kev ' +str(poli.channels[1]) + ' kev'
                )
        plt.plot(distances, poli.rates2[:,1],
                label = 'Channel 2 : '+ str(poli.channels[1])+ ' kev - ' + str(poli.channels[2]) + ' kev'
                )
        plt.plot(distances, poli.rates2[:,2],
                label = 'Channel 3 : '+ str(poli.channels[2])+ ' kev - ' + str(poli.channels[3])+ ' kev'
                )
        plt.plot(distances, poli.rates2[:,3],
                label = 'Channel 4 : '+ str(poli.channels[0])+ ' kev - ' + str(poli.channels[3])+ ' kev'
                )
        plt.legend(loc=1, borderaxespad=0.,frameon = False)
        plt.axis([300,600,0,1.1])
        plt.show()
    #plot prob to dectect 
    if prob:
        plt.figure() 
        plt.title('1 Second Integration Polimaster\nProbability of Detection in Less Than '+str(exp_sec)+' Seconds')
        plt.ylabel('Probability of Detection' )
        plt.xlabel('Distance CM From Ba133 Source' )

        ch1prob1 = (1 - (1 - poli.prob1[:,0])**exp_sec ) 
        ch2prob1 = (1 - (1 - poli.prob1[:,1])**exp_sec )
        ch3prob1 = (1 - (1 - poli.prob1[:,2])**exp_sec )
        ch4prob1 = (1 - (1 - poli.prob1[:,3])**exp_sec )


        #realx, realy = aux.splinyplt( poli.actual_dist[4:], poli.actual_prob[4:])
        #x_ch1p1, y_ch1pb1 = aux.splinyplt(distances,ch1prob1)
        #x_ch2p1, y_ch2pb1 = aux.splinyplt(distances,ch2prob1)
        #x_ch3p1, y_ch3pb1 = aux.splinyplt(distances,ch3prob1)
        #x_ch4p1, y_ch4pb1 = aux.splinyplt(distances,ch4prob1)
        
        all_errors = [err_xneg[4:], err_xpos[4:]]
#        plt.errorbar(poli.actual_dist, poli.actual_prob,yerr=all_errors, marker = 'x', linewidth = '3', c = 'r', label = 'Actual Result')
        plt.errorbar( poli.actual_dist[4:], poli.actual_prob[4:], yerr=all_errors, marker = 'x', linewidth = '3',c = 'r', label = 'Actual Result')

        plt.plot(distances, ch1prob1 ,marker = 'o', label = 'Channel 1 :' +str(poli.channels[0]) +'kev - '+str(poli.channels[1]) + ' kev')
        plt.plot(distances, ch2prob1, marker = 'o', label = 'Channel 2 : '+ str(poli.channels[1])+ ' kev - ' + str(poli.channels[2]) + ' kev')
        plt.plot(distances, ch3prob1, marker = 'o', label = 'Channel 3 : '+ str(poli.channels[2])+ ' kev - ' + str(poli.channels[3])+ ' kev')
        plt.plot(distances, ch4prob1, marker = 'o', label = 'Channel 4 : '+ str(poli.channels[0])+ ' kev - ' + str(poli.channels[3])+ ' kev')
        plt.legend(loc=1, borderaxespad=0.,frameon = False)
        plt.show()

        # 2 second Integration 
        plt.figure() 
        plt.title('2 Second Integration Polimaster\nProbability of Detection in Less Than '+str(exp_sec)+' Seconds')
        plt.ylabel('Probability of Detection' )
        plt.xlabel('Distance in CM to Ba133 Source' )
        ch1prob2 = 1 - (1 - poli.prob2[:,0])**2
        ch2prob2 = 1 - (1 - poli.prob2[:,1])**2
        ch3prob2 = 1 - (1 - poli.prob2[:,2])**2
        ch4prob2 = 1 - (1 - poli.prob2[:,3])**2

        #realx, realy = splinyplt( poli.actual_dist[4:], poli.actual_prob[4:])
        #x_ch1p2, y_ch1pb2 = splinyplt(distances,ch1prob2)
        #x_ch2p2, y_ch2pb2 = splinyplt(distances,ch2prob2)
        #x_ch3p2, y_ch3pb2 = splinyplt(distances,ch3prob2)
        #x_ch4p2, y_ch4pb2 = splinyplt(distances,ch4prob2)
        all_errors = [err_xneg, err_xpos]
        plt.errorbar(poli.actual_dist, poli.actual_prob,yerr=all_errors, marker = 'x', linewidth = '3', c = 'r', label = 'Actual Result')

        plt.plot( distances, ch1prob2, marker = 'o', label = 'Channel 1 :' +str(poli.channels[0]) +'kev - '+str(poli.channels[1]) + ' kev')
        plt.plot( distances, ch2prob2, marker = 'o', label = 'Channel 2 : '+ str(poli.channels[1])+ ' kev - ' + str(poli.channels[2]) + ' kev')
        plt.plot( distances, ch3prob2, marker = 'o', label = 'Channel 3 : '+ str(poli.channels[2])+ ' kev - ' + str(poli.channels[3])+ ' kev')
        plt.plot( distances, ch4prob2, marker = 'o', label = 'Channel 4 : '+ str(poli.channels[0])+ ' kev - ' + str(poli.channels[3])+ ' kev')
        plt.legend(loc=1, borderaxespad=0., frameon = False)
        plt.show()

# requires Ministore object whcih holds probs and rates of several runs
def plotmini(mini, rate = False, prob = True):
    err_xpos = mini.actualsig_pos
    err_xneg = mini.actualsig_neg
    distances = np.array(mini.dist_list)
    rates = np.array(mini.rates)
    probs = np.array(mini.prob)
    exp_sec = 2

    if True:
        plt.figure() 
        plt.title('Minirad D no res')
        plt.ylabel('Count Rate' )
        plt.xlabel('Distance CM' )
        print(rates)
        print(len(rates))

        plt.plot(distances, rates, label = 'minirad')
        plt.legend(loc=1, borderaxespad=0.,frameon = False)
        plt.show()
    if prob:
        plt.figure() 
        plt.title('Minirad\nProbability of Detection in Less Than '+str(exp_sec)+' Seconds')
        plt.ylabel('Probability of Detection' )
        plt.xlabel('Distance CM From Ba133 Source' )

        prob2alrm = (1 - (1 - probs)**exp_sec ) 

        all_errors = [err_xneg, err_xpos]
        plt.errorbar( mini.actual_dist, mini.actual_prob, yerr=all_errors, marker = 'x', linewidth = '3', c = 'r', label = 'Actual Result')
        plt.plot( distances, prob2alrm, 'o', label = 'Minirad')
        plt.legend(loc=1, borderaxespad=0., frameon = False)
        plt.show()


