'''Estimate the error in the metallicity calculation for the I06 method (Izotov 
2006)'''


################################################################################
#
#								IMPORT LIBRARIES
#
################################################################################


from metal_Izotov06_Cython import Izotov06
from random import normalvariate
import numpy as np
from multiprocessing import Process, cpu_count, Queue
from math import ceil

from time import time

import matplotlib.pyplot as plt


################################################################################
#
#							GENERATE FLUXES FUNCTION
#
################################################################################


cpdef rand_fluxes(galaxy_info, line_list, n, sample_method):

    flux_values = {}
    
    for line in line_list:
    
        mean = galaxy_info[line]
        sigma = galaxy_info[line + '_ERR']
        
        if sample_method == 'resample_normal':
        
            samples = np.random.normal(mean, sigma, size=n)
            
            while np.any(samples < 0):
                fillers = np.random.normal(mean, sigma, size=np.sum(samples < 0))
                samples[samples < 0] = fillers
                
        elif sample_method == 'gamma':
            
            shape = (mean*mean)/(sigma*sigma)
            scale = (sigma*sigma)/mean
            
            samples = np.random.gamma(shape, scale, size=n)
            
        elif sample_method == 'abs_normal':
        
            samples = np.abs(np.random.normal(mean, sigma, size=n))
            
        flux_values[line] = samples
        
    return flux_values
                

    


cpdef I06_error_process(galaxy, results_queue, n, sentinel):

    ############################################################################
    #   Initialization
    ############################################################################

    # Emission lines to generate
    em_lines = ['OII_3727_FLUX', 'NeIII_3869_FLUX', 'OIII_4363_FLUX', 'H_BETA_FLUX', 'OIII_4959_FLUX', 'OIII_5007_FLUX', 'NII_6548_FLUX', 'NII_6584_FLUX', 'SII_6717_FLUX', 'SII_6731_FLUX']
    
    if (galaxy['flagNeH'] == 0):
        em_lines.pop(1)
    
    if (galaxy['flagNH'] == 0):
        em_lines = em_lines.pop(-4)
        em_lines = em_lines.pop(-3)

    if (galaxy['flag3727'] == 0):
        em_lines = em_lines[1:]

    if (galaxy['flagOppSp'] == 0):
        em_lines = em_lines[:-2]

    # Initialize fake flux dictionary
    F = {'index':galaxy['index'], 'BPTclass':galaxy['BPTclass'], 'rabsmag':galaxy['rabsmag'], 'flag3727':galaxy['flag3727'], 'flagOppSp':galaxy['flagOppSp'], 'flagNH':galaxy['flagNH'], 'flagNeH':galaxy['flagNeH']}

    ############################################################################
    #   Generate fluxes
    ############################################################################
    
    sample_method = "resample_normal"
    
    flux_values = rand_fluxes(galaxy, em_lines, n, sample_method)
    
    ############################################################################
    #   Calculate abundances
    ############################################################################
    
    cdef int i, tries

    for i in xrange(n):
    
        ########################################################################
        #   Extract fluxes
        ########################################################################
        
        for line in em_lines:
            F[line] = flux_values[line][i]
            
        ########################################################################
        #   Calculate and store abundances
        ########################################################################
        
        met,_,_,_,logNH,_,logNeH,_ = Izotov06(F)
        
        if np.isfinite(met) and np.isfinite(logNH) and np.isfinite(logNeH):
            results_queue.put((met, logNH, logNeH))
            
        else:
            ####################################################################
            #   Re-run nan metallicities
            ####################################################################
    
            tries = 0
            
            while (np.isnan(met) or ((np.isnan(logNH) and galaxy['flagNH'] == 1) or (np.isnan(logNeH) and galaxy['flagNeH'] == 1))) and tries < 500:
            
                tries = tries + 1
        
                ################################################################
                # Generate flux
                ################################################################
        
                new_fluxes = rand_fluxes(galaxy, em_lines, 1, sample_method)
                
                for line in em_lines:
                    F[line] = new_fluxes[line]
            	
                ################################################################
                # Calculate abundances
                ################################################################
        
                met,_,_,_,logNH,_,logNeH,_ = Izotov06(F)
        
            ####################################################################
            #   Store abundances
            ####################################################################
            
            results_queue.put((met, logNH, logNeH))
            
            if tries == 500:
                print '****  Galaxy', galaxy['index'], 'failed to generate all real abundances in error calculation'
                
    
    # Return sentinel value (all done!)
    results_queue.put(sentinel)
    
    return None


################################################################################
#
#								ERROR FUNCTION
#
################################################################################


def I06_error(galaxy):
    '''Generates N fake data sets from a normal distribution with the mean 
    equalling the actual observation (fluxes) and the dispersion equalling the 
    reported estimate of the uncertainty in the observation (errors).  The 
    dispersion of the resulting metallicities serves as a fair estimate of the 
    uncertainty in the calculated metallicity.'''
    
    ############################################################################
    #	Initializations & Constants
    ############################################################################
    
    
    # number of data points
    N = 100000
    
    # initialize abundance arrays
    met = []
    logNH = []
    logNeH = []

    ############################################################################
    #	Initialize parallelization
    ############################################################################
    
    # Count number of logical cores available
    num_cores = cpu_count()
    
    # Number of active cores
    num_active = num_cores
    
    # Initialize process reference list
    processes = []
    
    # Initialize results queue
    results_queue = Queue()
    
    # Process termination sentinel
    sentinel = 'poop'
    
    # Number of runs per process
    runs = int(ceil(N/num_cores))
    
    # Initiate processes
    for _ in range(num_cores):
        P = Process(target=I06_error_process, args=(galaxy, results_queue, runs, sentinel))
        processes.append(P)
        
        P.start()
    
    ############################################################################
    #	Calculate metallicities
    ############################################################################
    
    while num_active > 0:
        
        try:
            # Extract next available output from the queue
            result_value = results_queue.get()
        except Queue.Empty:
            # Queue is empty (either all processes are finished or are all calculating)
            continue
        
        if result_value == sentinel:
            # Process has finished
            num_active -= 1
        else:
            # Unpack output tuple
            Z12logOH, N12logNH, Ne12logNeH = result_value
            
            # Append values to lists
            met.append(Z12logOH)
            logNH.append(N12logNH)
            logNeH.append(Ne12logNeH)
    
    ############################################################################
    #	Join processes
    ############################################################################
    
    for p in processes:
        p.join(None)
    
    ############################################################################
    #	Determine dispersion of abundance measurements
    ############################################################################
    
    '''if bad_temp[0]:
        print galaxy['index'], 'error'
    '''
    
    Zerr = np.std(met)
    NHerr = np.std(logNH)
    NeHerr = np.std(logNeH)
    
    '''
    ############################################################################
    #	Plot histogram of abundance measurements
    ############################################################################
    
    plt.hist(met, bins=50)
    plt.title('Galaxy #'+str(galaxy['index'])+' (I06)')
    '''
    
    ############################################################################
    #	Return error estimates
    ############################################################################
    
    return Zerr, NHerr, NeHerr
