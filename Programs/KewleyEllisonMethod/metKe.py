'''Import data file, run metallicity calculations on all galaxies included, and 
create text file of metallicity values.'''

################################################################################
#
#								IMPORT LIBRARIES
#
################################################################################


from metal_PP04 import PP04
#from metal_Izotov06 import Izotov06
from metal_D02 import D02
from metal_pyNeb import pyMet, pyMet_error
#from Izotov06_error import I06_error
from astropy.table import Table
import numpy as np
import pylab as P
#from time import time

# Cython code
import sys
sys.path.insert(1, '/Users/kellydouglass/Documents/Drexel/Research/Programs/KewleyEllisonMethod/cythonized_code/')
from metal_Izotov06_Cython import Izotov06
from Izotov06_error_Cython import I06_error


################################################################################
#
#								DEFINE FUNCTION
#
################################################################################


def metKE(method, outFile, extra_cols):
    '''method is the metallicity method to be used
    outFile is the name of the output file for the calculated metallicity values
    extra_cols is the number of columns of additional data in the orginal file 
    (before the beginning of the flux data)'''

    ############################################################################
    #
    #								INITIALIZATIONS
    #
    ############################################################################
    
    # Dummy variables (error code not yet written)
    PP04_error = 0
    D02_error = 0

    # Dictionary to relate method to function name
    method_options = {'I06':[Izotov06, I06_error], 'pyNeb':[pyMet, pyMet_error], 'PP04':[PP04, PP04_error], 'D02':[D02, D02_error]}
    
    # Output table
    out_table = Table()
    
    # Output lists
    Z_list = []
    Zerror_list = []
    logOppHp12_list = []
    logOpHp12_list = []
    t3_list = []
    logNH12_list = []
    NHerror_list = []
    logNpHp12_list = []
    logNeH12_list = []
    NeHerror_list = []
    logNeppHp12_list = []
    
    
    ############################################################################
    #
    #									INPUTS
    #
    ############################################################################
    
    '''# obtain method from user
    method = str(input('Which metallicity method? (I06, pyNeb, PP04, D02): '))
    outFile = str(input('Output file name (with extension): '))
    '''

    fileName = method + '.txt'
    method_calc = method_options[method][0]
    method_error = method_options[method][1]
    
    
    ############################################################################
    #
    #	IMPORT DATA
    #
    ############################################################################
    
    
    galaxy_data = Table.read(fileName, format='ascii.commented_header')
    
    
    ############################################################################
    #
    # ADD INITIAL COLUMNS TO OUTPUT TABLE
    #
    ############################################################################
    
    
    out_table['index'] = galaxy_data['index']
    
    for column in extra_cols:
        out_table[column] = galaxy_data[column]
    
    if method=='I06':
        out_table['flag3727'] = galaxy_data['flag3727']
        out_table['flagOppSp'] = galaxy_data['flagOppSp']
        #out_table['flagS2'] = galaxy_data['flagS2']
        #out_table['S2ratio'] = galaxy_data['S2ratio']
        #out_table['flagO3'] = galaxy_data['flagO3']
        #out_table['O3ratio'] = galaxy_data['O3ratio']
    
    
    ############################################################################
    #
    #							CALCULATE METALLICITY
    #
    ############################################################################
    
    
    # read through file and calculate metallicity for each data set
    for i_galaxy in range(len(galaxy_data)):
        
        if i_galaxy%100 == 0:
            print 'Calculating galaxy #', i_galaxy
        #print '****  Galaxy ID', galaxy_data['index'][i_galaxy]
        
        ########################################################################
        #							CALCULATE ABUNDANCES
        ########################################################################
        
        '''# Initialize temperature calculation flag
        bad_temp = [0]'''
        
        #t = time()
        # Calculate metallicity
        Z,logOppHp12,logOpHp12,t3,logNH12,logNpHp12,logNeH12,logNeppHp12 = method_calc(galaxy_data[i_galaxy])
        #print 'Abundance calculation:', time() - t

        '''if bad_temp[0]:
            print galaxy_data['index'][i_galaxy], 'calculation'
        '''
        
        # Calculate abundance errors
        '''
        Zerror = np.nan
        NHerror = np.nan
        '''
        
        #t = time()
        if np.isfinite(Z):
            Zerror, NHerror, NeHerror = method_error(galaxy_data[i_galaxy])
            if galaxy_data['flagNH'][i_galaxy] == 0:
                NHerror = np.nan
            if galaxy_data['flagNeH'][i_galaxy] == 0:
                NeHerror = np.nan
        else:
            Zerror = np.nan
            NHerror = np.nan
            NeHerror = np.nan
        #print '    Error calculation:', time() - t
        
        '''
        if galaxy_data['flag3727'][i_galaxy]:
            Zerror, NHerror = method_error(galaxy_data[i_galaxy])
            
            if galaxy_data['flagNH'][i_galaxy] == 0:
                NHerror = np.nan
        else:
            Zerror = np.nan
            NHerror = np.nan
        '''
        
        '''# Display histogram of metallicity values (only for complete estimates)
        if galaxy_data['flag_3727'][i_galaxy] == 1:
            fig, ax = P.subplots(1)
            n, bins, patches = ax.hist(met, 50, normed=1, histtype='stepfilled')
            y = P.normpdf(bins, Z, error)
            ax.plot(bins, y, 'k--', linewidth=1.5)
            ax.axvline(x=Z, linestyle='--', c='r', linewidth=1.5)
            txtstr = '$Z = %.2f$\n$\sigma = %.3f$'%(Z, error)
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.05, 0.95, txtstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
            P.show()
        '''
        
        ########################################################################
        #							SAVE RESULTS TO LISTS
        ########################################################################
        
        
        Z_list.append(Z)
        Zerror_list.append(Zerror)
        logOppHp12_list.append(logOppHp12)
        logOpHp12_list.append(logOpHp12)
        t3_list.append(t3)
        logNH12_list.append(logNH12)
        NHerror_list.append(NHerror)
        logNpHp12_list.append(logNpHp12)
        logNeH12_list.append(logNeH12)
        NeHerror_list.append(NeHerror)
        logNeppHp12_list.append(logNeppHp12)


    ############################################################################
    #
    #								BUILD OUTPUT TABLE
    #
    ############################################################################


    out_table['Z12logOH'] = Z_list
    out_table['Zerr'] = Zerror_list
    out_table['Z12logOppH'] = logOppHp12_list
    out_table['Z12logOpH'] = logOpHp12_list
    out_table['t3'] = t3_list
    out_table['N12logNH'] = logNH12_list
    out_table['NHerr'] = NHerror_list
    out_table['N12logNpH'] = logNpHp12_list
    out_table['Ne12logNeH'] = logNeH12_list
    out_table['NeHerr'] = NeHerror_list
    out_table['Ne12logNeppH'] = logNeppHp12_list
    
    
    ############################################################################
    #
    #	SAVE OUTPUT TABLE
    #
    ############################################################################
    
    
    out_table.write(outFile, format='ascii.commented_header', overwrite=True)
