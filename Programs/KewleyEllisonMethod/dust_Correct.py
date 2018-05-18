################################################################################
#
#   IMPORT LIBARIES
#
################################################################################


from astropy.table import Table
import pyneb as pn


################################################################################
#
#   DEFINE FUNCTION
#
################################################################################


def dustCorr_pyNeb(method, extra_cols):
    '''Correct emission lines for dust extinction using pyNeb's implementation 
    of the Cardelli89 extinction curve
    
    Assumptions:
      - Ha/Hb = 2.86 at 10,000 K at 100 cm^-3
    '''
    
    
    ############################################################################
    #   IMPORT DATA
    ############################################################################
    
    
    flux_table = Table.read(method+'.txt', format='ascii.commented_header')
    
    # Which lines need to be de-reddened?
    colNames = flux_table.colnames
    
    # Add 'index' to extra_cols
    extra_cols.append('index')
    if method == 'I06':
        extra_cols.append('flag3727')
        extra_cols.append('flagOppSp')
        extra_cols.append('flagNH')
        extra_cols.append('flagNeH')
        extra_cols.append('flagS2')
        extra_cols.append('S2ratio')
        extra_cols.append('flagO3')
        extra_cols.append('O3ratio')
    
    # Remove all column names that are not 'index' or in extra_cols
    emLines = [line for line in colNames if line not in extra_cols]
    
    
    ############################################################################
    #   INITIALIZE EXTINCTION CORRECTION
    ############################################################################
    
    
    RC = pn.RedCorr(law='CCM89')
    
    
    ############################################################################
    #   CORRECT FOR DUST EXTINCTION
    ############################################################################
    
    
    # observed H_alpha / H_beta ratio
    HaHb_obs = flux_table['H_ALPHA_FLUX']/flux_table['H_BETA_FLUX']
    
    # Calculate correction coefficient based on theoretical H_alpha/H_beta ratio
    RC.setCorr(HaHb_obs/2.86, 6563., 4861.)
    
    
    for line in emLines[::2]:
        # Extract wavelength from line
        waveLength = line.split('_')[1]
        if waveLength == 'BETA':
            waveLength = 4861.
        elif waveLength == 'ALPHA':
            waveLength = 6563.
        else:
            waveLength = float(waveLength)
        
        # De-redden emission line
        flux_table[line] = flux_table[line]*RC.getCorrHb(waveLength)
    
    
    ############################################################################
    #   SAVE DATA
    ############################################################################
    
    
    flux_table.write(method+'.txt', format='ascii.commented_header', overwrite=True)