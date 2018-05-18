'''Save off galaxies to file for calculations.  NO FILTERING'''

################################################################################
#
#								IMPORT LIBRARIES
#
################################################################################


import numpy as np
from astropy.table import Table


################################################################################
#
#								DEFINE FUNCTION
#
################################################################################


def totalCount_KE_noFilter(fileName, extra_cols):
    '''fileName is name of file with flux data
    extra_cols is the number of additional columns of data (before flux data begins)'''
    
    # fileName = 'I06_galaxyData.txt'
    
    ############################################################################
    #								IMPORT DATA
    ############################################################################
    
    
    galaxy_data = Table.read(fileName, format='ascii.commented_header')
    
    # Number of galaxies
    N_gal = len(galaxy_data)
    
    
    ############################################################################
    # 								OUTPUT TABLES
    ############################################################################
    
    
    # Output tables
    I06_table = Table()
    PP04_table = Table()
    D02_table = Table()
    pyNeb_table = Table()
    
    
    ############################################################################
    #	DETERMINE [OII] 3727 FLAG
    ############################################################################
    
    
    flag_3727_boolean = galaxy_data['OII_3727_FLUX'] > 0.
    flag_3727_list = flag_3727_boolean.astype(int)
    
    
    ############################################################################
    #	DETERMINE N/H FLAG
    ############################################################################
    
    
    flag_NH_boolean = galaxy_data['NII_6584_FLUX'] > 0.
    flag_NH_list = flag_NH_boolean.astype(int)
    
    
    ############################################################################
    #	DETERMINE Ne/H FLAG
    ############################################################################
    
    
    flag_NeH_boolean = galaxy_data['NeIII_3869_FLUX'] > 0.
    flag_NeH_list = flag_NeH_boolean.astype(int)
    
    
    ############################################################################
    #   INSERT DUMMY BPTCLASS COLUMN (IF NEEDED)
    ############################################################################
    
    
    if 'BPTclass' not in galaxy_data.colnames:
        galaxy_data['BPTclass'] = np.ones(len(galaxy_data))
    
    
    ############################################################################
    #
    #	BUILD OUTPUT TABLES
    #
    ############################################################################
    
    
    # Suffixes needed
    suffixes = ['_FLUX', '_FLUX_ERR']
    
    
    ################################	I06		################################
    
    # Lines needed
    I06_lines = ['OII_3727', 'NeIII_3869', 'OIII_4363', 'H_BETA', 'OIII_4959', 'OIII_5007', 'NII_6548', 'H_ALPHA', 'NII_6584']
    
    if 'OII_7325_FLUX' in galaxy_data.colnames:
        I06_lines.append('OII_7325')
    
    I06_table['index'] = galaxy_data['index']
    
    for column in extra_cols:
        I06_table[column] = galaxy_data[column]
    
    I06_table['flag3727'] = flag_3727_list
    I06_table['flagNH'] = flag_NH_list
    I06_table['flagNeH'] = flag_NeH_list
    
    for line in I06_lines:
        for ending in suffixes:
            I06_table[line + ending] = galaxy_data[line + ending]
    
    
    ################################	PP04	################################
    
    # Lines needed
    PP04_lines = ['H_BETA', 'OIII_5007', 'H_ALPHA', 'NII_6584']
    
    PP04_table['index'] = galaxy_data['index']
    
    for column in extra_cols:
        PP04_table[column] = galaxy_data[column]
    
    for line in PP04_lines:
        for ending in suffixes:
            PP04_table[line + ending] = galaxy_data[line + ending]
    
    
    ################################	D02		################################
    
    # Lines needed
    D02_lines = ['H_ALPHA', 'NII_6584']
    
    D02_table['index'] = galaxy_data['index']
    
    for column in extra_cols:
        D02_table[column] = galaxy_data[column]
    
    for line in D02_lines:
        for ending in suffixes:
            D02_table[line + ending] = galaxy_data[line + ending]
    
    	
    ################################	pyNeb	################################
    
    # Lines needed
    pyNeb_lines = ['OII_3727', 'NeIII_3869', 'OIII_4363', 'H_BETA', 'OIII_4959', 'OIII_5007', 'NII_6548', 'H_ALPHA', 'NII_6584']
    
    if 'SII_6717_FLUX' in galaxy_data.colnames:
        pyNeb_lines.append('SII_6717')
        pyNeb_lines.append('SII_6731')
    elif 'SII_6725_FLUX' in galaxy_data.colnames:
        pyNeb_lines.append('SII_6725')
    
    pyNeb_table['index'] = galaxy_data['index']
    
    for column in extra_cols:
        pyNeb_table[column] = galaxy_data[column]
    
    pyNeb_table['flag3727'] = flag_3727_list
    
    for line in pyNeb_lines:
        for ending in suffixes:
            pyNeb_table[line + ending] = galaxy_data[line + ending]
    
    
    
    ############################################################################
    #
    #							WRITE OUTPUT FILES
    #
    ############################################################################
    
    
    I06_table.write('I06.txt', format='ascii.commented_header', overwrite=True)
    PP04_table.write('PP04.txt', format='ascii.commented_header', overwrite=True)
    D02_table.write('D02.txt', format='ascii.commented_header', overwrite=True)
    pyNeb_table.write('pyNeb.txt', format='ascii.commented_header', overwrite=True)