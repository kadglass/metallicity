'''Main script to calculate metallicity of galaxy(ies)
Data comes from text file outputed by MPAJHU_extract_fromIndex.py'''


################################################################################
#
#	IMPORT LIBRARIES
#
################################################################################


from totalCount_KE_MPAJHU import totalCount_KE_MPAJHU
from totalCount_KE_Portsmouth import totalCount_KE_Portsmouth
from totalCount_KE_noFilter import totalCount_KE_noFilter
from dust_Correct import dustCorr_pyNeb
from metKe import metKE
from Z_calc_pyNeb import metpy
from time import time


################################################################################
#
#	INPUTS FROM USER
#
################################################################################
start_time = time()


# Metallicity method to use
method = str(input('Which metallicity method? (I06, pyNeb, PP04, D02): '))
#method = 'pyNeb'

# S/N restriction on [OIII] 4363
sigma_4363 = float(input('Minimum S/N of [OIII] 4363: '))

# File name of data
#fileName = input('File name of data (with extension): ')
fileName = '../../Data/dwarf_voidstatus_KE_MPAJHU_flux_oii.txt'
#fileName = '../../Data/dwarf_voidstatus_stellarMass_KE_MPAJHU_flux_oii.txt'
#fileName = '../../Data/dwarf_voidstatus_BPT_KE_MPAJHU_flux_oii.txt'
#fileName = '../../Data/dwarf_voidstatus_BPT_KE_Portsmouth_flux.txt'
#fileName = '../../Data/kias1033_5_KE_MPAJHU_flux_oii.txt'
#fileName = '../../Data/kias1033_5_BPT_KE_MPAJHU_flux_oii.txt'
#fileName = '../../Data/kias1033_5_BPT_KE_Portsmouth_flux.txt'
#fileName = '../../Data/kias1033_5_stellarMass_KE_MPAJHU_flux_oii.txt'
#fileName = '../../Data/I06_galaxyData.txt'
#fileName = '../../Data/Yin07_KE_MPAJHU_flux_oii.txt'
#fileName = '../../Data/Yin07_KE_MPAJHU_flux.txt'
#fileName = '../../Data/SHELS/Data/SHELSgalaxies_KE_flux.txt'

# List of columns of additional data (does not include index)
extra_cols = ['rabsmag', 'BPTclass', 'vflag', 'MPA_index']
#extra_cols = ['BPTclass']

################################################################################
#
#	CREATE OUTPUT FILE NAME
#
################################################################################


# Filter fileName for outFile
if fileName[0] == '.':
	for i in range(len(fileName)-1,0,-1):
		if fileName[i] == '/':
			break
			
	outFileName = fileName[i+1:]


# Build output file name
outFile = 'Cython_bisect_comp_Z_KE_' + method + 'ratioApprox_' + str(int(sigma_4363)) + 'sig_' + outFileName


################################################################################
#
#	FILTER GALAXIES
#
################################################################################
print 'Filtering galaxies'
print '-----------------------------'

'''
# Filter out galaxies that do not have valid spectra
if 'MPAJHU' in fileName:
	totalCount_KE_MPAJHU(fileName, extra_cols, sigma_4363)
elif 'Portsmouth' in fileName:
	totalCount_KE_Portsmouth(fileName, extra_cols, sigma_4363)
else:
	print 'Unknown or unidentified data type in file name.'
	exit()
'''
'''
# For testing code
totalCount_KE_noFilter(fileName, extra_cols)
outFile = 'Cython_bisect_comp_Z_KE_' + method + '_' + outFileName
'''


t = time()
print 'Filtering time:', t - start_time
################################################################################
#
#	CALCULATE METALLICITIES
#
################################################################################
print 'Calculating metallicities'
print '-----------------------------'


# Correct line fluxes for dust extinction
dustCorr_pyNeb(method, extra_cols)
outFile = outFile[:-4] + '_dustCorr.txt'


# Calculate metallicities of remaining galaxies
if method=='pyNeb':
	metpy(outFile, extra_cols)
else:
	metKE(method, outFile, extra_cols)


print 'Calculation time:', time() - t
print 'Total runtime:', time() - start_time
