'''Filter through database to determine if all line fluxes are present for 
various Kewley-Ellison metallicity calculation methods.  Designed for data from 
Portsmouth value-added catalog.'''

################################################################################
#
#								IMPORT LIBRARIES
#
################################################################################


from FluxCount import count
import numpy as np
from astropy.table import Table
import sys

sys.path.insert(1, '/Users/kellydouglass/Documents/Research/Data/')

from Reset_flag3727 import reset_flag3727


################################################################################
#
#								DEFINE FUNCTION
#
################################################################################


def totalCount_KE_Portsmouth(fileName, extra_cols, sigma_4363):
	'''fileName is name of file with flux data
	extra_cols is the number of additional columns of data (before flux data begins)
	sigma_4363 is the minimum S/N restriction to implement on [OIII] 4363'''
	
	
	############################################################################
	#								IMPORT DATA
	############################################################################
	# set up table with line intensities 3727,4363,H_beta(4863),4959,5007,6548,H_alpha(6565),6585,6717,6731
	# and line sigma 3727,4363,H_beta(4863),4959,5007,6548,H_alpha(6565),6585,6717,6731
	# and ew of 3727,4363,H_beta(4863),4959,5007,6548,H_alpha(6565),6585,6717,6731
	
	
	galaxy_data = Table.read(fileName, format='ascii.commented_header')
	
	# Number of galaxies
	N_gal = len(galaxy_data)
	
	
	############################################################################
	# 								OUTPUT TABLES
	############################################################################
	
	
	# Output tables
	I06_table = Table()
	pyNeb_table = Table()
	
	# Output index lists
	I06_indices = []
	pyNeb_indices = []

	# Initialize count dictionary
	C = {'I06':0, 'pyNeb':0}
	
	# Initialize total table
	tot_names = ['index', 'I06', 'pyNeb']
	tot = Table(names=tot_names, dtype=[np.int32]*len(tot_names))

	# Initialize line S/N ratio table
	StoN_rat_names = ['index', 'OII_3727', 'NeIII_3869', 'OIII_4363', 'H_BETA', 'OIII_4959', 'OIII_5007', 'NII_6548', 'H_ALPHA', 'NII_6584', 'SII_6717', 'SII_6731', 'S/N invalid']
	StoN_rat = Table(names=StoN_rat_names, dtype=[np.int32]*len(StoN_rat_names))
	
	# Initialize method matrix
	method = [['','I06','pyNeb'],
			  ['I06',0,0],
			  ['pyNeb',0,0]]

	# Initialize H_alpha flag total
	Ha = 0
	
	# Initialize flag lists
	flag_3727_list = -np.ones(N_gal)
	flag_NH_list = -np.ones(N_gal)
	flag_NeH_list = -np.ones(N_gal)
	flag_SII_list = -np.ones(N_gal)
	flag_OIII_list = -np.ones(N_gal)
	
	# Initialize ratio lists
	SII_ratio_list = -np.ones(N_gal)
	OIII_ratio_list = -np.ones(N_gal)
	
	# Initialize S/N ratio flag total dictionary
	SN = {}
	for name in StoN_rat_names[1:]:
		SN[name] = 0
	
	
	############################################################################
	#	FILTER GALAXIES
	############################################################################
	
	
	for i_galaxy in range(N_gal):

		########################################################################
		#						DETERMINE VALID METHOD(S)
		########################################################################

		# Determine to which methods this galaxy can be applied
		c, flag_3727, flag_NH, flag_NeH, StoN, flag_SII, SII_ratio, flag_OIII, OIII_ratio = count(galaxy_data[i_galaxy], sigma_4363)
		
		########################################################################
		#							RECORD RESULTS
		########################################################################
	
		# Update total number of galaxies excluded due to negative fluxes of emission lines
		for line in SN.keys():
			SN[line] = SN[line] + StoN[line]


		# Build total matrix
		tot.add_row(c)
	
		# Build line S/N ratio matrix
		StoN_rat.add_row(StoN)
		

		# Add to 3727 flag list
		flag_3727_list[i_galaxy] = flag_3727
		
		# Add to NH flag list
		flag_NH_list[i_galaxy] = flag_NH
		
		# Add to NeH flag list
		flag_NeH_list[i_galaxy] = flag_NeH
		
		# Add to SII flag list
		flag_SII_list[i_galaxy] = flag_SII
		
		# Add to OIII flag list
		flag_OIII_list[i_galaxy] = flag_OIII
		
		# Add to SII ratio list
		SII_ratio_list[i_galaxy] = SII_ratio
		
		# Add to OIII ratio list
		OIII_ratio_list[i_galaxy] = OIII_ratio

		# Just save off the index value of the galaxy for whichever method(s) it 
		# can be used for, so that we can build the tables at the end all at once.
		
		if c['I06']:
			# I06 can be used
			I06_indices.append(i_galaxy)
			
			# Increment count
			C['I06'] = C['I06'] + 1
			
		if c['pyNeb']:
			# pyNeb can be used
			pyNeb_indices.append(i_galaxy)
			
			# Increment count
			C['pyNeb'] = C['pyNeb'] + 1

		'''
		# track number of galaxies that overlap methods
		for j in range(i+1):
			if c[j]:
				method[j+1][i+1] += 1
		'''
		
	
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
	
	I06_table['index'] = galaxy_data['index'][I06_indices]
	
	for column in extra_cols:
		I06_table[column] = galaxy_data[column][I06_indices]
	
	I06_table['flag3727'] = flag_3727_list[I06_indices]
	I06_table['flagNH'] = flag_NH_list[I06_indices]
	I06_table['flagNeH'] = flag_NeH_list[I06_indices]
	I06_table['flagS2'] = flag_SII_list[I06_indices]
	I06_table['S2ratio'] = SII_ratio_list[I06_indices]
	I06_table['flagO3'] = flag_OIII_list[I06_indices]
	I06_table['O3ratio'] = OIII_ratio_list[I06_indices]
	
	for line in I06_lines:
		for ending in suffixes:
			I06_table[line + ending] = galaxy_data[line + ending][I06_indices]
			
			
	################################	pyNeb	################################
	
	# Lines needed
	pyNeb_lines = ['OII_3727', 'NeIII_3869', 'OIII_4363', 'H_BETA', 'OIII_4959', 'OIII_5007', 'NII_6548', 'H_ALPHA', 'NII_6584', 'SII_6717', 'SII_6731']
	
	pyNeb_table['index'] = galaxy_data['index'][pyNeb_indices]
	
	for column in extra_cols:
		pyNeb_table[column] = galaxy_data[column][pyNeb_indices]
	
	pyNeb_table['flag3727'] = flag_3727_list[pyNeb_indices]
	pyNeb_table['flagNH'] = flag_NH_list[pyNeb_indices]
	pyNeb_table['flagNeH'] = flag_NeH_list[pyNeb_indices]
	
	for line in pyNeb_lines:
		for ending in suffixes:
			pyNeb_table[line + ending] = galaxy_data[line + ending][pyNeb_indices]
			
			
	############################################################################
	#
	#								FIX FLAG_3727
	#
	############################################################################
	
	
	I06_table = reset_flag3727(I06_table)
	pyNeb_table = reset_flag3727(pyNeb_table)


	############################################################################
	#
	#							WRITE OUTPUT FILES
	#
	############################################################################


	I06_table.write('I06.txt', format='ascii.commented_header')
	pyNeb_table.write('pyNeb.txt', format='ascii.commented_header')
		
		
	############################################################################
	#
	#							PRINT OUTPUTS TO TERMINAL
	#
	############################################################################
	
	
	# print total count
	print 'I06:',C['I06']
	print 'pyNeb:',C['pyNeb']
	
	# print number of galaxies eliminated due to flux of [OII] 3727
	print 'Flux of [OII] 3727 < 0:',SN['OII_3727']
	
	# print number of galaxies eliminated due to flux of [NeIII] 3869
	print 'Flux of [NeIII] 3869 < 0:',SN['NeIII_3869']
	
	# print number of galaxies eliminated due to S/N of [OIII] 4363
	print 'Flux of [OIII] 4363 < 0:',SN['OIII_4363']
	
	# print number of galaxies eliminated due to S/N of H_beta
	print 'Flux of H_beta < 0:',SN['H_BETA']
	
	# print number of galaxies eliminated due to flux of [OIII] 4959
	print 'Flux of [OIII] 4959 < 0:',SN['OIII_4959']
	
	# print number of galaxies eliminated due to flux of [OIII] 5007
	print 'Flux of [OIII] 5007 < 0:',SN['OIII_5007']
	
	# print number of galaxies eliminated due to flux of [NII] 6548
	print 'Flux of [NII] 6548 < 0:',SN['NII_6548']
	
	# print number of galaxies eliminated due to flux of [NII] 6584
	print 'Flux of [NII] 6584 < 0:',SN['NII_6584']
	
	# print number of galaxies eliminated due to flux of [SII] 6716
	print 'Flux of [SII] 6717 < 0:',SN['SII_6717']
	
	# print number of galaxies eliminated due to flux of [SII] 6731
	print 'Flux of [SII] 6731 < 0:',SN['SII_6731']
	
	# print number of galaxies with non-physical SII ratios
	print '[SII] 6731/6717 < 0.68:',N_gal - sum(flag_SII_list)
	
	# print number of galaxies with non-physical OIII ratios
	print '[OIII] 4363/(4959+5007) > 0.05:',N_gal - sum(flag_OIII_list)
	

	'''
	# Print total matrix to file
	#tot_fileName = 'testData_matrix_KE_txt'
	tot_fileName = 'almostAllGalaxies_matrix_KE_MPAJHU.txt'
	tot.write(tot_fileName, format='ascii.commented_header')
	
	# print line S/N ratio matrix to file
	#StoN_rat_fileName = 'testData_SNratio_KE.txt'
	#StoN_rat_fileName = 'almostAllGalaxies_SNratio_KE_MPAJHU.txt'
	StoN_rat_fileName = 'SHELSgalaxies_SNratio_1sigI06.txt'
	StoN_rat.write(StoN_rat_fileName, format='ascii.commented_header')
	
	# print method matrix to file
	h = open('almostAllGalaxies_KEmethod_matrix_MPAJHU.txt','w')
	#h = open('testData_KEmethod_matrix.txt','w')
	for row in method:
		string = ''
		for element in row:
			string = string + str(element) + '\t'
		h.write(string+'\n')
	'''
