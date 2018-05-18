'''Set flag3727 to 0 for all galaxies who have invalid 3727 flux values.'''


################################################################################
#
#								IMPORT LIBRARIES
#
################################################################################


from numpy import ndarray
from astropy.table import Table


def reset_flag3727(data_table):


	################################################################################
	#
	#								USER INPUTS
	#
	################################################################################


	#data_file = '../Programs/KewleyEllisonMethod/comp_Z_KE_I06_2sig_kias1033_5_KE_MPAJHU_flux_oii.txt'
	#data_file = '../Programs/KewleyEllisonMethod/comp_Z_KE_I06_2sig_scaled_dwarf_voidstatus_KE_MPAJHU_flux_oii_NO.txt'
	#data_file = '../Programs/KewleyEllisonMethod/comp_Z_KE_I06_1sig_scaled_I06limits_dwarf_voidstatus_KE_MPAJHU_flux_oii_NO.txt'
	#data_file = '../Programs/KewleyEllisonMethod/comp_Z_KE_I06_1sigI06_dwarf_voidstatus_stellarMass_KE_MPAJHU_flux_oii_NO.txt'
	#data_file = '../Programs/KewleyEllisonMethod/comp_Z_KE_I06_0sig_scaled_I06limits_dwarf_voidstatus_KE_MPAJHU_flux_oii_NO.txt'
	#data_file = '../Programs/KewleyEllisonMethod/comp_Z_KE_I06_1sigI06_kias1033_5_KE_MPAJHU_flux_oii_NO.txt'

	bad_gal_file = '/Users/kellydouglass/Documents/Drexel/Research/Data/Invalid_oii_flux_galaxies.txt'


	################################################################################
	#
	#								IMPORT FILES
	#
	################################################################################


	# Data file
	#data = Table.read(data_file, format='ascii.commented_header')

	# List of galaxies with invalid oii_flux values
	bad_gal = Table.read(bad_gal_file, format='ascii.no_header')


	################################################################################
	#
	#								RESET FLAG3727
	#
	################################################################################


	for index in bad_gal['col1']:

		# Find galaxy in data table
		bool_array = data_table['index'] == index
	
		# Set flag3727 to 0
		data_table['flag3727'][bool_array] = 0
	
		#data[bool_array].pprint()
	
	
	'''################################################################################
	#
	#								SAVE DATA
	#
	################################################################################


	data.write(data_file, format='ascii.commented_header')
	'''
	
	
	############################################################################
	#
	#								RETURN DATA
	#
	############################################################################
	
	
	return data_table