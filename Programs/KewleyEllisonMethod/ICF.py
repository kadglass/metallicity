################################################################################
#
#	IMPORT MODULES
#
################################################################################


import numpy as np

import warnings
#warnings.filterwarnings('ignore', 'RuntimeWarning')


################################################################################
#
#	DEFINE FUNCTIONS
#
################################################################################


def N_ICF_I06(Z,v):
	'''Calculates the ICF for N based on that published by Izotov06'''
	
	# Number of galaxies
	N = len(Z)
	
	# Initialize ICF arrays
	ICF = np.ones(N)
	
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		############################################################################
		# Low metallicities
		boolarray = Z <= 7.2
	
		# Check to make sure there is at least one galaxy within this regime
		if np.any(boolarray):
			ICF[boolarray] = -0.825*v[boolarray] + 0.718 + 0.853/v[boolarray]
	
		############################################################################
		# 7.2 < Z < 7.6
		boolarray = np.logical_and(Z > 7.2, Z <= 7.6)
	
		ICF_low = -0.825*v + 0.718 + 0.853/v
		ICF_int = -0.809*v + 0.712 + 0.852/v
	
		# Check to make sure there is at least one galaxy within this regime
		if np.any(boolarray):
			# linear interpolate to find value
			ICF[boolarray] = ICF_low[boolarray] + (ICF_int[boolarray] - ICF_low[boolarray])*((Z[boolarray] - 7.2)/(7.6 - 7.2))
	
		############################################################################
		# 7.6 < Z < 8.2
		boolarray = np.logical_and(Z > 7.6, Z < 8.2)
	
		ICF_int = -0.809*v + 0.712 + 0.852/v
		ICF_high = -1.476*v + 1.752 + 0.688/v
	
		# Check to make sure there is at least one galaxy within this regime
		if np.any(boolarray):
			# linear interpolate to find value
			ICF[boolarray] = ICF_int[boolarray] + (ICF_high[boolarray] - ICF_int[boolarray])*((Z[boolarray] - 7.6)/(8.2 - 7.6))
	

		############################################################################
		# High metallicities
		boolarray = Z >= 8.2
	
		# Check to make sure there is at least one galaxy within this regime
		if np.any(boolarray):
			ICF[boolarray] = -1.476*v[boolarray] + 1.752 + 0.688/v[boolarray]
	
	return ICF
	
	
def Ne_ICF_I06(Z,w):
	'''Calculates the ICF for Ne based on that published by Izotov06'''
	
	# Number of galaxies
	N = len(Z)
	
	# Initialize ICF arrays
	ICF = np.ones(N)
	
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		############################################################################
		# Low metallicities
		boolarray = Z <= 7.2
	
		# Check to make sure there is at least one galaxy within this regime
		if np.any(boolarray):
			ICF[boolarray] = -0.385*w[boolarray] + 1.365 + 0.022/w[boolarray]
	
		############################################################################
		# 7.2 < Z < 7.6
		boolarray = np.logical_and(Z > 7.2, Z <= 7.6)
	
		ICF_low = -0.385*w + 1.365 + 0.022/w
		ICF_int = -0.405*w + 1.382 + 0.021/w
	
		# Check to make sure there is at least one galaxy within this regime
		if np.any(boolarray):
			# linear interpolate to find value
			ICF[boolarray] = ICF_low[boolarray] + (ICF_int[boolarray] - ICF_low[boolarray])*((Z[boolarray] - 7.2)/(7.6 - 7.2))
	
		############################################################################
		# 7.6 < Z < 8.2
		boolarray = np.logical_and(Z > 7.6, Z < 8.2)
	
		ICF_int = -0.405*w + 1.382 + 0.021/w
		ICF_high = -0.591*w + 0.927 + 0.546/w
	
		# Check to make sure there is at least one galaxy within this regime
		if np.any(boolarray):
			# linear interpolate to find value
			ICF[boolarray] = ICF_int[boolarray] + (ICF_high[boolarray] - ICF_int[boolarray])*((Z[boolarray] - 7.6)/(8.2 - 7.6))
	

		############################################################################
		# High metallicities
		boolarray = Z >= 8.2
	
		# Check to make sure there is at least one galaxy within this regime
		if np.any(boolarray):
			ICF[boolarray] = -0.591*w[boolarray] + 0.927 + 0.546/w[boolarray]
	
	return ICF