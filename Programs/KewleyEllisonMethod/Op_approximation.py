'''Calculate the O+ abundance based on linear relationship between O++ and O/H'''


################################################################################
#
#	IMPORT MODULES
#
################################################################################


import numpy as np

#import warnings


################################################################################
################################################################################


def Op_ratio_approx(BPTclass, rabsmag, logOppHp12, OppHp, OppSp):
	'''Includes the dependence of the ionization parameter as approximated by 
	O++/S+'''

	# Number of galaxies
	n_galaxies = len(BPTclass)

	############################################################################
	# Linear coefficients (depends on BPT classification)
	############################################################################
	
	# Initialize coefficient arrays
	a1 = np.zeros(n_galaxies)
	a2 = np.zeros(n_galaxies)
	a3 = np.zeros(n_galaxies)
	a4 = np.zeros(n_galaxies)
	b = np.zeros(n_galaxies)
	
	# Star-forming galaxies
	bool_SF = BPTclass == 1

	# dwarf galaxies
	bool_dwarf = rabsmag >= -17.

	# SF dwarf galaxies
	bool_SFdwarf = np.all([bool_SF, bool_dwarf], axis=0)
	a1[bool_SFdwarf] = 0.8148
	a2[bool_SFdwarf] = -0.2715
	a3[bool_SFdwarf] = 0.03064
	a4[bool_SFdwarf] = 0.00056
	b[bool_SFdwarf] = 1.927


	# SF non-dwarf galaxies
	bool_SFbright = np.all([bool_SF, np.logical_not(bool_dwarf)], axis=0)
	a1[bool_SFbright] = 0.7786
	a2[bool_SFbright] = -0.402
	a3[bool_SFbright] = 0.04511
	a4[bool_SFbright] = 0.0008964
	b[bool_SFbright] = 2.292

	'''
	# NOT UPDATED VALUES FOR INCLUSION WITH O++/S+
	
	# AGN or composite galaxies
	bool_comp = BPTclass == 3
	bool_AGN = BPTclass == 4
	bool_compAGN = np.all([bool_comp, bool_AGN], axis=0)
	a[bool_compAGN] = 0.7643
	b[bool_compAGN] = 2.053
	'''

	############################################################################
	# Full metallicity
	############################################################################
	
	Z = b + a1*logOppHp12 + a2*OppSp + a3*logOppHp12*OppSp + a4*OppSp**2

	OH = 10**(Z - 12.)

	############################################################################
	# O+ abundance
	############################################################################
	
	OpHp = OH - OppHp
	
	# Initialize
	logOpHp12 = np.nan*np.ones(n_galaxies)
	
	bool_abund = OpHp > 0
	# If OpHp < 0, must return nan values for all abundances
	Z[np.invert(bool_abund)] = np.nan
	OpHp[np.invert(bool_abund)] = np.nan
	
		
	return Z,logOpHp12,OpHp






def Op_approx(BPTclass, rabsmag, logOppHp12, OppHp):

	# Number of galaxies
	n_galaxies = len(BPTclass)

	############################################################################
	# Linear coefficients (depends on BPT classification)
	############################################################################
	
	# Initialize coefficient arrays
	a = np.zeros(n_galaxies)
	b = np.zeros(n_galaxies)
	
	# Star-forming galaxies
	bool_SF = BPTclass == 1
	# dwarf galaxies
	bool_dwarf = rabsmag >= -17.
	# SF dwarf galaxies
	bool_SFdwarf = np.all([bool_SF, bool_dwarf], axis=0)
	a[bool_SFdwarf] = 0.8366
	b[bool_SFdwarf] = 1.608
	# SF non-dwarf galaxies
	bool_SFbright = np.all([bool_SF, np.logical_not(bool_dwarf)], axis=0)
	a[bool_SFbright] = 0.8238
	b[bool_SFbright] = 1.787
	
	# AGN or composite galaxies
	bool_comp = BPTclass == 3
	bool_AGN = BPTclass == 4
	bool_compAGN = np.all([bool_comp, bool_AGN], axis=0)
	a[bool_compAGN] = 0.7643
	b[bool_compAGN] = 2.053

	############################################################################
	# Full metallicity
	############################################################################
	
	Z = a*logOppHp12 + b

	OH = 10**(Z - 12.)

	############################################################################
	# O+ abundance
	############################################################################
	
	OpHp = OH - OppHp
	
	# Initialize
	logOpHp12 = np.nan*np.ones(n_galaxies)
	
	bool_abund = OpHp > 0
	logOpHp12[bool_abund] = 12. + np.log10(OpHp[bool_abund])
	# If OpHp < 0, must return nan values for all abundances
	Z[np.invert(bool_abund)] = np.nan
	OpHp[np.invert(bool_abund)] = np.nan
	
		
	return Z,logOpHp12,OpHp