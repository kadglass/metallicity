'''Calculates metallicity based on pyNeb outlined in Luridiana et al. (2015)'''


################################################################################
#
#								IMPORT LIBRARIES
#
################################################################################


import pyneb as pn
import numpy as np
from astropy.table import Table
from Op_approximation import Op_approx
from ICF import N_ICF_I06, Ne_ICF_I06

import matplotlib.pyplot as plt


################################################################################
#
#					DEFINE FUNCTION FOR METALLICITY CALCULATION
#
################################################################################


def pyMet(flux):
	'''flux is a dictionary of arrays of the emission line flux values.
	flagNH is boolean - 0 if N abundance cannot be calculated, 1 if it can'''
	
	# Number of galaxies
	n_galaxies = len(flux['flag3727'])
	
	############################################################################
	#	Define ions
	############################################################################
	
	# [OIII]
	O3 = pn.Atom('O',3)
	
	# [OII]
	O2 = pn.Atom('O',2)
	
	# [NII]
	N2 = pn.Atom('N',2)
	
	# [NeIII]
	Ne3 = pn.Atom('Ne',3)
	
	'''
	# [SII]
	S2 = pn.Atom('S',2)
	'''
	
	############################################################################
	#	Calculate electron number density
	############################################################################
	
	Ne = 100.*np.ones(n_galaxies)	# cm^-3
	
	############################################################################
	#	Calculate electron temperature and density in [OIII] region
	############################################################################
	
	# (4959+5007)/4363
	O3ratio = (flux['OIII_4959_FLUX'] + flux['OIII_5007_FLUX'])/flux['OIII_4363_FLUX']
	
	'''
	# 6731/6716
	S2ratio = flux['SII_6731_FLUX']/flux['SII_6717_FLUX']
	'''
	'''for i in range(5):
		# T [K]
		T3 = O3.getTemDen(O3ratio, den=Ne, to_eval="(L(4959) + L(5007))/L(4363)", maxIter=100)
		
		# Ne [cm^-3]
		Ne = S2.getTemDen(S2ratio, tem=T3, to_eval='L(6731)/L(6716)', maxIter=100)
		
		print 'Temp:', T3, 'Density:', Ne
	'''
	'''
	# T [K], Ne [cm^-3]
	diags = pn.Diagnostics()
	T3, Ne = diags.getCrossTemDen('[OIII] 4363/5007+', '[SII] 6731/6716', 1./O3ratio, S2ratio, max_iter=10000)
	'''
	# T [K]
	T3 = O3.getTemDen(O3ratio, den=Ne, to_eval='(L(4959) + L(5007))/L(4363)', maxIter=100)
	
	# Convert T3 to numpy array if output is scalar
	if not isinstance(T3, np.ndarray):
		T3 = np.array([T3])
	
	############################################################################
	#	Calculate ion abundance of O++
	############################################################################
	
	# Abundance ratio of O++
	OppH = O3.getIonAbundance(int_ratio=(flux['OIII_4959_FLUX'] + flux['OIII_5007_FLUX']), tem=T3, den=Ne, to_eval="L(4959) + L(5007)", Hbeta=flux['H_BETA_FLUX'])
	
	# Metallicity
	Z12logOppH = 12. + np.log10(OppH)
	
	
	############################################################################
	#
	#						Calculate ion abundance of O+
	#
	############################################################################


	########################################################################
	#	Calculate electron temperature in [OII] region
	########################################################################
	
	T2 = 0.7*T3 + 3000.
	
	########################################################################
	#	Calculate ion abundance of O+
	########################################################################
	
	# Abundance ratio of O+
	OpH = O2.getIonAbundance(int_ratio=flux['OII_3727_FLUX'], tem=T2, den=Ne, to_eval="L(3726) + L(3729)", Hbeta=flux['H_BETA_FLUX'])
	
	# Approximation for O+ (if [OII] 3727 is not valid AND galaxy is SF, composite, or AGN)
	bool_3727 = np.invert(np.array(flux['flag3727'], dtype=np.bool))
	bool_SF = flux['BPTclass'] == 1
	bool_comp = flux['BPTclass'] == 3
	bool_AGN = flux['BPTclass'] == 4
	bool_approx = np.all([bool_3727, np.any([bool_SF, bool_comp, bool_AGN], axis=0)], axis=0)
	
	# Linear approximation for O+ abundance
	_,_,OpH[bool_approx] = Op_approx(flux['BPTclass'][bool_approx], flux['rabsmag'][bool_approx], Z12logOppH[bool_approx], OppH[bool_approx])
		
	
	# Metallicity
	Z12logOpH = 12. + np.log10(OpH)
	Z = 12. + np.log10(OppH + OpH)
	
	
	############################################################################
	#
	#						Calculate ion abundance of N+
	#
	############################################################################
	
	
	#if flux['flagNH']:
	# Abundance ratio of N+
	NpH = N2.getIonAbundance(int_ratio=(flux['NII_6548_FLUX'] + flux['NII_6584_FLUX']), tem=T2, den=Ne, to_eval="L(6548) + L(6584)", Hbeta=flux['H_BETA_FLUX'])
	N12logNpH = 12. + np.log10(NpH)
		
	# Nitrogen abundance
	v = OpH/(OpH + OppH)
	ICF = N_ICF_I06(Z,v)
	
	NH = ICF*NpH
	N12logNH = 12. + np.log10(NH)
	
	
	############################################################################
	#
	#						Calculate ion abundance of Ne++
	#
	############################################################################
	
	
	#if flux['flagNeH']:
	# Abundance ratio of Ne++
	NeppH = Ne3.getIonAbundance(int_ratio=flux['NeIII_3869_FLUX'], tem=T3, den=Ne, to_eval="L(3869)", Hbeta=flux['H_BETA_FLUX'])
	Ne12logNeppH = 12. + np.log10(NeppH)
		
	# Nitrogen abundance
	w = OppH/(OpH + OppH)
	ICF = Ne_ICF_I06(Z,w)
	
	NeH = ICF*NeppH
	Ne12logNeH = 12. + np.log10(NeH)
	
	############################################################################
	#	Return metallicity estimates
	############################################################################
	
	return Z, Z12logOppH, Z12logOpH, Ne, T3, N12logNH, N12logNpH, Ne12logNeH, Ne12logNeppH
	
	
################################################################################
#
#								  ERROR FUNCTION
#
################################################################################


def pyMet_error(galaxy):
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
	
	# Emission lines to generate
	em_lines = ['OII_3727_FLUX', 'NeIII_3869_FLUX', 'OIII_4363_FLUX', 'H_BETA_FLUX', 'OIII_4959_FLUX', 'OIII_5007_FLUX', 'NII_6548_FLUX', 'NII_6584_FLUX']
	
	# initialize fake flux table
	F = Table()
	F['BPTclass'] = np.ones(N)*galaxy['BPTclass']
	F['flag3727'] = np.ones(N)*galaxy['flag3727']
	
	if (galaxy['flag3727'] == 0):
		F['OII_3727_FLUX'] = galaxy['OII_3727_FLUX']*np.ones(N)
		em_lines = em_lines[1:]
	
	############################################################################
	#	Generate fluxes
	############################################################################
	
	for line in em_lines:
		if galaxy[line + '_ERR'] <= 0:
			F[line] = np.ones(N)*galaxy[line]
		else:
			F[line] = np.abs(np.random.normal(galaxy[line], galaxy[line + '_ERR'], N))
			
		'''
		F[line] = -np.ones(N)
		for i in range(N):
			# generate flux from normal distributions
			while F[line][i] <= 0.:
				F[line][i] = normalvariate(galaxy[line], galaxy[line + '_ERR'])
		'''
					
	############################################################################
	# Calculate and store metallicity
	############################################################################
	
	met,_,_,_,_,logNH,_,logNeH,_ = pyMet(F)
	
	############################################################################
	# Re-run nan metallicities
	############################################################################
	
	Z_boolarray = np.any([np.isnan(met), np.isnan(logNH), np.isnan(logNeH)], axis=0)
	
	tries = 0
	while sum(Z_boolarray.astype(int)) > 0 and tries < 500:
		# At least one estimate needs to be re-run
		for line in em_lines:
			if galaxy[line + '_ERR'] <= 0:
				F[line][Z_boolarray] = np.ones(sum(Z_boolarray.astype(int)))*galaxy[line]
			else:
				F[line][Z_boolarray] = np.abs(np.random.normal(galaxy[line], galaxy[line + '_ERR'], sum(Z_boolarray.astype(int))))
				
		# Re-run nan metallicities with new fluxes
		met[Z_boolarray],_,_,_,_,logNH[Z_boolarray],_,logNeH[Z_boolarray],_ = pyMet(F[Z_boolarray])
		Z_boolarray = np.isnan(met)
		tries += 1
		
	############################################################################
	# Determine dispersion of metallicity measurements
	############################################################################
	
	if tries == 500:
		# Failed to generate all real metallicity values
		Zerr = np.nan
		NHerr = np.nan
		NeHerr = np.nan
	else:
		Zerr = np.std(met)
		NHerr = np.std(logNH)
		NeHerr = np.std(logNeH)
		
	'''
	############################################################################
	#	Histogram of "fake" metallicities
	############################################################################
	
	plt.hist(met, bins=50)
	plt.title('Galaxy #'+str(galaxy['index'])+' (pyNeb)')
	'''
	############################################################################
	#	Return error estimate
	############################################################################

	return Zerr, NHerr, NeHerr, met
