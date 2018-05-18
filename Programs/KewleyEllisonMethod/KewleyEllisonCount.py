'''Filter through database to determine if all line fluxes are present for 
various metallicity calculation methods.

10 methods listed from Kewley & Ellison (2008):
	0: M91: [OII]3727, H_b, [OIII]4959,5007, [NII]6584
	1: KD02: [OII]3727, H_b, [OIII]4959,5007, [NII]6584
	2: KK04: [OII]3727, H_b, [OIII]4959,5007,  [NII]6584
	3: Z94: [OII]3727, H_b, [OIII]4959,5007
	4: T04: [OII]3727, H_b, [OIII]4959,5007, [NII]6548,6584, H_a, [SII]6716,6731
0:	5: I06: ([OII]3727,) [OIII]4363,4959,5007, H_b, ([SII]6716,6731)
1:	6: PP04: (H_b, [OIII]5007,) H_a, [NII]6584
	7: P05: [OII]3727, H_b, [OIII]4959,5007, [NII]6584
2:	8: D02: H_a, [NII]6584

3:  pyNeb: [OII3727], [OIII]4363,4959,5007, H_b, [SII]6716,6731
'''


################################################################################
#
#								DEFINE FUNCTION
#
################################################################################


def count(galaxy, sigma_4363):
	'''Based on the provided line fluxes, determines which methods can be used 
	to calculate the galaxy's metallicity.
	Conditions: must be emission line, and H_beta must have at least a 5-sigma 
	detection.
	Addendum: [OIII] 4363 must only have a sigma_4363 detection.'''
	
	
	############################################################################
	#
	#							INITIALIZE OUTPUTS
	#
	############################################################################
	
	
	# Flag for [OII] 3727: 1 - pass, 0 - fail
	flag_3727 = 1

	# Flag for [SII] 6717 or [SII] 6731: 1 - pass, 0 - fail
	flag_OppSp = 1
	
	# Flag for N/H ratio: 1 - pass, 0 - fail
	flag_NH = 1
	
	# Flag for Ne/H ratio: 1 - pass, 0 - fail
	flag_NeH = 1
	
	# Flag dictionary for signal-to-noise of lines: 1 - fail, 0 - pass
	StoN = {'index':galaxy['index'], 'S/N invalid':0}
	
	# Flag for [SII] 6731/6717 validity: 1 - pass, 0 - fail
	flag_SII = 0
	
	# Flag for [OIII] 4363/(4959+5007) validity: 1 - pass, 0 - fail
	flag_OIII = 0
	
	# Flag dictionary for methods: 1 - can use, 0 - cannot use
	C = {'index':galaxy['index'], 'I06':1, 'PP04':1, 'D02':1, 'pyNeb':1}
	
	
	############################################################################
	#
	#							EVALUATE FLUX QUALITY
	#
	############################################################################
	# Restrictions: H_beta > 5sigma, [OIII] 4363 > 1sigma
	# T04 restrictions: H_beta > 5sigma, H_alpha > 5sigma, [NII] 6584 > 5sigma
	# I06 restrictions: [OIII] 4363 > 1sigma, H_beta > 10^-14 erg/s/cm^2, [OIII] 4959/H_beta > 0.7 and [OII] 3727/H_beta < 1.0

	############################################################################
	# [OII] 3727
	############################################################################
	
	if (galaxy['OII_3727_FLUX'] <= 0.):
		flag_3727 = 0
		StoN['OII_3727'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['OII_3727'] = 0
		
	############################################################################
	# [NeIII] 3869
	############################################################################
	
	if (galaxy['NeIII_3869_FLUX'] <= 0.):
		flag_NeH = 0
		StoN['NeIII_3869'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['NeIII_3869'] = 0
		
	############################################################################
	# [OIII] 4363
	############################################################################
	
	if (galaxy['OIII_4363_FLUX'] <= 0.) or (galaxy['OIII_4363_FLUX'] < sigma_4363*galaxy['OIII_4363_FLUX_ERR']):
		# cannot use I06, pyNeb
		C['I06'] = 0
		C['pyNeb'] = 0
		StoN['OIII_4363'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['OIII_4363'] = 0
		
	############################################################################
	# H_beta
	############################################################################
	
	if (galaxy['H_BETA_FLUX'] <= 0.) or (galaxy['H_BETA_FLUX'] < 5.*galaxy['H_BETA_FLUX_ERR']):
		# cannot use I06, pyNeb
		C['I06'] = 0
		C['pyNeb'] = 0
		StoN['H_BETA'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['H_BETA'] = 0
			
	############################################################################
	# [OIII] 4959
	############################################################################
	
	if (galaxy['OIII_4959_FLUX'] <= 0.):
		# cannot use I06, pyNeb
		C['I06'] = 0
		C['pyNeb'] = 0
		StoN['OIII_4959'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['OIII_4959'] = 0
			
	############################################################################
	# [OIII] 5007
	############################################################################
	
	if (galaxy['OIII_5007_FLUX'] <= 0.):
		# cannot use I06, pyNeb
		C['I06'] = 0
		C['pyNeb'] = 0
		StoN['OIII_5007'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['OIII_5007'] = 0
		
	############################################################################
	# [NII] 6548
	############################################################################
	
	if (galaxy['NII_6548_FLUX'] <= 0.):
		flag_NH = 0
		StoN['NII_6548'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['NII_6548'] = 0
		
	############################################################################
	# H_alpha
	############################################################################
	
	if (galaxy['H_ALPHA_FLUX'] <= 0.):# or (galaxy['H_ALPHA_FLUX'] < 5.*galaxy['H_ALPHA_FLUX_ERR']):
		# cannot use PP04, D02, (I06)
		C['PP04'] = 0
		C['D02'] = 0
		StoN['H_ALPHA'] = 1
	else:
		StoN['H_ALPHA'] = 0
		
	############################################################################
	# [NII] 6584
	############################################################################
	
	if (galaxy['NII_6584_FLUX'] <= 0.):# or (galaxy['NII_6584_FLUX'] < 5.*galaxy['NII_6584_FLUX_ERR']):
		# Cannot use PP04, D02, (I06)
		C['PP04'] = 0
		C['D02'] = 0
		flag_NH = 0
		StoN['NII_6584'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['NII_6584'] = 0
		
	############################################################################
	# [SII] 6716
	############################################################################
	
	if (galaxy['SII_6717_FLUX'] <= 0.):
		# Cannot use O+ approximation
		flag_OppSp = 0
		StoN['SII_6717'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['SII_6717'] = 0
		
	############################################################################
	# [SII] 6731
	############################################################################
	
	if (galaxy['SII_6731_FLUX'] <= 0.):
		# Cannot use O+ approximation
		flag_OppSp = 0
		StoN['SII_6731'] = 1
		StoN['S/N invalid'] = 1
	else:
		StoN['SII_6731'] = 0
		
	############################################################################
	# [SII] 6731/6717
	############################################################################
	
	if (StoN['SII_6731'] == 0) and (StoN['SII_6717'] == 0):
	
		SIIratio = galaxy['SII_6731_FLUX']/galaxy['SII_6717_FLUX']
	
		if (SIIratio > 0.68):
			# Not physically possible
			flag_SII = 1
	
	else:
		SIIratio = 0.
		
	############################################################################
	# [OIII] 4363/(4959+5007)
	############################################################################
	
	if (StoN['OIII_4363'] == 0) and (StoN['OIII_4959'] == 0) and (StoN['OIII_5007'] == 0):
	
		OIIIratio = galaxy['OIII_4363_FLUX']/(galaxy['OIII_4959_FLUX'] + galaxy['OIII_5007_FLUX'])
	
		if (OIIIratio < 0.05):
			# Not physically possible outside this range
			flag_OIII = 1
	
	else:
		OIIIratio = 0.

	
	'''	
	############################################################################
	# I06 RESTRICTIONS
	############################################################################
	# Only if it has made it through everything else
	
	if C['I06']:
		if (galaxy['OIII_4959_FLUX']/galaxy['H_BETA_FLUX'] < 0.7) and (galaxy['OII_3727_FLUX']/galaxy['H_BETA_FLUX'] > 1.):
			C['I06'] = 0
			C['pyNeb'] = 0
	'''	
	
	
	############################################################################
	#
	#								RETURN OUTPUT
	#
	############################################################################
	
	
	return C, flag_3727, flag_OppSp, flag_NH, flag_NeH, StoN, flag_SII, SIIratio, flag_OIII, OIIIratio
