'''Calculates metallicity based on method outlined in Denicolo, Terlevich, and
Terlevich (2002)'''

import numpy as np

def D02(flux,flag):
	'''flux[0] is H_alpha (6563)
	flux[1] is [NII] 6584'''
	
	# N2 ratio
	N2 = np.log10(flux[1]/flux[0])
	
	# Calculate metallicity
	logOH12 = 9.12 + 0.73*N2
	'''
	# Calibrate to T04 (coefficients from KE08)
	a = 193.9000
	b = -64.87895
	c = 7.411102
	d = -0.2756653
	y = a + b*logOH12 + c*logOH12**2 + d*logOH12**3
	'''
	return logOH12
