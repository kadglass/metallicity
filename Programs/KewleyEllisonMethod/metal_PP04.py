'''Calculates metallicity based on N2 method outlined in Pettini and Pagel (2004)'''

import numpy as np

def PP04(flux,flag):
	'''flux[0] is H_alpha (6563)
	flux[1] is [NII] 6583'''
	
	# N2 ratio
	N2 = np.log10(flux[1]/flux[0])
	
	# Calculate metallicity
	logOH12 = 9.37 + 2.03*N2 + 1.26*N2**2 + 0.32*N2**3
	'''
	# Calibrate to T04 (coefficients from KE08)
	a = -1661.9380
	b = 585.17650
	c = -68.471750
	d = 2.6766690
	y = a + b*logOH12 + c*logOH12**2 + d*logOH12**3
	'''
	return logOH12
