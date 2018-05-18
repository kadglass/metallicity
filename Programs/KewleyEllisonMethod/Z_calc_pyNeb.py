'''Calculate metallicity of galaxies using pyNeb package.'''


################################################################################
#
#								IMPORT LIBRARIES
#
################################################################################


from metal_pyNeb import pyMet
from astropy.table import Table
import numpy as np
from time import time
from multiprocessing import Process, cpu_count, Queue
from math import ceil


################################################################################
#
#	METALLICITY ERROR FUNCTION (for parallelization)
#
################################################################################

def error_parallel(galaxies, Z12logOH, index_queue, results_queue, sentinel):
	
	while not index_queue.empty():
		i = index_queue.get()
		if not np.isnan(Z12logOH[i]):
			Z_err,NH_err,NeH_err = pyMet_error(galaxies[i])
			results_queue.put((i, Z_err, NH_err, NeH_err))
	
	# Return sentinel value (all done!)
	results_queue.put(sentinel)


################################################################################
#
#	METALLICITY PIPELINE FUNCTION
#
################################################################################


def metpy(outFile, extra_cols):
	'''outFile is the name of the output file
	extra_cols is the number of additional columns of data preceeding the flux 
	data in the original flux file'''
	start_time = time()
	############################################################################
	#	DATA
	############################################################################


	#fileName = '../data_comp/kias1033_5_MPAJHUvSDSSdr12_flux_comp.txt'


	############################################################################
	#	IMPORT DATA
	############################################################################
	print 'Importing data'
	t = time()


	#fluxData = Table.read(fileName, format='ascii.commented_header')
	fluxData = Table.read('pyNeb.txt', format='ascii.commented_header')


	print 'Data imported', time() - t
	############################################################################
	#							CALCULATE METALLICITY
	############################################################################
	print 'Calculating metallicity'
	t = time()


	'''flag3727 is boolean - 0 if 3727 cannot be used, 1 if it can
	flagNH is boolean - 0 if N abundance cannot be calculated, 1 if it can'''
	'''
	flux = {}

	flux['OII_3727'] = fluxData['OII_3727_FLUX']
	flux['NeIII_3869'] = fluxData['NeIII_3869_FLUX']
	flux['OIII_4363'] = fluxData['OIII_4363_FLUX']
	flux['H_BETA'] = fluxData['H_BETA_FLUX']
	flux['OIII_4959'] = fluxData['OIII_4959_FLUX']
	flux['OIII_5007'] = fluxData['OIII_5007_FLUX']
	flux['NII_6548'] = fluxData['NII_6548_FLUX']
	flux['NII_6584'] = fluxData['NII_6584_FLUX']
	flux['SII_6717'] = fluxData['SII_6717_FLUX']
	flux['SII_6731'] = fluxData['SII_6731_FLUX']
	flux['flag3727'] = fluxData['flag3727']
	flux['flagNH'] = fluxData['flagNH']
	flux['flagNeH'] = fluxData['flagNeH']
	flux['BPTclass'] = fluxData['BPTclass']
	'''
	# Calculate metallicity
	Z12logOH, Z12logOppH, Z12logOpH, Ne, T3, N12logNH, N12logNpH, Ne12logNeH, Ne12logNeppH = pyMet(fluxData)
	
	
	print 'Metallicity calculated', time() - t
	############################################################################
	#
	#								ESTIMATE ERRORS
	#
	############################################################################
	print 'Estimating errors'
	t = time()
	
	'''
	fluxError = {}
	
	fluxError['OII_3727'] = fluxData['OII_3727_FLUX_ERR']
	fluxError['NeIII_3869'] = fluxData['NeIII_3869_FLUX_ERR']
	fluxError['OIII_4363'] = fluxData['OIII_4363_FLUX_ERR']
	fluxError['H_BETA'] = fluxData['H_BETA_FLUX_ERR']
	fluxError['OIII_4959'] = fluxData['OIII_4959_FLUX_ERR']
	fluxError['OIII_5007'] = fluxData['OIII_5007_FLUX_ERR']
	fluxError['NII_6548'] = fluxData['NII_6548_FLUX_ERR']
	fluxError['NII_6584'] = fluxData['NII_6584_FLUX_ERR']
	fluxError['SII_6717'] = fluxData['SII_6717_FLUX_ERR']
	fluxError['SII_6731'] = fluxData['SII_6731_FLUX_ERR']
	'''
	Zerr = np.empty(len(fluxData))
	Zerr[:] = np.nan
	
	NHerr = np.empty(len(fluxData))
	NHerr[:] = np.nan
	
	NeHerr = np.empty(len(fluxData))
	NeHerr[:] = np.nan
	
	############################################################################
	#	INITIALIZE PARALLELIZATION
	############################################################################
	
	# Count number of logical cores available
	num_cores = cpu_count()
	
	# Number of active cores
	num_active = num_cores
	
	# Initialize process reference list
	processes = []
	
	# Initialize results queue
	results_queue = Queue()
	
	# Initialize index queue
	index_queue = Queue()
	for i in xrange(len(fluxData)):
		index_queue.put(i)
	
	# Process termination sentinel
	sentinel = 'poop'
	
	# Number of runs per process
	runs = int(ceil(len(fluxData)/num_cores))
	
	# Initiate processes
	for _ in range(num_cores):
		P = Process(target=error_parallel, args=(fluxData, Z12logOH, index_queue, results_queue, sentinel))
		processes.append(P)
		
		P.start()
	
	############################################################################
	############################################################################
	
	while num_active > 0:
		
		try:
			# Extract next available output from the queue
			result_value = results_queue.get()
		except Queue.Empty:
			# Queue is empty (either all processes are finished or all are calculating)
			continue
		
		if result_value == sentinel:
			# Process has finished
			num_active -= 1
		else:
			# Unpack output tuple
			i, Z_err, NH_err, NeH_err = result_value
			
			# Save values to arrays
			Zerr[i] = Z_err
			NHerr[i] = NH_err
			NeHerr[i] = NeH_err
			
			
	############################################################################
	############################################################################
	
	# Join processes
	for p in processes:
		p.join(None)
	
	
	'''
	for i in range(len(fluxData)):
		# Calculate metallicity error if metallicity is not nan
		if not np.isnan(Z12logOH[i]):
			Zerr[i],NHerr[i],NeHerr[i] = pyMet_error(fluxData[i])
	'''

	print 'Errors calculated:', time() - t
	################################################################################
	#
	#									OUTPUT
	#
	################################################################################
	print 'Saving output'
	t = time()


	Zdata = Table()

	Zdata['index'] = fluxData['index']

	for column in extra_cols:
		Zdata[column] = fluxData[column]

	Zdata['Z12logOH'] = Z12logOH
	Zdata['Zerr'] = Zerr
	Zdata['Z12logOppH'] = Z12logOppH
	Zdata['Z12logOpH'] = Z12logOpH
	Zdata['N12logNH'] = N12logNH
	Zdata['NHerr'] = NHerr
	Zdata['N12logNpH'] = N12logNpH
	Zdata['Ne12logNeH'] = Ne12logNeH
	Zdata['NeHerr'] = NeHerr
	Zdata['Ne12logNeppH'] = Ne12logNeppH
	Zdata['T3'] = T3
	Zdata['Ne'] = Ne
	Zdata['flag3727'] = fluxData['flag3727']


	################################################################################
	#
	#									SAVE OUTPUT
	#
	################################################################################
	
	
	Zdata.write(outFile, format='ascii.commented_header')


	print 'Output saved', time() - t
	print 'Total run time', time() - start_time
