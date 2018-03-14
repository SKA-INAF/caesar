#!/usr/bin/env python

##################################################
###          MODULE IMPORT
##################################################
## STANDARD MODULES
import os
import sys
import subprocess
import string
import time
import signal
from threading import Thread
import datetime
import numpy as np
import random
import math
import errno

## ASTRO MODULES
from scipy import constants
from astropy.io import fits

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

def str2bool(v):
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")
	
	# - MANDATORY OPTIONS
	parser.add_argument('-vis', '--vis', dest='vis', required=True, type=str,action='store',help='Input visibility CASA table name to be deconvolved')
	
	# OPTIONAL OPTIONS	
	parser.add_argument('-outvis', '--outvis', dest='outvis', required=False, type=str, default='vis_deconv.ms',action='store',help='Output visibility CASA table name stored after deconvolution')
	parser.add_argument('-bmaj', '--bmaj', dest='bmaj', required=True, type=float, default=10, action='store',help='Beam bmaj in arcsec (default=5)')
	parser.add_argument('-bmin', '--bmin', dest='bmin', required=True, type=float, default=5, action='store',help='Beam bmin in arcsec (default=5)')
	parser.add_argument('-uvdist_min', '--uvdist_min', dest='uvdist_min', required=False, type=float, default='0',action='store',help='uvdist min value used in flagging (default=0)')
	parser.add_argument('-uvdist_max', '--uvdist_max', dest='uvdist_max', required=False, type=float, default='1000000',action='store',help='uvdist max value used in flagging (default=0)')
	parser.add_argument('-flagdata','--flagdata', dest='flagdata', action='store_true')	
	parser.set_defaults(flagdata=False)
	parser.add_argument('-c', dest='scriptname', required=False, type=str, default='',action='store',help='Script name')

	args = parser.parse_args()	

	return args

def sigma2fwhm():
	f= 2.*np.sqrt(2*np.log(2.))
	return f

def gauss(x,sigma,A=1,mu=0):
	return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def get_deconv_gaus_sigma(bmaj,bmin,freq):
	""" Get gaussian sigma of deconvolution gaussian (in fourier plane)"""
	""" bmaj/bmin in arcsec """
	""" freq in GHz """
	c= constants.c
	f= sigma2fwhm()
	sigmaX= 91000./bmaj*c/freq*2/f
	sigmaY= 91000./bmin*c/freq*2/f
	return (sigmaX,sigmaY)


def deconvolve(vis,visout,bmaj,bmin):
	""" Deconvolve visibility """
	
	## Create a new visibility set and work with that
	print('INFO: Copying input visibility set in %s for later modification...' % str(visout))
	concat(vis=vis,concatvis=visout)

	###################################
	## Open SPECTRAL_WINDOW table and compute sigmaX/sigmaY in Fourier plane for all frequency channels
	###################################
	print('INFO: Opening SPECTRAL_WINDOW table and computing sigmaX/sigmaY in Fourier plane for all frequency channels present...')
	tb.open(visout + '/SPECTRAL_WINDOW')
	freq_channels= tb.getcol('CHAN_FREQ')
	sigmaU_list= []
	sigmaV_list= []
	for freq in freq_channels:
		print('INFO: Computing gaus sigma for frequency %s ...' % str(freq) )
		(sigmaU,sigmaV)= get_deconv_gaus_sigma(bmaj,bmin,freq)
		sigmaU_list.append(sigmaU)
		sigmaV_list.append(sigmaV)
	
	print('INFO: sigmaU')
	print(sigmaU_list)	
	print('INFO: sigmaV')
	print(sigmaV_list)

	nchannels= len(freq_channels) 
	print('INFO: %s channels present' % str(nchannels) )
	if nchannels <= 0:
		print('ERROR: Empty number of channels (check table reading or data integrity!)')
		return -1

	tb.close()
	
	#############################
	## Open vis MAIN table
	#############################
	print ('Opening vis file %s' % str(vis))
	tb.open(visout,nomodify=False)

	## Compute uvdistance
	uvdist_list= np.sqrt(tb.getcol('UVW')[0]**2+tb.getcol('UVW')[1]**2+tb.getcol('UVW')[2]**2)

	## Compute corrected u & v
	data= tb.getcol('DATA')
	data_corr= data.copy()
	for freq_index in range(data.shape[1]):
		for uv_index in range(data.shape[2]):
			ucorr= gauss(uvdist_list[uv_index],sigmaU_list[freq_index])
			vcorr= gauss(uvdist_list[uv_index],sigmaV_list[freq_index])
			data_corr[0,freq_index,uv_index]/= ucorr[0]
			data_corr[1,freq_index,uv_index]/= vcorr[0]
			#data_corr[0,freq_index,uv_index]= data[0,freq_index,uv_index]/ucorr[0]
			#data_corr[1,freq_index,uv_index]= data[1,freq_index,uv_index]/vcorr[0]

	## Write corrected u & v to table
	tb.putcol('CORRECTED_DATA',data_corr)

	## Close table
	tb.close()

	return 0

	
def flag(vis,uvdist_min,uvdist_max):
	""" Flag visibility data by uvdist range """ 
	uvrange= str(uvdist_min) + '~' + str(uvdist_max)
	flagdata(vis=vis,uvrange=uvrange)


##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	#===========================
	#==   Get script args
	#===========================
	print('INFO: Get script args')
	try:
		args= get_args()
	except Exception as ex:
		print("Failed to get and parse options (err=%s)",str(ex))
		return 1

	vis= args.vis
	outvis= args.outvis
	Bmaj= args.bmaj
	Bmin= args.bmin
	flagdata= args.flagdata	
	uvdist_min= args.uvdist_min
	uvdist_max= args.uvdist_max
	
	print("*** ARGS ***")
	print("vis: %s" % vis)
	print("outvis: %s" % outvis)
	print("Beam (Bmaj/Bmin): (%s,%s)" % (Bmaj, Bmin))
	print("flag data? %s uvdist min/max=%s/%s" % (flagdata,uvdist_min,uvdist_max) )
	print("************")

	#===========================
	#==   Deconvolve vis
	#===========================
	print('INFO: Running visibility deconvolution...')
	deconvolve(vis=vis,visout=outvis,bmaj=Bmaj,bmin=Bmin)

	#===========================
	#==   Flagging data
	#===========================
	if flagdata:
		print('INFO: Flagging deconvolvedibilities by uvdist range...')
		flag(outvis,uvdist_min,uvdist_max)


###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	#main()
	sys.exit(main())

