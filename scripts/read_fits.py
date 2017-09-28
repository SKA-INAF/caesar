#!/usr/bin/python

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

## Graphics modules
import matplotlib.pyplot as plt
import datetime
import numpy as np

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## ASTRO
from astropy.io import fits
import scipy.stats as stats
#from scipy.stats import kurtosis, skew
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats 
##################################################



#### GET SCRIPT ARGS ####
def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")
	parser.add_argument('-i', '--file', dest='filename', required=True, type=str, action='store',help='file name')
	
	args = parser.parse_args()	

	return args

def compute_stats(x):
	""" Compute stats """
	
	## Compute stats
	npixels= np.size(x)
	pixel_min= np.min(x)
	pixel_max= np.max(x)
	mean= np.mean(x)
	stddev= np.std(x,ddof=1)
	median= np.median(x)
	mad= mad_std(x)
	skewness= stats.skew(x)
	kurtosis= stats.kurtosis(x)
	
	## Compute robust stats
	niter= 1
	sigmaclip= 3
	[mean_clipped, median_clipped, stddev_clipped] = sigma_clipped_stats(x, sigma=sigmaclip, iters=niter, std_ddof=1)

	print '*** IMG STATS ***'	
	print 'n=',npixels
	print 'min/max=',pixel_min,'/',pixel_max
	print 'mean=',mean
	print 'stddev=',stddev
	print 'median=',median
	print 'mad=',mad
	print 'skew=',skewness
	print 'kurtosis=',kurtosis
	print 'mean_clipped=',mean_clipped
	print 'median_clipped=',median_clipped
	print 'stddev_clipped=',stddev_clipped
	print '*****************'
	
	
##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	## Get script args
	print('INFO: Get script args')
	try:
		args= get_args()
	except Exception as ex:
		print("ERROR: Failed to get and parse options (err=%s)",str(ex))
		return 1

	filename= args.filename

	print("*** ARGS ***")
	print("filename: %s" % filename)
	print("************")
		

	## Read FITS	
	print ('Reading FITS image %s ...' % filename)
	start_read = time.time()
	hdu= fits.open(filename, memmap=False)
	img_data = hdu[0].data
	end_read = time.time()
	elapsed_read = end_read - start_read

	## Get pixel list without zeros & nan
	#print ('Getting pixel list ...')
	#x= np.ravel(img_data)
	#print(x)

	## Compute stats
	print ('INFO: Computing stats...')
	start_stats = time.time()
	x = img_data[np.logical_and(np.isfinite(img_data),img_data!=0)]
	compute_stats(x)
	end_stats = time.time()
	elapsed_stats = end_stats - start_stats

	## Print performance info
	print ('*** PERFORMANCE INFO***')
	print ('t_read (s): %s' % str(elapsed_read))
	print ('t_stats (s): %s' % str(elapsed_stats))

			
###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

