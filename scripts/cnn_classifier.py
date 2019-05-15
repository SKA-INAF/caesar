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

## ASTRO
from astropy.io import fits

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections


#### GET SCRIPT ARGS ####
def str2bool(v):
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

###########################
##     ARGS
###########################
def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")

	## INPUT DATA
	parser.add_argument('-filelist_bkg', '--filelist_bkg', dest='filelist_bkg', required=True, type=str,action='store',help='Filename with list of train bkg data')
	parser.add_argument('-filelist_source', '--filelist_source', dest='filelist_source', required=True, type=str,action='store',help='Filename with list of train source data')
	parser.add_argument('-nx', '--nx', dest='nx', required=False, type=int, default=101, action='store',help='Image width in pixels (default=101)')
	parser.add_argument('-ny', '--ny', dest='ny', required=False, type=int, default=101, action='store',help='Image height in pixels (default=101)')	
	parser.add_argument('-nsamples', '--nsamples', dest='nsamples', required=False, type=int, default=10, action='store',help='Number of train images extracted from input maps (default=10)')	
	parser.add_argument('--saveimg', dest='saveimg', action='store_true')	
	parser.set_defaults(saveimg=False)

	args = parser.parse_args()	

	return args

###########################
##     READ INPUT DATA
###########################
# Read filelist
def read_list(filename):
	""" Read a file list line by line """
	
	try:
		f = open(filename, 'r')
	except IOError:
		errmsg= 'Could not read file: ' + filename
		print "ERROR: " + errmsg
		raise IOError(errmsg)

	fields= []
	for line in f:
		line = line.strip()
		line_fields = line.split()
		fields.append(line_fields)

	f.close()	

	return fields




# Read bkg data
def read_bkg_data(filename):
	""" Read input bkg data """

	# - Read list with files
	filelist_data= []
	try:
		filelist_data= read_list(filename)
	except IOError:
		errmsg= 'Cannot read file: ' + filename
		print "ERROR: " + errmsg
		raise IOError(errmsg)

	filelist_data_shape= np.shape(filelist_data)
	print 'filelist_data=',filelist_data	

	# - Check list
	nfiles_img= filelist_data_shape[0]
	if nfiles_img<=0:
		errmsg= 'Empty file: ' + filename
		print "ERROR: " + errmsg
		raise IOError(errmsg)
	print ('INFO: #%d bkg image data found in list...' % nfiles_img)

	# - Read data in list
	#for item in filelist_data:
	#	filename= item[0]
	#	print ('INFO: Reading file %s ...' % filename) 

	#	# - Read bkg img
	#	try:
	#		data= read_img(filename)
	#	except Exception as ex:
	#		errmsg= 'Failed to read image data (err=' + str(ex) + ')'
	#		print "ERROR: " + errmsg
	#		raise IOError(errmsg)
		
	return filelist_data
	

# Read source data
def read_source_data(filename):
	""" Read input source data """
	filelist_data= []
	try:
		filelist_data= read_list(filename)
	except IOError:
		errmsg= 'Cannot read file: ' + filename
		print "ERROR: " + errmsg
		raise IOError(errmsg)

	filelist_data_shape= np.shape(filelist_data)

	print 'filelist_data=',filelist_data		

	# - Check list
	if len(filelist_data_shape) != 2:
		errmsg= 'Invalid number of columns in file ' + filename + ' (2 expected)!'
		print "ERROR: " + errmsg
		raise IOError(errmsg)

	nfiles_img= filelist_data_shape[0]
	if nfiles_img<=0:
		errmsg= 'Empty file: ' + filename
		print "ERROR: " + errmsg
		raise IOError(errmsg)
	print ('INFO: #%d source image data found in list...' % nfiles_img)

	# - Read data in list
	for item in filelist_data:
		filename= item[0]
		filename_spar= item[1]
		print ('INFO: Reading files: %s, %s ...' % (filename,filename_spar) ) 


###########################
##     PREPARE TRAIN DATA
###########################
# - Write FITS image
def write_fits(data,filename):
	""" Read data to FITS image """
	hdu= fits.PrimaryHDU(data)
	hdul= fits.HDUList([hdu])
	hdul.writeto(filename,overwrite=True)
	
# - Read FITS image
def read_img(filename):
	""" Read FITS image and return data """

	try:
		hdu= fits.open(filename)
	except Exception as ex:
		errmsg= 'Cannot read image file: ' + filename
		print "ERROR: " + errmsg
		raise IOError(errmsg)

	data= hdu[0].data
	data_size= np.shape(data)
	nchan= len(data.shape)
	if nchan==4:
		output_data= data[0,0,:,:]
	elif nchan==2:
		output_data= data	
	else:
		errmsg= 'Invalid/unsupported number of channels found in file ' + filename + ' (nchan=' + str(nchan) + ')!'
		print "ERROR: " + errmsg
		raise IOError(errmsg)

	return output_data

# - Get cropped image data 
def crop_img(data,x0,y0,dx,dy):
	""" Extract sub image of size (dx,dy) around pixel (x0,y0) """

	#- Extract crop data
	xmin= int(x0-dx/2)
	xmax= int(x0+dx/2)
	ymin= int(y0-dy/2)
	ymax= int(y0+dy/2)		
	crop_data= data[ymin:ymax,xmin:xmax]
	
	print ('DEBUG: (xmin,xmax)=(%d,%d), (ymin,ymax)=(%d,%d)' % (xmin,xmax,ymin,ymax))

	#- Replace NAN with zeros and inf with large numbers
	np.nan_to_num(crop_data)

	return crop_data

# - Set bkg train data
def make_bkg_train_data(imglist,nsamples,imgsizex,imgsizey,writeimg=False):
	""" Prepare bkg train data """

	# - Read data in list
	imgcounter= 0
	for item in imglist:
		imgcounter+= 1
		filename= item[0]
		print ('INFO: Reading file %s ...' % filename) 

		# - Read main bkg img
		try:
			data= read_img(filename)
		except Exception as ex:
			errmsg= 'Failed to read image data (err=' + str(ex) + ')'
			print "ERROR: " + errmsg
			raise IOError(errmsg)
	
		imgsize= np.shape(data)
		nx= imgsize[1]
		ny= imgsize[0]
		marginx= imgsizex/2
		marginy= imgsizey/2
		print ('INFO: Bkg image no. %d has size (%d,%d)' % (imgcounter,nx,ny) )	

		# - Extract nsamples per img
		index= 0
		while index < nsamples:
			if index%100==0 :
				print ("INFO: Generating sample image no. %s/%s from image %d ..." % (index+1,nsamples,imgcounter))
	
			# - Generate crop img center randomly
			x0= int(np.random.uniform(marginx,nx-marginx-1))
			y0= int(np.random.uniform(marginy,ny-marginy-1))
			print ('DEBUG: (x0,y0)=(%d,%d)' % (x0,y0))

			# - Extract crop img data
			data_crop= crop_img(data,x0,y0,imgsizex,imgsizey)
			imgcropsize= np.shape(data_crop)
			print ('INFO: Sample image no. %s/%s from image %d has size (%d,%d)' % (index+1,nsamples,imgcounter,imgcropsize[1],imgcropsize[0]) )	
 
			# - Check data integrity (skip if all zeros or nan/inf)
			n_nonzero= np.count_nonzero(data_crop)
			n_finite= (np.isfinite(data_crop)).sum()
			if n_nonzero<=0 or n_finite<=0:
				print ('WARN: Skip sample image (all pixels NaN/inf/zero)...')
				continue

			# - Save crop img to file?
			outfilename= 'train_bkg_' + str(index+1) + '-RUN' + str(imgcounter) + '.fits'
			if writeimg:
				write_fits(data_crop,outfilename)

			# ...
			# ...

			index+= 1

##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	#===========================
	#==   PARSE ARGS
	#===========================
	print('INFO: Get script args ...')
	try:
		args= get_args()
	except Exception as ex:
		print("Failed to get and parse options (err=%s)",str(ex))
		return 1

	filelist_bkg= args.filelist_bkg
	filelist_source= args.filelist_source
	nx= args.nx
	ny= args.ny
	nsamples= args.nsamples
	saveimg= args.saveimg
	
	#===========================
	#==   READ BKG/SOURCE DATA
	#===========================
	print ('INFO: Reading bkg data from file %s ...' % filelist_bkg)
	try:
		filenames_bkg= read_bkg_data(filelist_bkg)
	except Exception as ex:
		print("Failed to read bkg data (err=%s)",str(ex))
		return 1
		
	print ('INFO: Reading source data from file %s ...' % filelist_source)
	try:
		read_source_data(filelist_source)
	except Exception as ex:
		print("Failed to read source data (err=%s)",str(ex))
		return 1

	#===========================
	#==   GENERATE TRAIN DATA
	#===========================
	print ('INFO: Generating bkg train data ...')
	make_bkg_train_data(
		imglist=filenames_bkg,
		nsamples=nsamples,
		imgsizex=nx,imgsizey=ny,
		writeimg=saveimg
	)



###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

