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
from astropy.units import Quantity
#from astropy import wcs
from astropy import units as u
from astropy.wcs import WCS
import radio_beam

## OPENCV
import cv2 as cv

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## LOGGER
import logging
import logging.config
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s - %(message)s",datefmt='%Y-%m-%d %H:%M:%S')
logger= logging.getLogger(__name__)
logger.setLevel(logging.INFO)

##################################################

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
	parser.add_argument('-filenames', '--filenames', dest='filenames', required=True, type=str, action='store',help='List of FITS image names (separated by commas) to be convolved to common resolution')
	
	parser.add_argument('--userbeam', dest='userbeam', action='store_true')	
	parser.set_defaults(userbeam=False)
	parser.add_argument('-bmaj', '--bmaj', dest='bmaj', required=False, type=float, default=None, action='store',help='Desired common beam bmaj in arcsec')
	parser.add_argument('-bmin', '--bmin', dest='bmin', required=False, type=float, default=None,action='store',help='Desired common beam bmin in arcsec')
	parser.add_argument('-pa', '--pa', dest='pa', required=False, type=float, default=None, action='store',help='Desired common beam position angle in deg')
	

	args = parser.parse_args()	

	return args

#######################
##   UTILITY CLASS   ##
#######################
class Utils(object):
	""" Class collecting utility methods """

	def __init__(self):
		""" Return a Utils object """

	@classmethod
	def getBaseFileNoExt(cls, filename):
		""" Get basefilename without extension """
		filename_base = os.path.basename(filename)
		filename_base_noext = os.path.splitext(filename_base)[0]
		return filename_base_noext

	@classmethod
	def read_fits(cls, filename):
		""" Read FITS image and return data """

		# - Open file
		try:
			hdu = fits.open(filename, memmap=False)
		except Exception as ex:
			errmsg = 'Cannot read image file: ' + filename
			logger.error(errmsg)
			raise IOError(errmsg)

		# - Read metadata
		header = hdu[0].header
		output_header= header

		# - Read data
		data = hdu[0].data
		data_size = np.shape(data)
		nchan = len(data.shape)
		if nchan == 4:
			output_data = data[0, 0, :, :]
			#output_data = data
		
		elif nchan == 2:
			output_data = data
		else:
			errmsg = 'Invalid/unsupported number of channels found in file ' + filename + ' (nchan=' + str(nchan) + ')!'
			logger.error(errmsg)
			hdu.close()
			raise IOError(errmsg)

			logger.info("Deleting NAXIS3/NAXIS4 from header...")
			if 'NAXIS3' in output_header:
				del output_header['NAXIS3']
			if 'NAXIS4' in output_header:
				del output_header['NAXIS4']
			if 'CTYPE3' in output_header:
				del output_header['CTYPE3']
			if 'CRVAL3' in output_header:
				del output_header['CRVAL3']
			if 'CDELT3' in output_header:
				del output_header['CDELT3']
			if 'CRPIX3' in output_header:
				#del output_header['CRPIX3']
				output_header.pop('CRPIX3', None)
			if 'CROTA3' in output_header:
				del output_header['CROTA3']
			if 'CTYPE4' in output_header:
				del output_header['CTYPE4']
			if 'CRVAL4' in output_header:
				del output_header['CRVAL4']
			if 'CDELT4' in output_header:
				del output_header['CDELT4']
			if 'CRPIX4' in output_header:
				del output_header['CRPIX4']
			if 'CROTA4' in output_header:
				del output_header['CROTA4']

		#print(output_header)

		#print(output_header["NAXIS3"])
		#print(output_header["NAXIS4"])
		print(output_header["CRPIX3"])
		#print(output_header["CRPIX4"])


		# - Close file
		hdu.close()

		return output_data, output_header

	@classmethod
	def hasBeamInfo(cls, header):
		""" Check if header has beam information """
		hasBmaj = ('BMAJ' in header)
		hasBmin = ('BMIN' in header)
		hasBeamInfo = (hasBmaj and hasBmin)
		if not hasBeamInfo:
			return False

		hasBeamVal = (header['BMAJ'] and header['BMIN'])
		if not hasBeamVal:
			return False

		return True

	@classmethod
	def write_fits(cls, data, filename, header=None):
		""" Read data to FITS image """
		if header:
			hdu = fits.PrimaryHDU(data, header)
		else:
			hdu = fits.PrimaryHDU(data)

		hdul = fits.HDUList([hdu])
		hdul.writeto(filename, overwrite=True)

##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	#===========================
	#==   PARSE ARGS
	#===========================
	logger.info("Get script args ...")
	try:
		args= get_args()
	except Exception as ex:
		logger.error("Failed to get and parse options (err=%s)",str(ex))
		return 1

	filenames= [x.strip() for x in args.filenames.split(',')]
	userbeam= args.userbeam
	bmaj_user= args.bmaj
	bmin_user= args.bmin
	pa_user= args.pa

	print("=== ARGS ===")
	print("filenames")
	print(filenames)
	print("userbeam? %d" % userbeam)
	if userbeam:
		print("Beam (bmaj=%f arcsec,bmin=%f arcsec,pa=%f deg)" % (bmaj_user,bmin_user,pa_user))
	print("============")
	

	#===========================
	#==   READ IMAGES
	#===========================	
	filenames_out= []
	data_list= []
	header_list= []
	pixsize_x= []
	pixsize_y= []
	beam_list= []

	for filename in filenames:
		# - Store output filename
		filename_out= Utils.getBaseFileNoExt(filename) + '_conv.fits'
		filenames_out.append(filename_out)

		# - Get image data & header
		#data_alldim, header= Utils.read_fits(filename)
		#nchan = len(data_alldim.shape)
		
		data, header= Utils.read_fits(filename)
		nchan = len(data.shape)
		
		logger.info("nchan=%d" % nchan)

		#if nchan==4:
		#	data= data_alldim[0, 0, :, :]
		#else:	
		#	data= data_alldim
		data[np.isnan(data)]= 0.0 # replace all NAN pixels with 0
		data_list.append(data)
		header_list.append(header)

		# - Get beam and WCS info		
		wcs = WCS(header)
		hasBeamInfo= Utils.hasBeamInfo(header)
		xc= header['CRPIX1']
		yc= header['CRPIX2']
		if nchan==4:
			ra, dec = wcs.all_pix2world(xc,yc,0,0,0,ra_dec_order=True)
		else:	
			ra, dec = wcs.all_pix2world(xc,yc,0,ra_dec_order=True)
		print("ra=%f, dec=%f" % (ra,dec))
		dx= abs(header['CDELT1']) # in deg
		dy= abs(header['CDELT2']) # in deg
		pixsize_x.append(dx)
		pixsize_y.append(dy)

		if hasBeamInfo:
			bmaj= header['BMAJ'] # in deg
			bmin= header['BMIN'] # in deg
			pa= header['BPA'] if 'BPA' in header else 0 # in deg
			beam= radio_beam.Beam(bmaj*u.deg,bmin*u.deg,pa*u.deg)
			beam_list.append(beam)
		else:
			logger.error("No BMAJ/BMIN keyword present in file " + filename + "!")
			return -1


	#===========================
	#==   SET COMMON BEAM
	#===========================	
	if userbeam:
		#common_beam_bmaj= bmaj_user.to(u.arcsec).value
		#common_beam_bmin= bmin_user.to(u.arcsec).value
		#common_beam_pa= pa_user.to(u.deg).value
		common_beam_bmaj= bmaj_user
		common_beam_bmin= bmin_user
		common_beam_pa= pa_user
		common_beam= radio_beam.Beam(bmaj_user*u.arcsec,bmin_user*u.arcsec,pa_user*u.deg)
	else:
		beams= radio_beam.Beams(beams=beam_list)
		common_beam= radio_beam.commonbeam.common_manybeams_mve(beams)
		common_beam_bmaj= common_beam.major.to(u.arcsec).value
		common_beam_bmin= common_beam.minor.to(u.arcsec).value
		common_beam_pa= common_beam.pa.to(u.deg).value
		
	logger.info("Convolving images to common beam size (bmaj,bmin,pa)=(%s,%s,%s) ..." % (str(common_beam_bmaj),str(common_beam_bmin),str(common_beam_pa)))
		

	#===========================
	#==   CONVOLVE IMAGES
	#===========================	
	for index in range(0,len(data_list)):

		# - Find convolving beam for this image
		bmaj, bmin, pa= radio_beam.utils.deconvolve(common_beam,beam_list[index])
		bmaj_deg= bmaj.to(u.deg).value
		bmin_deg= bmin.to(u.deg).value
		pa_deg= pa.to(u.deg).value
		conv_beam= radio_beam.Beam(bmaj_deg*u.deg,bmin_deg*u.deg,pa_deg*u.deg)
			
		bmaj_arcsec= bmaj.to(u.arcsec).value
		bmin_arcsec= bmin.to(u.arcsec).value
		ny= data_list[index].shape[0]
		nx= data_list[index].shape[1]
		
		# - Create convolution kernel
		dx= pixsize_x[index]
		dy= pixsize_y[index]
		pixsize= max(dx,dy)
		conv_kernel= conv_beam.as_kernel(pixsize*u.deg)
		conv_kernel.normalize()
		kernel=conv_kernel.array
		logger.info("Convolution kernel size: %d x %d" % (kernel.shape[0],kernel.shape[1]))
		
		# - Convolve image 
		logger.info("Convolving image %d (size=%d,%d) by beam (bmaj,bmin,pa)=(%s,%s,%s) ..." % (index+1,nx,ny,str(bmaj_arcsec),str(bmin_arcsec),str(pa_deg)))
		data_conv= cv.filter2D(np.float64(data_list[index]),-1,kernel,borderType=cv.BORDER_CONSTANT)
		
		# - Write output FITS	
		logger.info("Saving convolved image to file %s ..." % filenames_out[index])
		Utils.write_fits(data_conv,filenames_out[index],header_list[index])

		

	return 0


###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

