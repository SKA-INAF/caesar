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
##from ctypes import *

## ASTRO
#from scipy import ndimage
#from astropy.io import fits
#from astropy.units import Quantity
#from astropy.modeling.parameters import Parameter
#from astropy.modeling.core import Fittable2DModel
#from astropy.modeling.models import Box2D, Gaussian2D, Ring2D, Ellipse2D, TrapezoidDisk2D, Disk2D, AiryDisk2D, Sersic2D
#from astropy import wcs

## ROOT
#import ROOT
#from ROOT import gSystem, TFile, TTree, gROOT, AddressOf

## CAESAR
#gSystem.Load('libCaesar')
#from ROOT import Caesar

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## Graphics modules
#import matplotlib.pyplot as plt
#import pylab
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

	## Needed to trick CASA when running in batch	
	parser.add_argument('-c', dest='scriptname', required=False, type=str, default='',action='store',help='Script name')

	# - RUN
	parser.add_argument('--no-convolve', dest='do_convolve', action='store_false')	
	parser.set_defaults(do_convolve=True)

	# - INPUTS
	parser.add_argument('-mosaic', '--mosaic', dest='mosaic', required=False, type=str, default='', action='store',help='Mosaic Image in units Jy/beam (.fits)')
	parser.add_argument('-residual', '--residual', dest='residual', required=False, type=str, default='', action='store',help='Residual Image in units Jy/pixel (.fits)')
	parser.add_argument('-model', '--model', dest='model', required=False, type=str, default='', action='store',help='Sky Model Image (if given override model generation) in units Jy/pixel (.fits)')
	parser.add_argument('-img', '--img', dest='img', required=False, type=str, default='', action='store',help='Restored Image (.fits)')
		

	# - POINT SOURCE GENERATION OPTIONS
	parser.add_argument('-marginx', '--marginx', dest='marginx', required=False, type=int, default=0,action='store',help='Image x margin in pixels')
	parser.add_argument('-marginy', '--marginy', dest='marginy', required=False, type=int, default=0,action='store',help='Image y margin in pixels')
	parser.add_argument('-nsources', '--nsources', dest='nsources', required=False, type=int, default=0, action='store',help='Number of point-source to be generated (if >0 overrides the density generation) (default=0)')
	parser.add_argument('-source_density', '--source_density', dest='source_density', required=False, type=float, default=1000, action='store',help='Compact source density (default=1000)')
	parser.add_argument('-Smin', '--Smin', dest='Smin', required=False, type=float, default=1.e-6, action='store',help='Minimum source flux in Jy (default=1.e-6)')
	parser.add_argument('-Smax', '--Smax', dest='Smax', required=False, type=float, default=1, action='store',help='Maximum source flux in Jy (default=1)')
	parser.add_argument('-bmaj', '--bmaj', dest='bmaj', required=False, type=float, default=6.5, action='store',help='Beam bmaj in arcsec (default=6.5)')
	parser.add_argument('-bmin', '--bmin', dest='bmin', required=False, type=float, default=6.5, action='store',help='Beam bmin in arcsec (default=6.5)')
	parser.add_argument('-pa', '--pa', dest='pa', required=False, type=float, default=0, action='store',help='Beam position angle in deg (default=0)')

	parser.add_argument('--extsources', dest='extsources', action='store_true')	
	parser.set_defaults(extsources=False)
	parser.add_argument('-nsources_ext', '--nsources_ext', dest='nsources_ext', required=False, type=int, default=0, action='store',help='Number of extended source to be generated (if >0 overrides the density generation) (default=0)')
	parser.add_argument('-source_density_ext', '--source_density_ext', dest='source_density_ext', required=False, type=float, default=50, action='store',help='Extended source density (default=50)')
	parser.add_argument('-Smin_ext', '--Smin_ext', dest='Smin_ext', required=False, type=float, default=1.e-6, action='store',help='Minimum extended source flux in Jy (default=1.e-6)')
	parser.add_argument('-Smax_ext', '--Smax_ext', dest='Smax_ext', required=False, type=float, default=1, action='store',help='Maximum extended source flux in Jy (default=1)')
	parser.add_argument('-bmaj_min', '--bmaj_min', dest='bmaj_min', required=False, type=float, default=4, action='store',help='Gaussian components min bmaj in arcsec (default=4)')
	parser.add_argument('-bmaj_max', '--bmaj_max', dest='bmaj_max', required=False, type=float, default=10, action='store',help='Gaussian components max bmaj in arcsec (default=10)')
	parser.add_argument('-bmin_min', '--bmin_min', dest='bmin_min', required=False, type=float, default=4, action='store',help='Gaussian components  min bmin in arcsec (default=4)')
	parser.add_argument('-bmin_max', '--bmin_max', dest='bmin_max', required=False, type=float, default=10, action='store',help='Gaussian components  max bmin in arcsec (default=10)')
	parser.add_argument('-pa_min', '--pa_min', dest='pa_min', required=False, type=float, default=-90, action='store',help='Gaussian components  min position angle in deg (default=0)')
	parser.add_argument('-pa_max', '--pa_max', dest='pa_max', required=False, type=float, default=90, action='store',help='Gaussian components  max position angle in deg (default=180)')
	

	# - OUTPUT FILE OPTIONS
	parser.add_argument('-outfile_img', '--outfile_img', dest='outfile_img', required=False, type=str, default='simmap.fits',action='store',help='Output filename')
	parser.add_argument('-outfile_model', '--outfile_model', dest='outfile_model', required=False, type=str, default='skymodel.fits', action='store',help='Model filename')
	parser.add_argument('-outfile_model_ext', '--outfile_model_ext', dest='outfile_model_ext', required=False, type=str, default='extsourcemodel.fits', action='store',help='Extended source model filename')
	parser.add_argument('-outfile_modelconv', '--outfile_modelconv', dest='outfile_modelconv', required=False, type=str, default='skymodel_conv.fits', action='store',help='Model convolved filename')
	parser.add_argument('-outfile_sources', '--outfile_sources', dest='outfile_sources', required=False, type=str, default='sources.dat',action='store',help='Ascii file with list of generated sources')
	parser.add_argument('-outfile_sourcepars', '--outfile_sourcepars', dest='outfile_sourcepars', required=False, type=str, default='source_pars.dat',action='store',help='Ascii file with list of generated source pars')
	parser.add_argument('-outfile_ds9region', '--outfile_ds9region', dest='outfile_ds9region', required=False, type=str, default='dsregion.reg',action='store',help='DS9 source region filename')
	
	args = parser.parse_args()	

	return args




###########################
##     SIMULATOR CLASS
###########################
SIGMA_TO_FWHM= np.sqrt(8*np.log(2))

class SkyMapSimulator(object):

	""" Sky map simulator class

			Attributes:
				mosaic_file: mosaic image (.fits)
				ny: image height in pixels
				pixsize: pixel size in arcsec (default=1)
	"""

	def __init__(self, mosaic_file):
		""" Return a SkyMapGenerator object """

		## Input maps
		self.mosaic_file= mosaic_file
		self.mosaic_im= None
		self.mosaic_cs= None

		self.model_data= None
		self.model_im= None

		## Convolved map
		self.model_conv_data= None

		## Ext source map
		self.model_data_ext= None
		self.model_ext_im= None

		## Image parameters
		self.nx= 0
		self.ny= 0
		self.marginx= 0 # in pixels (no margin)
		self.marginy= 0 # in pixels (no margin)
		self.pixsize= 0 # in arcsec
		self.gridx= []
		self.gridy= []
		self.beam_bmaj= 6.5 # in arcsec
		self.beam_bmin= 6.5 # in arcsec
		self.beam_bpa= 0 # in deg
		self.beam_area= self.compute_beam_area(self.beam_bmaj,self.beam_bmin) # in pixels
		
		## Compact source parameters
		self.generate_skymodel= True
		self.ps_list= []
		self.nsources= 0 # default is density generator
		self.source_density= 2000. # in sources/deg^2
		self.Smin= 1.e-6 # in Jy 
		self.Smax= 1 # in Jy

		## Extended source parameters
		self.add_ext_sources= False
		self.truncate_models= False
		self.exts_list= []
		self.trunc_thr= 0.01 # 1% flux truncation at maximum
		self.nsources_ext= 0 # default is density generator
		self.source_density_ext= 50. # in sources/deg^2
		self.beam_bpa_min= -90 # deg
		self.beam_bpa_max= 90 # deg
		self.beam_bmaj_min= 4	# arcsec
		self.beam_bmaj_max= 10 # arcsec
		self.beam_bmin_min= 4	# arcsec 
		self.beam_bmin_max= 10 # arcsec
		self.Smin_ext= 1.e-6 # in Jy 
		self.Smax_ext= 1 # in Jy
		
		## Map output file
		self.model_outfile= 'skymodel.fits'
		self.model_ext_outfile= 'extsources.fits'
		self.img_outfile= 'simmap.fits'
		self.img_im= None	

		## DS9 & ascii output file
		self.source_outfile= 'sources.dat'
		self.ds9_outfile= 'ds9region.reg'
		self.source_par_outfile= 'source_pars.dat'
	
	def set_margins(self,marginx,marginy):
		""" Set margin in X & Y """
		#if (marginx<0 or marginy<0 or marginx>=self.nx/2 or marginy>=self.ny/2) :
		#	raise ValueError('Invalid margin specified (<0 or larger than image half size!')
		self.marginx= marginx
		self.marginy= marginy

	def set_ds9region_outfile(self,filename):
		""" Set the output DS9 region filename """
		self.ds9_outfile= filename

	def set_source_outfile(self,filename):
		""" Set the output source list filename """
		self.source_outfile= filename

	def set_source_par_outfile(self,filename):
		""" Set the output source par list filename """
		self.source_par_outfile= filename

	def set_model_outfile(self,filename):
		""" Set the output model filename """
		self.model_outfile= filename

	def set_ext_model_outfile(self,filename):
		""" Set the output extended source model filename """
		self.model_ext_outfile= filename

	def set_img_outfile(self,filename):
		""" Set the output image filename """
		self.img_outfile= filename

	
	def set_source_flux_range(self,Smin,Smax):
		""" Set source flux range """
		self.Smin= Smin
		self.Smax= Smax

	def set_ext_source_flux_range(self,Smin,Smax):
		""" Set extended source flux range """
		self.Smin_ext= Smin
		self.Smax_ext= Smax

	def set_nsources(self,n):
		""" Set number of sources to be generated """
		if n<0:	
			raise ValueError('Invalid number of sources specified (shall be >=0)')
		self.nsources= n

	def set_nsources_ext(self,n):
		""" Set number of extended sources to be generated """
		if n<0:	
			raise ValueError('Invalid number of sources specified (shall be >=0)')
		self.nsources_ext= n

	def set_source_density(self,density):
		""" Set compact source density in deg^-2 """
		self.source_density= density

	def set_ext_source_density(self,density):
		""" Set extended source density in deg^-2 """
		self.source_density_ext= density

	def compute_beam_area(self,Bmaj,Bmin):
		""" Compute beam area """
		A= np.pi*Bmaj*Bmin/(4*np.log(2)) #2d gaussian area with FWHM=fx,fy (in arcsec^2)
		pixelArea= np.fabs(self.pixsize*self.pixsize) # in arcsec^2
		beam_area= A/pixelArea # in pixels
		return beam_area

	def set_beam(self,bmaj,bmin,pa):
		""" Set beam pars """
		self.beam_bmaj= bmaj
		self.beam_bmin= bmin
		self.beam_bpa= pa
		self.beam_area= self.compute_beam_area(bmaj,bmin)

	def set_beam_bmaj_range(self,bmaj_min,bmaj_max):
		""" Set beam bmaj range """
		self.beam_bmaj_min= bmaj_min
		self.beam_bmaj_max= bmaj_max	
	
	def set_beam_bmin_range(self,bmin_min,bmin_max):
		""" Set beam bmin range """
		self.beam_bmin_min= bmin_min
		self.beam_bmin_max= bmin_max	

	def set_beam_pa_range(self,pa_min,pa_max):
		""" Set beam pa range """
		self.beam_bpa_min= pa_min
		self.beam_bpa_max= pa_max	

	def add_ext_sources(self,choice):
		""" Set beam randomization """
		self.add_ext_sources= choice

	def write_source_list(self):
		""" Write source list to file """

		## Open file	
		fout = open(self.source_outfile, 'wb')
	
		## Write point source data
		for i in range(len(self.ps_list)):
			name= self.ps_list[i][0]
			x= self.ps_list[i][1]
			y= self.ps_list[i][2]
			S= self.ps_list[i][3]
			x_wcs= self.model_im.toworld([x,y])['numeric'][0]
			y_wcs= self.model_im.toworld([x,y])['numeric'][1]
			x_wcs= np.rad2deg(x_wcs)
			y_wcs= np.rad2deg(y_wcs)
			
			data= (("%s %s %s %s %s %s") % (name,x,y,x_wcs,y_wcs,S) )

			fout.write(data)
			fout.write('\n')	

		## Write ext sources
		for i in range(len(self.exts_list)):
			name= self.exts_list[i][0]
			x= self.exts_list[i][1]
			y= self.exts_list[i][2]
			S= self.exts_list[i][3]
			x_wcs= self.model_im.toworld([x,y])['numeric'][0]
			y_wcs= self.model_im.toworld([x,y])['numeric'][1]
			x_wcs= np.rad2deg(x_wcs)
			y_wcs= np.rad2deg(y_wcs)

			data= (("%s %s %s %s %s %s") % (name,x,y,x_wcs,y_wcs,S) )

			fout.write(data)
			fout.write('\n')	

		fout.close();


	def write_source_par_list(self):
		""" Write source par list to file """
		## Open file	
		fout = open(self.self.source_par_outfile, 'wb')
	
		## Write header
		header= ("# name x(pix) y(pix) S(Jy/beam) sigmax(pix) sigmay(pix) theta(deg)")
		fout.write(header)
		fout.write('\n')	

		## Write point-source pars
		sigmax= self.beam_bmaj/(self.pixsize*SIGMA_TO_FWHM) # in pixels
		sigmay= self.beam_bmin/(self.pixsize*SIGMA_TO_FWHM) # in pixels
		theta= 90 + self.beam_bpa # in deg
		scaleFactor= self.beam_area

		for i in range(len(self.ps_list)):
			name= self.ps_list[i][0]
			x= self.ps_list[i][1]
			y= self.ps_list[i][2]
			S= self.ps_list[i][3]*scaleFactor # Convert from Jy/pixel to Jy/beam
			
			data= (("%s %s %s %s %s %s %s") % (name,x,y,S,sigmax,sigmay,theta) )

			fout.write(data)
			fout.write('\n')
			fout.flush()			

		## Write ext sources
		for i in range(len(self.exts_list)):
			name= self.exts_list[i][0]
			x= self.exts_list[i][1]
			y= self.exts_list[i][2]
			S= self.exts_list[i][3] # Fluxes are already in Jy/beam
			sigmax= self.exts_list[i][4] 
			sigmay= self.exts_list[i][5]
			theta= self.exts_list[i][6]
			
			data= (("%s %s %s %s %s %s %s") % (name,x,y,S,sigmax,sigmay,theta) )

			fout.write(data)
			fout.write('\n')	
			fout.flush()

		fout.close();



	def write_ds9_regions(self):
		""" Write DS9 region """
	
		## Open file	
		fout = open(self.ds9_outfile, 'wb')

		ds9_outfile_dir= os.path.dirname(self.ds9_outfile)
		ds9_outfile_base= os.path.basename(self.ds9_outfile)
		wcs_ds9_outfile= ds9_outfile_dir + '/wcs_' + ds9_outfile_base

		fout_wcs= open(wcs_ds9_outfile, 'wb')
		
		## Write file header
		fout.write('global color=red font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n')
		fout.write('image\n')

		cs= self.model_im.coordsys() 
		wcs_type= cs.referencecode()[0]
		wcs_type= wcs_type.lower()
		wcs_type_str= (("%s\n") % wcs_type)

		fout_wcs.write('global color=red font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n')		
		fout_wcs.write(wcs_type_str)

		## Write point sources
		for i in range(len(self.ps_list)):
			name= self.ps_list[i][0]
			x= self.ps_list[i][1]
			y= self.ps_list[i][2]
			x_wcs= self.model_im.toworld([x,y])['numeric'][0]
			y_wcs= self.model_im.toworld([x,y])['numeric'][1]
			x_wcs= np.rad2deg(x_wcs)
			y_wcs= np.rad2deg(y_wcs)
			
			x_ds9= x+1 # NB: DS9 starts index from 1
			y_ds9= y+1
			S= self.ps_list[i][3]
			
			# point=[circle|box|diamond|cross|x|arrow|boxcircle] [size]
			#region= (("circle point %s %s # %s") % (x,y,name) )
			bmaj_pix= self.beam_bmaj/self.pixsize # in pixels
			bmin_pix= self.beam_bmin/self.pixsize # in pixels
			sigmax_pix= self.beam_bmaj/(self.pixsize*SIGMA_TO_FWHM)
			sigmay_pix= self.beam_bmin/(self.pixsize*SIGMA_TO_FWHM)
			ell_semimajor_axis_pix= bmaj_pix/2
			ell_semiminor_axis_pix= bmin_pix/2
			theta= self.beam_bpa - 90 ## to convert to DS9 format 
			region= (("ellipse %s %s %s %s %s # text={%s}") % (x_ds9,y_ds9,ell_semimajor_axis_pix,ell_semiminor_axis_pix,theta,name) )
			
			fout.write(region)
			fout.write('\n')	
			fout.flush()	

			## Write WCS DS9
			bmaj_wcs= self.beam_bmaj
			bmin_wcs= self.beam_bmin
			ell_semimajor_axis_wcs= bmaj_wcs/2
			ell_semiminor_axis_wcs= bmin_wcs/2
			theta_wcs= self.beam_bpa - 90 ## to convert to DS9 format 
			region_wcs= (("ellipse(%s %s %s\" %s\" %s) # text={%s}") % (x_wcs,y_wcs,ell_semimajor_axis_wcs,ell_semiminor_axis_wcs,theta_wcs,name) )
			fout_wcs.write(region_wcs)
			fout_wcs.write('\n')	
			fout_wcs.flush()
		

		## Write ext sources
		for i in range(len(self.exts_list)):
			name= self.exts_list[i][0]
			x= self.exts_list[i][1]
			y= self.exts_list[i][2]
			x_wcs= self.model_im.toworld([x,y])['numeric'][0]
			y_wcs= self.model_im.toworld([x,y])['numeric'][1]
			x_wcs= np.rad2deg(x_wcs)
			y_wcs= np.rad2deg(y_wcs)
			
			x_ds9= x+1 # NB: DS9 starts index from 1
			y_ds9= y+1
			S= self.exts_list[i][3]
			sigmax= self.exts_list[i][4] # in pix
			sigmay= self.exts_list[i][5] # in pix
			theta= self.exts_list[i][6]
			
			bmaj_pix= sigmax*SIGMA_TO_FWHM # in pixels
			bmin_pix= sigmay*SIGMA_TO_FWHM # in pixels
			ell_semimajor_axis_pix= bmaj_pix/2
			ell_semiminor_axis_pix= bmin_pix/2
			region= (("ellipse %s %s %s %s %s # text={%s}") % (x_ds9,y_ds9,ell_semimajor_axis_pix,ell_semiminor_axis_pix,theta,name) )
			
			fout.write(region)
			fout.write('\n')	
			fout.flush()

			## Write WCS DS9
			sigmax_wcs= sigmax*self.pixsize # in arcsec
			sigmay_wcs= sigmay*self.pixsize # in arcsec
			bmaj_wcs= sigmax_wcs*SIGMA_TO_FWHM 
			bmin_wcs= sigmay_wcs*SIGMA_TO_FWHM 
			ell_semimajor_axis_wcs= bmaj_wcs/2
			ell_semiminor_axis_wcs= bmin_wcs/2
			theta_wcs= theta
			region_wcs= (("ellipse(%s %s %s\" %s\" %s) # text={%s}") % (x_wcs,y_wcs,ell_semimajor_axis_wcs,ell_semiminor_axis_wcs,theta_wcs,name) )
			fout_wcs.write(region_wcs)
			fout_wcs.write('\n')		
			fout_wcs.flush()

		fout.close();
		fout_wcs.close();

	
	def write_data_to_fits(self,im,outputfile):
		""" Write casa image to FITS """
		im.tofits(outputfile,overwrite=True)
		

	def init(self):
		""" Initialize data """
		
		## Check if mosaic map is empty
		if not self.mosaic_file:
			raise ValueError('Missing mosaic map filename!')

		## Read mosaic map from file
		file_ext= os.path.splitext(self.mosaic_file)[1]
		if file_ext=='.fits':
			print 'INFO: Importing mosaic FITS file %s as CASA image ...' % self.mosaic_file
			self.mosaic_im= ia.newimagefromfits(infile=self.mosaic_file)
		else:
			print 'INFO: Importing mosaic CASA image %s ' % self.mosaic_file
			self.mosaic_im= ia.newimagefromfile(infile=self.mosaic_file)

		## Check image read
		if not self.mosaic_im:
			errmsg= 'ERROR: Failed to import given mosaic input file!'
			print(errmsg)
			raise ValueError(errmsg)
			
		## Read and store mosaic image info 
		self.mosaic_cs= self.mosaic_im.coordsys() 	
		self.mosaic_data= self.mosaic_im.getregion()
		data_shape= self.mosaic_im.shape()
		self.nx= data_shape[0]
		self.ny= data_shape[1]
		dx= np.abs(np.rad2deg(self.mosaic_cs.increment(type='direction')['numeric'][0]))*3600 # in arcsec
		dy= np.abs(np.rad2deg(self.mosaic_cs.increment(type='direction')['numeric'][1]))*3600 # in arcsec
		self.pixsize= min(dx,dy)
		beam= self.mosaic_im.restoringbeam()
		self.beam_bmaj= beam['major']['value'] # in arcsec
		self.beam_bmin= beam['minor']['value'] # in arcsec
		self.beam_bpa= beam['positionangle']['value'] # in deg	
		print 'INFO: Mosaic image shape=', data_shape
		print ('INFO: Restored image beam= (%s,%s,%s)') % (self.beam_bmaj,self.beam_bmin,self.beam_bpa)
					
		## Check margins
		if (self.marginx<0 or self.marginy<0 or self.marginx>=self.nx/2 or self.marginy>=self.ny/2) :
			raise ValueError('Invalid margin specified (<0 or larger than image half size!')
					
		## Initialize grid
		self.gridy, self.gridx = np.mgrid[0:self.ny, 0:self.nx]


	def generate_compact_sources(self):
		""" Generate list of compact sources in the map.
					- Uniform spatial distribution
					- Uniform flux distribution

				Arguments:	
					density: source density in #sources/deg^2 (e.g. 2000)
		"""
		
		# Get mosaic pixels with non-nan values
		mask= ~np.isnan(self.mosaic_data)
		npix_good= np.count_nonzero(mask)

		# Compute number of sources to be generated given map area in pixels
		area= ((self.nx-2*self.marginx)*(self.ny-2*self.marginy))*self.pixsize/(3600.*3600.) # in deg^2
		 
		if self.nsources>0:
			nsources= self.nsources
		else: # density generator
			nsources= int(round(self.source_density*area))
		S_min= self.Smin
		S_max= self.Smax
		lgS_min= np.log(S_min)
		lgS_max= np.log(S_max)
		randomize_flux= False
		if self.Smin<self.Smax:
			randomize_flux= True

		print 'INFO: Generating #',nsources,' compact sources in map...'

		
		## Start generation loop
		self.model_data= ia.makearray(0,[self.nx,self.ny,1,1])
		##self.model_data= Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.nx, y_width=2*self.ny)(self.gridx, self.gridy)
		self.ps_list= []
		index= 0

		#for index in range(0,nsources):	
		while index < nsources:
			if index%100==0 :
				print ("INFO: Generating compact source no. %s/%s" % (index+1,nsources))
			
			## Generate random coordinates
			x0= np.random.uniform(self.marginx,self.nx-self.marginx-1)
			y0= np.random.uniform(self.marginy,self.ny-self.marginy-1)

			## Compute amplitude given significance level and bkg
			## Generate flux uniform in log
			if randomize_flux:
				lgS= np.random.uniform(lgS_min,lgS_max)
				S= np.exp(lgS)
			else:
				S= S_min
				
			## Fill sky model map
			ix= int(np.round(x0))
			iy= int(np.round(y0))
			##if self.model_data[iy,ix]>0:
			if self.model_data[ix][iy][0][0]>0 or not mask[ix][iy][0][0]:
				print ("INFO: Skip source generated @ pixel (%s,%s) as already filled or corresponding to nan value in mosaic ..." % (ix,iy) ) 
				continue

			index+= 1
			source_name= 'S' + str(index)
			#self.model_data[iy,ix]+= S
			self.model_data[ix][iy][0][0]+= S
			self.ps_list.append([source_name,x0,y0,S])
			
			print ('INFO: Source %s: Pos(%s,%s), ix=%s, iy=%s, S=%s' % (source_name,str(x0),str(y0),str(ix),str(iy),str(S)))

		## Create skymodel image
		self.model_im= ia.newimagefromarray(pixels=self.model_data, csys=self.mosaic_cs.torecord())  
		self.model_im.setbrightnessunit('Jy/pixel')
		print ('INFO: model units=%s') % self.model_im.brightnessunit()	
		print 'INFO: model shape=', self.model_im.shape()

		print ('INFO: End point source generation')

		
	#####################################
	###         CONVOLVE SKYMODEL     ##
	#####################################	
	def convolve_skymodel(self):
		""" Convolve sky model with restored beam """
				
		## Convolving skymodel with restoring beam
		## NB: If input map is in Jy/pixel units, the convolve task automatically scale results to Jy/beam (no further scaling needed)
		bmaj= str(self.beam_bmaj) + 'arcsec'
		bmin= str(self.beam_bmin) + 'arcsec'
		bpa= str(self.beam_bpa) + 'deg'
		print ('INFO: Convolving skymodel with restored beam (%s,%s,%s)') % (bmaj,bmin,bpa)
		self.model_im_conv= self.model_im.convolve2d(type='gauss',major=bmaj,minor=bmin,pa=bpa)

		## Get pixel data from convolved image
		self.model_conv_data= self.model_im_conv.getregion()

	#####################################
	###         GAUS 2D MODEL          ##
	#####################################	
	def generate_blob(self,ampl,x0,y0,sigmax,sigmay,theta,trunc_thr=0.01):
		""" Generate a blob 
				Arguments: 
					ampl: peak flux in Jy
					x0, y0: gaussian means in pixels
					sigmax, sigmay: gaussian sigmas in pixels
					theta: rotation in degrees
					trunc_thr: truncation significance threshold
		"""
		data= Gaussian2D(ampl,x0,y0,sigmax,sigmay,theta=math.radians(theta))(self.gridx, self.gridy)
		
		## Truncate data such that sum(data)_trunc/sum(data)<f
		f= trunc_thr 
		if self.truncate_models:
			totFlux= (float)(np.sum(data,axis=None))
			print('Gaus2D total flux=%s' % str(totFlux))

			data_vect_sorted= np.ravel(data)
			data_csum= np.cumsum(data_vect_sorted)/totFlux
			fluxThr= data_vect_sorted[np.argmin(data_csum<f)]
			print('Gaus2D fluxThr=%s' % str(fluxThr))
			data[data<fluxThr] = 0		

		return data

	#################################################
	###   ADD EXTENDED SOURCES (GAUS COMPONENTS)   ##
	#################################################
	def generate_ext_sources(self):
		""" Add extended gaus components to convolved map """

		# Get mosaic pixels with non-nan values
		mask= ~np.isnan(self.mosaic_data)
		npix_good= np.count_nonzero(mask)

		# Compute number of sources to be generated given map area in pixels
		area= ((self.nx-2*self.marginx)*(self.ny-2*self.marginy))*self.pixsize/(3600.*3600.) # in deg^2
		 
		if self.nsources_ext>0:
			nsources= self.nsources_ext
		else: # density generator
			nsources= int(round(self.source_density_ext*area))

		## Set flux generation range
		fluxScaleFactor= self.beam_area
		S_min= self.Smin_ext * fluxScaleFactor # set flux in Jy/beam units
		S_max= self.Smax_ext * fluxScaleFactor # set flux in Jy/beam units
		lgS_min= np.log(S_min)
		lgS_max= np.log(S_max)
		randomize_flux= False
		if self.Smin_ext<self.Smax_ext:
			randomize_flux= True

		## Set gaus pars generation
		randomize_gaus= False
		Bmaj_min= self.beam_bmaj_min
		Bmaj_max= self.beam_bmaj_max
		Bmin_min= self.beam_bmin_min
		Bmin_max= self.beam_bmin_max
		Pa_min= self.beam_bpa_min
		Pa_max= self.beam_bpa_max	
		if self.beam_bmaj_min<self.beam_bmaj_max:
			randomize_gaus= True
		if self.beam_bmin_min<self.beam_bmin_max:
			randomize_gaus= True
		if self.beam_bpa_min<self.beam_bpa_max:
			randomize_gaus= True

		## Generate sources
		print 'INFO: Generating #',nsources,' extended sources in map...'
		index= 0

		while index < nsources:
			if index%100==0 :
				print ("INFO: Generating extended source no. %s/%s" % (index+1,nsources))
			
			## Generate random coordinates
			x0= np.random.uniform(self.marginx,self.nx-self.marginx-1)
			y0= np.random.uniform(self.marginy,self.ny-self.marginy-1)

			## Compute amplitude given significance level and bkg
			## Generate flux uniform in log
			if randomize_flux:
				lgS= np.random.uniform(lgS_min,lgS_max)
				S= np.exp(lgS)
			else:
				S= S_min

			## Generate gaus pars
			bmin= random.uniform(Bmin_min,Bmin_max)
			bmaj= random.uniform(bmin,Bmaj_max)
			pa= random.uniform(Pa_min,Pa_max)
			sigmax= bmaj/(self.pixsize * SIGMA_TO_FWHM)
			sigmay= bmaj/(self.pixsize * SIGMA_TO_FWHM)
			theta = 90 + pa							

			## Fill sky model map
			ix= int(np.round(x0))
			iy= int(np.round(y0))
			if self.model_data[ix][iy][0][0]>0 or not mask[ix][iy][0][0]:
				print ("INFO: Skip ext source generated @ pixel (%s,%s) as already filled or corresponding to nan value in mosaic ..." % (ix,iy) ) 
				continue

			## Generate gaus 2D data
			source_data= self.generate_blob(ampl=S,x0=x0,y0=y0,sigmax=sigmax,sigmay=sigmay,theta=theta,trunc_thr=self.trunc_thr)
			if source_data is None:
				print('Failed to generate Gaus2D (hint: too large trunc threshold), skip and regenerate...')
				continue

			## Add generated source to image and list
			index+= 1
			source_name= 'Sext' + str(index)
			self.model_data_ext+= source_data
			self.exts_list.append([source_name,x0,y0,S,sigmax,sigmay,theta])
			print ('INFO: Ext source %s: Pos(%s,%s), ix=%s, iy=%s, S=%s' % (source_name,str(x0),str(y0),str(ix),str(iy),str(S)))

		## Create image with ext sources
		self.model_ext_im= ia.newimagefromarray(pixels=self.model_data_ext, csys=self.mosaic_cs.torecord())  
		self.model_ext_im.setbrightnessunit('Jy/beam')
		print ('INFO: model units=%s') % self.model_ext_im.brightnessunit()	
		print 'INFO: model shape=', self.model_ext_im.shape()

		print ('INFO: End ext source generation')


	#####################################
	###         RUN           ##
	#####################################
	def run(self):
		""" Generate sky map and convolve with beam """
		
		## == INITIALIZE DATA ==
		print ('INFO: Initializing simulator data...')
		self.init()

		## == GENERATE SKYMODEL WITH POINT SOURCES ==
		print ('INFO: Generating compact sources...')
		self.generate_compact_sources()

		## == CONVOLVE SKY MODEL WITH BEAM ==
		print ('INFO: Convolving skymodel with beam ...')
		self.convolve_skymodel()

		## == ADD EXTENDED SOURCES ==
		if self.add_ext_sources:
			print ('INFO: Generating extended sources...')
			self.generate_ext_sources()	

		## == CREATE FINAL MAP = MOSAIC + CONVOLVED SKY MODEL + EXT SOURCE (IF ENABLED) ===
		data= self.mosaic_data + self.model_conv_data
		if self.add_ext_sources:
			data+= self.model_data_ext			

		cs = self.mosaic_im.coordsys()
		bmaj= str(self.beam_bmaj) + 'arcsec'
		bmin= str(self.beam_bmin) + 'arcsec'
		bpa= str(self.beam_bpa) + 'deg'
		self.img_im= ia.newimagefromarray(pixels=data, csys=cs.torecord())
		self.img_im.setbrightnessunit('Jy/beam')		
		self.img_im.setrestoringbeam(major=bmaj,minor=bmin,pa=bpa) 
			
		## == WRITE MAPS TO FITS FILES ==
		print ('INFO: Writing model to FITS...')
		self.write_data_to_fits(self.model_im,self.model_outfile)

		print ('INFO: Writing ext source map to FITS...')
		self.write_data_to_fits(self.model_ext_im,self.model_ext_outfile)
		
		print ('INFO: Writing image to FITS...')
		self.write_data_to_fits(self.img_im,self.img_outfile)
		

		## == WRITE DS9 REGION FILE ==
		print ('INFO: Writing DS9 regions ...')
		self.write_ds9_regions()

		print ('INFO: Writing source list ...')
		self.write_source_list()

		## Close open CASA images	
		print ('INFO: Closing CASA images ...')
		self.mosaic_im.done()
		self.model_im.done()
		self.img_im.done()
		if self.model_ext_im:
			self.model_ext_im.done()



###########################
##     IMAGER CLASS
###########################

class SkyMapImager(object):

	""" Sky map imager class

			Attributes:
				residual: beam residual image (.fits)
				model: skymodel image (.fits)
				img_file: beam restored image (.fits)
	"""

	def __init__(self, residual_file, model_file, img_file):
		""" Return a SkyMapImager object """

		## Input maps
		self.residual_file= residual_file
		self.residual_im= None
		self.residual_data= None
		
		self.model_file= model_file
		self.model_im= None
		self.model_data= None

		self.img_file= img_file
		self.img_im= None
		
		self.model_im_conv= None
		self.model_conv_data= None
		
		self.img_conv_im= None

		
		## Initialize restoring beam
		self.beam_bmaj= 6.5 # in arcsec
		self.beam_bmin= 6.5 # in arcsec
		self.beam_bpa= 0 # in deg
		
		## Map output file
		self.img_outfile= 'simmap.fits'
		self.img_outfile_casa= 'simmap'
		self.model_conv_outfile= 'skymodel_conv.fits'
		
		
	def set_map_outfile(self,filename):
		""" Set the output map filename """
		filename_base= os.path.splitext(filename)[0]
		filename_ext= os.path.splitext(filename)[1]
		if filename_ext!='.fits':
			raise ValueError('Invalid outfile specified (missing .fits extension)!')
		self.img_outfile= filename
		self.img_outfile_casa= str(filename_base)

	
	def set_modelconv_outfile(self,filename):
		""" Set the output model colvolved filename """
		self.model_conv_outfile= filename

	
	def write_data_to_fits(self,im,outputfile):
		""" Write casa image to FITS """
		im.tofits(outputfile,overwrite=True)
		

	def init(self):
		""" Initialize imager """

		## Check mandatory filename inputs		
		if not self.residual_file:
			raise ValueError('Missing residual map filename!')
		
		if not self.img_file:			
			raise ValueError('Missing restored image input filename!')
			
		if not self.model_file:
			raise ValueError('Missing skymodel image input filename!')



		
		#####################################
		## Read residual map from file
		#####################################
		file_ext= os.path.splitext(self.residual_file)[1]
		if file_ext=='.fits':
			print 'INFO: Importing residual FITS file %s as CASA image ...' % self.residual_file
			self.residual_im= ia.newimagefromfits(infile=self.residual_file)
		else:
			print 'INFO: Importing residual CASA image %s ' % self.residual_file
			self.residual_im= ia.newimagefromfile(infile=self.residual_file)

		self.residual_data= self.residual_im.getregion()
		print 'INFO: Residual image shape=', self.residual_im.shape()	

		## Get bounding box (used to extract submap from model)
		bb= self.residual_im.boundingbox()
		bl_coords= bb['blcf'].split(',')
		tr_coords= bb['trcf'].split(',')
		x0_wcs= bl_coords[0]
		x1_wcs= tr_coords[0]
		y0_wcs= bl_coords[1]
		y1_wcs= tr_coords[1]
		cs= self.residual_im.coordsys()
		xref_wcs= np.rad2deg(cs.referencevalue(type='direction')['numeric'][0]) # in deg
		yref_wcs= np.rad2deg(cs.referencevalue(type='direction')['numeric'][1]) # in deg
		nx= self.residual_im.shape()[0]
		ny= self.residual_im.shape()[1]
		print('INFO: (nx,ny)=(%s,%s), pixref_wcs(%s,%s)') % (nx,ny,xref_wcs,yref_wcs)

			
		#####################################
		## Read restored image from file
		#####################################
		file_ext= os.path.splitext(self.img_file)[1]
		if file_ext=='.fits':
			print 'INFO: Importing restored image FITS file %s as CASA image ...' % self.img_file
			self.img_im= ia.newimagefromfits(infile=self.img_file)
		else:
			print 'INFO: Importing restored CASA image %s ' % self.img_file
			self.img_im= ia.newimagefromfile(infile=self.img_file)

		## Get beam info from restored image
		beam= self.img_im.restoringbeam()
		units= beam['major']['unit']
		bmaj= beam['major']['value'] # expected in rad in CASA
		bmin= beam['minor']['value'] # expected in rad in CASA
		bpa= beam['positionangle']['value'] # expected in deg in CASA
		if units=='rad':
			bmaj= np.rad2deg(bmaj)	
			bmin= np.rad2deg(bmin)

		self.beam_bmaj= bmaj*3600 # in arcsec
		self.beam_bmin= bmin*3600 # in arcsec
		self.beam_bpa= bpa # in deg	
		print ('INFO: Restored image beam= (%s,%s,%s)') % (self.beam_bmaj,self.beam_bmin,self.beam_bpa)
		print 'INFO: Restored image shape=', self.img_im.shape()	
			
		#####################################
		## Read skymodel map from file
		#####################################
		print 'INFO: Importing model map from file %s as CASA image ' % self.model_file
		ia.open(self.model_file)

		## Create region from bounding box of beam residual map
		#print ('INFO: Creating region from bounding box of beam residual map [(x0,y0),(x1,y1)]=[(%s,%s)(%s,%s)]') % (x0_wcs,y0_wcs,x1_wcs,y1_wcs)
		#x0= ia.topixel([x0_wcs,y0_wcs])['numeric'][0]
		#y0= ia.topixel([x0_wcs,y0_wcs])['numeric'][1]
		#x1= ia.topixel([x1_wcs,y1_wcs])['numeric'][0]
		#y1= ia.topixel([x1_wcs,y1_wcs])['numeric'][1]
		xref= ia.topixel([str(xref_wcs)+'deg',str(yref_wcs)+'deg'])['numeric'][0]
		yref= ia.topixel([str(xref_wcs)+'deg',str(yref_wcs)+'deg'])['numeric'][1]
		x0= xref-nx/2
		x1= xref+nx/2-1
		y0= yref-ny/2
		y1= yref+ny/2-1
		print ('INFO: Creating region from bounding box of beam residual map pixref_wcs(%s,%s), pixref(%s,%s), [(x0,y0),(x1,y1)]=[(%s,%s)(%s,%s)]') % (xref_wcs,yref_wcs,xref,yref,x0,y0,x1,y1)
	
		reg= rg.box(blc=[x0,y0],trc=[x1,y1])
		
		#################################	
		## Extract model subimage
		#################################
		print ('INFO: Extracting model subimage using region [(x0,y0),(x1,y1)]=[(%s,%s)(%s,%s)]') % (x0,y0,x1,y1)
		self.model_im= ia.subimage(region=reg)
		self.model_data= self.model_im.getregion()
		print ('INFO: model units=%s') % self.model_im.brightnessunit()	
		print 'INFO: model shape=', self.model_im.shape()	

		## Close model image
		ia.done

		
		
	#####################################
	###         CONVOLVE SKYMODEL     ##
	#####################################	
	def convolve_skymodel(self):
		""" Convolve sky model with restored beam """
				
		## Convolving skymodel with restoring beam
		## NB: If input map is in Jy/pixel units, the convolve task automatically scale results to Jy/beam (no further scaling needed)
		bmaj= str(self.beam_bmaj) + 'arcsec'
		bmin= str(self.beam_bmin) + 'arcsec'
		bpa= str(self.beam_bpa) + 'deg'
		print ('INFO: Convolving skymodel with restored beam (%s,%s,%s)') % (bmaj,bmin,bpa)
		self.model_im_conv= self.model_im.convolve2d(type='gauss',major=bmaj,minor=bmin,pa=bpa)

		## Get pixel data from convolved image
		self.model_conv_data= self.model_im_conv.getregion()


	#####################################
	###         RUN                    ##
	#####################################
	def run(self):
		""" Convolve skymodel with beam """
		
		## == INITIALIZE DATA ==
		print ('INFO: Init and read input maps ...')
		self.init()
	
		## == CONVOLVE SKY MODEL WITH BEAM ==
		print ('INFO: Convolving skymodel with beam ...')
		self.convolve_skymodel()		

		
		## == MAKE FINAL MAP ==
		print ('INFO: Creating final map I=res + model_convolved ...')
		print 'INFO: residual data size=',self.residual_data.shape
		print 'INFO: model conv data size=',self.model_conv_data.shape	
		data= self.residual_data + self.model_conv_data
		cs = self.img_im.coordsys()
		bmaj= str(self.beam_bmaj) + 'arcsec'
		bmin= str(self.beam_bmin) + 'arcsec'
		bpa= str(self.beam_bpa) + 'deg'
		self.img_conv_im= ia.newimagefromarray(pixels=data, csys=cs.torecord(), outfile=self.img_outfile_casa, overwrite=True)
		self.img_conv_im.setbrightnessunit('Jy/beam')		
		self.img_conv_im.setrestoringbeam(major=bmaj,minor=bmin,pa=bpa)  
		self.img_conv_im.calcmask('T', name='mask1')
		self.img_conv_im.maskhandler('set', ['mask1'])

		## == WRITE MAPS TO FITS FILES ==
		print ('INFO: Writing images to FITS...')
		self.write_data_to_fits(self.img_conv_im,self.img_outfile)
		self.write_data_to_fits(self.model_im_conv,self.model_conv_outfile)

		## Close open CASA images
		self.model_im.done()
		self.img_im.done()
		self.residual_im.done()
		self.img_conv_im.done()



##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	#===========================
	#==   Get script args
	#===========================
	print('Get script args')
	try:
		args= get_args()
	except Exception as ex:
		print("Failed to get and parse options (err=%s)",str(ex))
		return 1
	
	# - Run options
	do_convolve= args.do_convolve	

	# - Input maps
	mosaic_file= args.mosaic
	residual_file= args.residual
	img_file= args.img
	model_file= args.model
	
	# - Compact source args
	marginX= args.marginx
	marginY= args.marginy
	Smin= args.Smin
	Smax= args.Smax
	source_density= args.source_density
	nsources= args.nsources

	bmaj= args.bmaj
	bmin= args.bmin
	pa= args.pa

	# - Extended source args
	add_ext_sources= args.extsources
	nsources_ext= args.nsources_ext
	source_density_ext= args.source_density_ext
	Smin_ext= args.Smin_ext
	Smax_ext= args.Smax_ext
	bmaj_min= args.bmaj_min
	bmaj_max= args.bmaj_max
	bmin_min= args.bmin_min
	bmin_max= args.bmin_max
	pa_min= args.pa_min
	pa_max= args.pa_max	

	# - Output args
	outfile_img= args.outfile_img
	outfile_sources= args.outfile_sources
	outfile_sourcepars= args.outfile_sourcepars
	outfile_ds9region= args.outfile_ds9region
	outfile_model= args.outfile_model
	outfile_modelconv= args.outfile_modelconv
	outfile_model_ext= args.outfile_model_ext

	print("*** ARGS ***")
	print("Margin X: %s" % marginX)
	print("Margin Y: %s" % marginY)
	print("Source flux range: (%s,%s)" % (Smin, Smax))
	print("Source density (deg^-2): %s" % source_density)	
	print("Image output filename : %s " % outfile_img)
	print("Beam (%s,%s,%s)" % (str(bmaj),str(bmin),str(pa)))
	print("Ext source flux range: (%s,%s)" % (Smin_ext, Smax_ext))
	print("Ext source density (deg^-2): %s" % source_density_ext)	
	print("Ext source gaus beam (%s-%s,%s-%s,%s-%s)" % (bmaj_min,bmaj_max,bmin_min,bmin_max,pa_min,pa_max))
	print("************")

	## Check if generate sky model
	if model_file:
		generate_model= False
		model_file_input= model_file
	else:
		generate_model= True
		model_file_input= outfile_model


	## Generate simulated sky map
	if generate_model:
		print ('INFO: Generate simulated skymodel...')
		simulator= SkyMapSimulator(mosaic_file)
		simulator.set_margins(marginX,marginY)
		simulator.set_model_outfile(outfile_model)
		simulator.set_ext_model_outfile(outfile_model_ext)
		simulator.set_img_outfile(outfile_img)
		simulator.set_ds9region_outfile(outfile_ds9region)
		simulator.set_source_outfile(outfile_sources)
		simulator.set_source_par_outfile(outfile_sourcepars)
		simulator.set_nsources(nsources)
		simulator.set_source_flux_range(Smin,Smax)
		simulator.set_source_density(source_density)
		simulator.set_beam(bmaj,bmin,pa)
		simulator.add_ext_sources(add_ext_sources)
		simulator.set_nsources_ext(nsources_ext)
		simulator.set_ext_source_flux_range(Smin_ext,Smax_ext)
		simulator.set_beam_bmaj_range(bmaj_min,bmaj_max)
		simulator.set_beam_bmin_range(bmin_min,bmin_max)
		simulator.set_beam_pa_range(pa_min,pa_max)
		simulator.run()

	## Generate restored image by convolving model with restored beam
	if do_convolve:
		print ('INFO: Convolve skymodel with beam ...')
		imager= SkyMapImager(residual_file,model_file_input,img_file)
		imager.set_map_outfile(outfile_img)
		imager.set_modelconv_outfile(outfile_modelconv)
		imager.run()
	
	print 'INFO: End run'
	
###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

