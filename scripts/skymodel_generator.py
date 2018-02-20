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
from scipy import ndimage
##import pyfits
from astropy.io import fits
from astropy.units import Quantity
from astropy.modeling.parameters import Parameter
from astropy.modeling.core import Fittable2DModel
from astropy.modeling.models import Box2D, Gaussian2D, Ring2D, Ellipse2D, TrapezoidDisk2D, Disk2D, AiryDisk2D, Sersic2D
from photutils.datasets import make_noise_image
from astropy import wcs

## ROOT
import ROOT
from ROOT import gSystem, TFile, TTree, gROOT, AddressOf

## CAESAR
gSystem.Load('libCaesar')
from ROOT import Caesar

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## Graphics modules
import matplotlib.pyplot as plt
import pylab
##################################################

#### GET SCRIPT ARGS ####
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
	
	# - GENERAL IMAGE OPTIONS
	parser.add_argument('-nx', '--nx', dest='nx', required=True, type=int, action='store',help='Image width in pixels')
	parser.add_argument('-ny', '--ny', dest='ny', required=True, type=int, action='store',help='Image height in pixels')
	parser.add_argument('-marginx', '--marginx', dest='marginx', required=False, type=int, default=0,action='store',help='Image x margin in pixels')
	parser.add_argument('-marginy', '--marginy', dest='marginy', required=False, type=int, default=0,action='store',help='Image y margin in pixels')
	
	parser.add_argument('-pixsize', '--pixsize', dest='pixsize', required=True, type=float, action='store',help='Map pixel size in arcsec')
	parser.add_argument('-bmaj', '--bmaj', dest='bmaj', required=True, type=float, default=10, action='store',help='Beam bmaj in arcsec (default=5)')
	parser.add_argument('-bmin', '--bmin', dest='bmin', required=True, type=float, default=5, action='store',help='Beam bmin in arcsec (default=5)')
	parser.add_argument('-bpa', '--bpa', dest='bpa', required=False, type=float, default=0, action='store',help='Beam bpa in deg (default=0)')
	parser.add_argument('-crpix1', '--crpix1', dest='crpix1', required=False, type=float, default=1, action='store',help='CRPIX1 fits keyword (default=1)')
	parser.add_argument('-crpix2', '--crpix2', dest='crpix2', required=False, type=float, default=1, action='store',help='CRPIX2 fits keyword (default=1)')
	parser.add_argument('-crval1', '--crval1', dest='crval1', required=False, type=float, default=254.851041667, action='store',help='CRVAL1 fits keyword (default=1)')
	parser.add_argument('-crval2', '--crval2', dest='crval2', required=False, type=float, default=-41.4765888889, action='store',help='CRVAL2 fits keyword (default=1)')
	#parser.add_argument('-ctype1', '--ctype1', dest='ctype1', required=False, type=str, default='RA---NCP', action='store',help='CTYPE1 fits keyword (default=1)')
	#parser.add_argument('-ctype2', '--ctype2', dest='ctype2', required=False, type=str, default='DEC--NCP', action='store',help='CTYPE2 fits keyword (default=1)')
	parser.add_argument('-ctype1', '--ctype1', dest='ctype1', required=False, type=str, default='RA---SIN', action='store',help='CTYPE1 fits keyword (default=1)')
	parser.add_argument('-ctype2', '--ctype2', dest='ctype2', required=False, type=str, default='DEC--SIN', action='store',help='CTYPE2 fits keyword (default=1)')

	# - BKG OPTIONS
	parser.add_argument('--bkg', dest='enable_bkg', action='store_true')	
	parser.add_argument('--no-bkg', dest='enable_bkg', action='store_false')	
	parser.set_defaults(enable_bkg=True)
	parser.add_argument('-bkg_level', '--bkg_level', dest='bkg_level', required=False, type=float, default=10e-6, action='store',help='Bkg level (default=0)')
	parser.add_argument('-bkg_rms', '--bkg_rms', dest='bkg_rms', required=False, type=float, default=100e-6, action='store',help='Bkg rms (default=0)')

	# - COMPACT SOURCE OPTIONS
	parser.add_argument('-npixels_min', '--npixels_min', dest='npixels_min', required=False, type=int, default=5, action='store',help='Minimum number of pixels for a generated source (default=5)')
	parser.add_argument('--compactsources', dest='enable_compactsources', action='store_true')	
	parser.add_argument('--no-compactsources', dest='enable_compactsources', action='store_false')	
	parser.set_defaults(enable_compactsources=True)
	parser.add_argument('-zmin', '--zmin', dest='zmin', required=False, type=float, default=1, action='store',help='Minimum source significance level in sigmas above the bkg (default=1)')
	parser.add_argument('-zmax', '--zmax', dest='zmax', required=False, type=float, default=30, action='store',help='Maximum source significance level in sigmas above the bkg (default=30)')
	parser.add_argument('-source_density', '--source_density', dest='source_density', required=False, type=float, default=1000, action='store',help='Compact source density (default=1000)')

	# - EXTENDED SOURCES
	parser.add_argument('--extsources', dest='enable_extsources', action='store_true')	
	parser.add_argument('--no-extsources', dest='enable_extsources', action='store_false')	
	parser.set_defaults(enable_extsources=True)
	parser.add_argument('-ext_source_density', '--ext_source_density', dest='ext_source_density', required=False, type=float, default=100, action='store',help='Extended source density (default=1000)')
	parser.add_argument('-zmin_ext', '--zmin_ext', dest='zmin_ext', required=False, type=float, default=0.1, action='store',help='Minimum extended source significance level in sigmas above the bkg (default=0.1)')
	parser.add_argument('-zmax_ext', '--zmax_ext', dest='zmax_ext', required=False, type=float, default=2, action='store',help='Maximum extended source significance level in sigmas above the bkg (default=2)')
	parser.add_argument('-ext_scale_min', '--ext_scale_min', dest='ext_scale_min', required=False, type=float, default=10, action='store',help='Minimum extended source size in arcsec (default=10)')
	parser.add_argument('-ext_scale_max', '--ext_scale_max', dest='ext_scale_max', required=False, type=float, default=3600, action='store',help='Maximum extended source size in arcsec (default=3600)')
	parser.add_argument('-ext_source_type', '--ext_source_type', dest='ext_source_type', required=False, type=int, default=-1, action='store',help='Extended source type to generate (-1=all types from available models, 1=ring, 2=ellipse, 3=bubble+shell, 4=airy disk (default=-1)')

	# - SOURCE MODEL OPTIONS
	parser.add_argument('-ring_rmin', '--ring_rmin', dest='ring_rmin', required=False, type=float, default=0.5, action='store',help='Minimum ring radius in arcsec (default=1)')
	parser.add_argument('-ring_rmax', '--ring_rmax', dest='ring_rmax', required=False, type=float, default=10, action='store',help='Maximum ring radius in arcsec (default=10)')
	parser.add_argument('-ring_wmin', '--ring_wmin', dest='ring_wmin', required=False, type=float, default=5, action='store',help='Minimum ring width in arcsec (default=1)')
	parser.add_argument('-ring_wmax', '--ring_wmax', dest='ring_wmax', required=False, type=float, default=20, action='store',help='Maximum ring width in arcsec (default=10)')
	parser.add_argument('-ellipse_rmin', '--ellipse_rmin', dest='ellipse_rmin', required=False, type=float, default=0.5, action='store',help='Ellipse bmaj in arcsec (default=1)')
	parser.add_argument('-ellipse_rmax', '--ellipse_rmax', dest='ellipse_rmax', required=False, type=float, default=10, action='store',help='Ellipse bmin in arcsec (default=10)')
	
	parser.add_argument('-disk_shell_ampl_ratio_min', '--disk_shell_ampl_ratio_min', dest='disk_shell_ampl_ratio_min', required=False, type=float, default=0.1, action='store',help='Disk/shell amplitude ratio min (default=0.1)')
	parser.add_argument('-disk_shell_ampl_ratio_max', '--disk_shell_ampl_ratio_max', dest='disk_shell_ampl_ratio_max', required=False, type=float, default=0.8, action='store',help='Disk/shell amplitude ratio max (default=0.8)')
	parser.add_argument('-disk_shell_radius_ratio_min', '--disk_shell_radius_ratio_min', dest='disk_shell_radius_ratio_min', required=False, type=float, default=0.6, action='store',help='Disk/shell radius ratio min (default=0.6)')
	parser.add_argument('-disk_shell_radius_ratio_max', '--disk_shell_radius_ratio_max', dest='disk_shell_radius_ratio_max', required=False, type=float, default=0.9, action='store',help='Disk/shell radius ratio max (default=0.8)')
		
	parser.add_argument('-zmin_model', '--zmin_model', dest='model_trunc_zmin', required=False, type=float, default=1, action='store',help='Minimum source significance level in sigmas above the bkg below which source data are set to 0 (default=1)')
	parser.add_argument('-mask_boxsize', '--mask_boxsize', dest='mask_boxsize', required=False, type=float, default=10, action='store',help='Mask box size in pixels (default=10)')
	parser.add_argument('-trunc_thr', '--trunc_thr', dest='trunc_thr', required=False, type=float, default=0.01, action='store',help='Source model truncation thr (default=0.01)')

	# - OUTPUT FILE OPTIONS
	parser.add_argument('-outputfile', '--outputfile', dest='outputfile', required=False, type=str, default='simmap.fits',action='store',help='Output filename')
	parser.add_argument('-outputfile_model', '--outputfile_model', dest='outputfile_model', required=False, type=str, default='skymodel.fits', action='store',help='Model filename')
	parser.add_argument('-outputfile_sources', '--outputfile_sources', dest='outputfile_sources', required=False, type=str, default='sources.root',action='store',help='Skymodel source ROOT Output filename')
	parser.add_argument('-outputfile_ds9region', '--outputfile_ds9region', dest='outputfile_ds9region', required=False, type=str, default='dsregion.reg',action='store',help='DS9 source region filename')
	parser.add_argument('-outputfile_casaregion', '--outputfile_casaregion', dest='outputfile_casaregion', required=False, type=str, default='casa_mask.dat',action='store',help='CASA source region filename')
	
	args = parser.parse_args()	

	return args





###########################
##     MODELS
###########################

class RingSector2D(Fittable2DModel):
	""" Two dimensional radial symmetric Ring model """

	amplitude = Parameter(default=1)
	x_0 = Parameter(default=0)
	y_0 = Parameter(default=0)
	r_in = Parameter(default=1)
	width = Parameter(default=1)
	theta_min = Parameter(default=-np.pi)
	theta_max = Parameter(default=np.pi)

	def __init__(self, amplitude=amplitude.default, x_0=x_0.default, y_0=y_0.default, r_in=r_in.default, width=width.default, theta_min=theta_min.default, theta_max=theta_max.default, **kwargs):
		# If outer radius explicitly given, it overrides default width.
		if width is None:
			width = self.width.default
		if theta_min is None:
			theta_min = self.theta_min.default
		if theta_max is None:
			theta_max = self.theta_max.default

		super(RingSector2D, self).__init__(amplitude=amplitude, x_0=x_0, y_0=y_0, r_in=r_in, width=width, theta_min=theta_min, theta_max=theta_max, **kwargs)

	@staticmethod
	def evaluate(x, y, amplitude, x_0, y_0, r_in, width, theta_min, theta_max):
		"""Two dimensional Ring sector model function."""

		rr = (x - x_0) ** 2 + (y - y_0) ** 2
		theta = np.arctan2(x-x_0,y-y_0)
		r_range = np.logical_and(rr >= r_in ** 2, rr <= (r_in + width) ** 2)
		theta_range= np.logical_and(theta>=theta_min, theta<=theta_max)
		sector_range = np.logical_and(r_range,theta_range)
		result = np.select([sector_range], [amplitude])

		if isinstance(amplitude, Quantity):
			return Quantity(result, unit=amplitude.unit, copy=False)
		else:
			return result

	@property
	def bounding_box(self):
		"""
			Tuple defining the default ``bounding_box``.
				``((y_low, y_high), (x_low, x_high))``
		"""

		dr = self.r_in + self.width
		return ( (self.y_0 - dr, self.y_0 + dr), (self.x_0 - dr, self.x_0 + dr) )

	@property
	def input_units(self):
		if self.x_0.unit is None:
			return None
		else:
			return {'x': self.x_0.unit, 'y': self.y_0.unit}

	def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
		# Note that here we need to make sure that x and y are in the same
		# units otherwise this can lead to issues since rotation is not well
		# defined.
		if inputs_unit['x'] != inputs_unit['y']:
			raise UnitsError("Units of 'x' and 'y' inputs should match")
		return OrderedDict([('x_0', inputs_unit['x']), ('y_0', inputs_unit['x']), ('r_in', inputs_unit['x']), ('width', inputs_unit['x']), ('amplitude', outputs_unit['z'])])


###########################
##     SIMULATOR CLASS
###########################

class SkyMapSimulator(object):

	""" Sky map simulator class

			Attributes:
				nx: image width in pixels
				ny: image height in pixels
				pixsize: pixel size in arcsec (default=1)
	"""

	def __init__(self, nx, ny, pixsize=1):
		""" Return a SkyMapGenerator object """

		## Image parameters
		self.nx = nx #in pixels
		self.ny = ny # in pixels
		self.marginx= 0 # in pixels (no margin)
		self.marginy= 0 # in pixels (no margin)
		
		self.pixsize= pixsize # in arcsec
		self.gridy, self.gridx = np.mgrid[0:ny, 0:nx]
		self.crpix1= 1
		self.crpix2= 1
		self.crval1= 254.851041667
		self.crval2= -41.4765888889
		self.ctype1= 'RA---SIN'
		self.ctype2= 'DEC--SIN'

		## Source model
		self.truncate_models= True
		self.trunc_thr= 0.01 # 1% flux truncation at maximum
		self.trunc_model_zmin= 1

		## Mask box size
		self.mask_boxsize= 10 # in pixels

		## Bkg parameters
		self.simulate_bkg= True
		self.bkg_level= 0 # in Jy
		self.bkg_rms= 10.e-6 # in Jy

		## Compact source parameters
		self.simulate_compact_sources= True
		self.source_density= 2000. # in sources/deg^2
		self.beam_bmaj= 6.5 # in arcsec
		self.beam_bmin= 6.5 # in arcsec
		self.beam_bpa= 0 # in deg
		self.beam_area= self.compute_beam_area(self.beam_bmaj,self.beam_bmin) # in pixels
		self.zmin= 1 # in sigmas 
		self.zmax= 30 # in sigmas
		self.npixels_min= 5
		
		## Extended source parameters
		self.simulate_ext_sources= True
		self.ext_source_type= -1 # all source models generated
		self.ext_source_density= 10 # in sources/deg^2
		self.zmin_ext= 0.5 # in sigmas 
		self.zmax_ext= 5	 # in sigmas 
		self.ring_rmin= 2. # in arcsec
		self.ring_rmax= 10. # in arcsec
		self.ring_width_min= 5 # in arcsec
		self.ring_width_max= 10 # in arcsec 	 
		self.ellipse_rmin= 1 # in arcsec
		self.ellipse_rmax= 10 # in arcsec
		self.ellipse_rratiomin= 0.7 # ratio rmin/rmax
		self.disk_rmin= 2 # in arcsec
		self.disk_rmax= 10 # in arcsec
		self.shell_disk_ampl_ratio_min= 0.1
		self.shell_disk_ampl_ratio_max= 0.8
		self.shell_disk_radius_ratio_min= 0.6
		self.shell_disk_radius_ratio_max= 0.9
		self.sersic_radius= 10 # in arcsec
		self.sersic_ellipticity= 0.5
		self.sersic_index= 4		

		## Map output file
		self.mapfilename= 'simmap.fits'
		self.modelfilename= 'skymodel.fits'
		
		## DS9 output file
		self.ds9filename= 'ds9region.reg'

		## CASA region output file	
		self.casafilename= 'casamask.dat'

		## Caesar img & sources
		self.outfilename= 'SimOutput.root'	
		self.outfile= None
		self.outtree= None
		self.cs = None
		self.caesar_sources= []	
		self.caesar_img= None

	def init(self):
		""" Initialize data """
		## Initialize output tree & file
		self.outfile= ROOT.TFile(self.outfilename,'RECREATE')
		self.outtree= ROOT.TTree('SourceInfo','SourceInfo')
		self.cs = Caesar.Source()
		self.outtree.Branch('Source',self.cs)		

	def set_mask_box_size(self,boxsize):
		""" Set mask box size """
		if boxsize<=0:	
			raise ValueError('Invalid boxsize specified (shall be larger than 0')
		self.mask_boxsize= boxsize	

	def set_margins(self,marginx,marginy):
		""" Set margin in X & Y """
		if (marginx<0 or marginy<0 or marginx>=self.nx/2 or marginy>=self.ny/2) :
			raise ValueError('Invalid margin specified (<0 or larger than image half size!')
		self.marginx= marginx
		self.marginy= marginy

	def set_ref_pix(self,x,y):
		""" Set reference pixel (CRPIX1,CRPIX2) in FITS output """
		self.crpix1= x
		self.crpix2= y

	def set_ref_pix_coords(self,x,y):
		""" Set reference pixel coords (CRPIX1,CRPIX2) in FITS output """
		self.crval1= x
		self.crval2= y

	def set_coord_system_type(self,x,y):
		""" Set coord system type (CTYPE1,CTYPE2) in FITS output """
		self.ctype1= x
		self.ctype2= y
	
	def enable_compact_sources(self,choice):
		""" Enable/disable compact source generation """
		self.simulate_compact_sources= choice

	def enable_extended_sources(self,choice):
		""" Enable/disable extended source generation """
		self.simulate_extended_sources= choice

	def enable_bkg(self,choice):
		""" Enable/disable bkg generation """
		self.simulate_bkg= choice

	def set_npixels_min(self,value):
		""" Set the minimum number of pixels for a generated source"""
		self.npixels_min= value

	def truncate_models(self,choice):
		""" Enable/disable continuous model truncation (gaussian, airy disk, ...) """
		self.truncate_models= choice

	def set_model_trunc_significance(self,value):
		""" Set the significance level below which source model data are truncated """
		self.trunc_model_zmin= value

	def set_model_trunc_thr(self,value):
		""" Set the flux percentage level for source model truncation """
		self.trunc_thr= value

	def set_ext_source_type(self,value):
		""" Set the extended source type to be generated (-1=all, 1=ring, 2=ellipse, 3=bubble+shell, 4=airy)"""
		self.ext_source_type= value

	def set_ds9region_filename(self,filename):
		""" Set the output DS9 region filename """
		self.ds9filename= filename

	def set_casaregion_filename(self,filename):
		""" Set the output CASA region filename """
		self.casafilename= filename

	def set_map_filename(self,filename):
		""" Set the output map filename """
		self.mapfilename= filename

	def set_model_filename(self,filename):
		""" Set the output model filename """
		self.modelfilename= filename

	def set_source_filename(self,filename):
		""" Set the output source ROOT filename """
		self.outfilename= filename

	def set_source_significance_range(self,zmin,zmax):
		""" Set source significance range """
		self.zmin= zmin
		self.zmax= zmax

	def set_ext_source_significance_range(self,zmin,zmax):
		""" Set source significance range """
		self.zmin_ext= zmin
		self.zmax_ext= zmax

	def set_source_density(self,density):
		""" Set compact source density in deg^-2 """
		self.source_density= density

	def set_ext_source_density(self,density):
		""" Set extended source density in deg^-2 """
		self.ext_source_density= density

	def set_ring_pars(self,rmin,rmax,wmin,wmax):
		""" Set ring model parameters"""
		self.ring_rmin= rmin
		self.ring_rmax= rmax
		self.ring_width_min= wmin
		self.ring_width_max= wmax 	 

	def set_sersic_pars(self,radius,ell,index):
		""" Set Sersic model pars"""
		self.sersic_radius= radius
		self.sersis_ellipticity= ell
		self.sersic_index= index

	def set_disk_pars(self,rmin,rmax):
		""" Set disk model parameters"""
		self.disk_rmin= rmin
		self.disk_rmax= rmax
		
	def set_disk_shell_pars(self,ampl_ratio_min,ampl_ratio_max,radius_ratio_min,radius_ratio_max):
		""" Set disk shell model parameters"""
		self.shell_disk_ampl_ratio_min= ampl_ratio_min
		self.shell_disk_ampl_ratio_max= ampl_ratio_max
		self.shell_disk_radius_ratio_min= radius_ratio_min
		self.shell_disk_radius_ratio_max= radius_ratio_max
		

	def set_ellipse_pars(self,rmin,rmax):
		""" Set ring model parameters"""
		self.ellipse_rmin= rmin
		self.ellipse_rmax= rmax

	def set_bkg_pars(self,bkg_level,bkg_rms):	
		""" Set bkg parameters """
		self.bkg_level= bkg_level
		self.bkg_rms= bkg_rms	

	def set_beam_info(self,Bmaj,Bmin,Bpa):
		""" Set beam info """
		self.beam_bmaj= Bmaj
		self.beam_bmin= Bmin
		self.beam_bpa= Bpa
		self.beam_area= self.compute_beam_area(Bmaj,Bmin)
		
	def compute_beam_area(self,Bmaj,Bmin):
		""" Compute beam area """
		A= np.pi*Bmaj*Bmin/(4*np.log(2)) #2d gaussian area with FWHM=fx,fy (in arcsec^2)
		pixelArea= np.fabs(self.pixsize*self.pixsize) # in arcsec^2
		beam_area= A/pixelArea # in pixels
		return beam_area

	def compute_beam_sigma(self,fwhm):
		""" """
		sigma= fwhm/(2.*np.sqrt(2.*np.log(2.)))
		return sigma

	def generate_bkg(self):
		""" Generate bkg data """
		shape = (self.ny, self.nx)
		bkg_data = make_noise_image(shape, type='gaussian', mean=self.bkg_level, stddev=self.bkg_rms)
		return bkg_data
    
	def generate_blob(self,ampl,x0,y0,sigmax,sigmay,theta,trunc_thr=0.01):
		""" Generate a blob 
				Arguments: 
					ampl: peak flux in Jy
					x0, y0: gaussian means in pixels
					sigmax, sigmay: gaussian sigmas in pixels
					theta: rotation in degrees
					trunc_thr: truncation significance threshold
		"""
		#modelFcn= Gaussian2D(ampl,x0,y0,sigmax,sigmay,theta=math.radians(theta)) 
		data= Gaussian2D(ampl,x0,y0,sigmax,sigmay,theta=math.radians(theta))(self.gridx, self.gridy)
		totFlux= (float)(np.sum(data,axis=None))
		print('totFlux=%s' % str(totFlux))

		## Truncate data such that sum(data)_trunc/sum(data)<f
		f= trunc_thr 
		if self.truncate_models:
			data_vect_sorted= np.ravel(data)
			data_csum= np.cumsum(data_vect_sorted)/totFlux
			fluxThr= data_vect_sorted[np.argmin(data_csum<f)]
			print('fluxThr=%s' % str(fluxThr))
			data[data<fluxThr] = 0		

		## Truncate data at minimum significance
		#ampl_min= (trunc_thr*self.bkg_rms) + self.bkg_level
		#if self.truncate_models:
		#	data[data<ampl_min] = 0		

		return data

	def generate_ring(self,ampl,x0,y0,radius,width):
		""" Generate a ring 
				Arguments: 
					ampl: peak flux in Jy
					x0, y0: means in pixels
					radius: ring radius in pixels
					width: ring width in pixels
		"""
		data= Ring2D(ampl,x0,y0,radius,width)(self.gridx, self.gridy) 
		return data

	def generate_ring_sector(self,ampl,x0,y0,radius,width,theta_min,theta_max):
		""" Generate a ring 
				Arguments: 
					ampl: peak flux in Jy
					x0, y0: means in pixels
					radius: ring radius in pixels
					width: ring width in pixels
					theta_min, theta_max: sector theta min/max in degrees
		"""
		data= RingSector2D(ampl,x0,y0,radius,width,np.radians(theta_min),np.radians(theta_max))(self.gridx, self.gridy) 
		return data

	def generate_bubble(self,ampl,x0,y0,radius,shell_ampl,shell_radius,shell_width,shell_theta_min,shell_theta_max):
		""" Generate a bubble with a shell """
		disk_data= Disk2D(ampl,x0,y0,radius)(self.gridx, self.gridy) 
		shell_data= self.generate_ring_sector(shell_ampl,x0,y0,shell_radius,shell_width,shell_theta_min,shell_theta_max)
		data= disk_data + shell_data
		return data
		

	def generate_ellipse(self,ampl,x0,y0,a,b,theta):
		""" Generate ellipse """
		data= Ellipse2D(ampl,x0,y0,a,b,math.radians(theta))(self.gridx, self.gridy) 
		return data

	def generate_airy_disk(self,ampl,x0,y0,radius,trunc_thr=0.01):
		""" Generate Airy disk """
		data= AiryDisk2D(amplitude=ampl,x_0=x0,y_0=y0,radius=radius)(self.gridx, self.gridy) 
		totFlux= (float)(np.sum(data,axis=None))

		## Truncate data such that sum(data)_trunc/sum(data)<f
		f= trunc_thr 
		if self.truncate_models:
			data_vect_sorted= np.ravel(data)
			data_csum= np.cumsum(data_vect_sorted)/totFlux
			fluxThr= data_vect_sorted[np.argmin(data_csum<f)]
			data[data<fluxThr] = 0

		## Truncate data at minimum significance
		#ampl_min= (self.zmin_ext*self.bkg_rms) + self.bkg_level
		#if self.truncate_models:
		#	data[data<ampl_min] = 0			

		return data

	def generate_sersic(self,ampl,x0,y0,radius,ell,index,theta,trunc_thr=0.01):
		""" Generate Sersic model """
		data= Sersic2D(amplitude=ampl,x_0=x0,y_0=y0,r_eff=radius,n=index,ellip=ell,theta=math.radians(theta))(self.gridx, self.gridy) 
		totFlux= (float)(np.sum(data,axis=None))

		## Truncate data such that sum(data)_trunc/sum(data)<f
		f= trunc_thr 
		if self.truncate_models:
			data_vect_sorted= np.ravel(data)
			data_csum= np.cumsum(data_vect_sorted)/totFlux
			fluxThr= data_vect_sorted[np.argmin(data_csum<f)]
			data[data<fluxThr] = 0	

		## Truncate data at minimum significance
		#ampl_min= (self.zmin_ext*self.bkg_rms) + self.bkg_level
		#if self.truncate_models:
		#	data[data<ampl_min] = 0			

		return data

	def make_caesar_source(self,source_data,source_name,source_id,source_type,source_sim_type,ampl=None,x0=None,y0=None,source_max_scale=None):
		""" Create Caesar source from source data array """

		# Create Caesar source
		source= Caesar.Source()

		# Get source indexes and fill pixels in Caesar source
		source_indexes= np.column_stack(np.where(source_data!=0))
		nRows= (source_data.shape)[0]	
		nCols= (source_data.shape)[1]	
		for index in source_indexes:
			rowId= index[0]
			colId= index[1]
			S= source_data[rowId,colId]
			ix= colId
			iy= rowId
			#iy= nRows-1-rowId
			gbin= ix + iy*nCols
			pixel= Caesar.Pixel(gbin,ix,iy,ix,iy,S)
			source.AddPixel(pixel)

			# Is at edge
			if (ix==0) or (ix==nCols-1) or (iy==0) or (iy==nRows-1):
				source.SetEdgeFlag(True)

		# Retun None if npixels is too small	
		nPix= source.GetNPixels()
		if nPix<self.npixels_min:
			print ('INFO: Too few pixels (%s) for this source, return None!' % str(nPix))
			return None
			
		# If true info are not given compute them
		#    - S= count integral
		#    - baricenter of binary map
		if x0 is None or y0 is None:
			print ('INFO: No source true pos given, computing it from data...')
			data_binary= np.where(source_data!=0,1,0)
			#print 'INFO: data_binary sum=%d', sum(data_binary)
			[y0,x0]= ndimage.measurements.center_of_mass(data_binary)

		if ampl is None:
			print ('INFO: No source true flux given, computing integral from data...')
			ampl= np.sum(source_data,axis=None)

		# Set some flags
		source.SetName(source_name)
		source.SetId(source_id)
		source.SetType(source_type)
		source.SetFlag(Caesar.Source.eFake)
		source.SetSimType(source_sim_type)
		if source_max_scale is not None:
			source.SetSimMaxScale(source_max_scale)
		source.SetTrueInfo(ampl,x0,y0)

		# Set flux correction factor
		fluxCorrection= self.beam_area
		source.SetBeamFluxIntegral(fluxCorrection)

		# Compute stats & morph pars
		source.ComputeStats();
		source.ComputeMorphologyParams();

		return source

	def make_caesar_image(self,data):
		""" Make Caesar image from array data """
		
		# Get source indexes and fill pixels in Caesar source
		img_indexes= np.column_stack(np.where(data!=0))
		nRows= (data.shape)[0]	
		nCols= (data.shape)[1]	

		# Set metadata
		metadata= Caesar.ImgMetaData()
		metadata.Nx= self.nx
		metadata.Ny= self.ny
		metadata.Cx= (int)(self.crpix1)
		metadata.Cy= (int)(self.crpix2)
		metadata.Xc= self.crval1
		metadata.Yc= self.crval2
		metadata.dX= -self.pixsize/3600.
		metadata.dY= self.pixsize/3600.
		metadata.CoordTypeX= self.ctype1
		metadata.CoordTypeY= self.ctype2
		metadata.BUnit= 'JY/PIXEL'
		metadata.Bmaj= self.beam_bmaj/3600.
		metadata.Bmin= self.beam_bmin/3600.
		metadata.Bpa= self.beam_bpa

		# Create Caesar image
		img= Caesar.Image(nCols,nRows,"img")
		img.SetMetaData(metadata)
		
		for index in img_indexes:
			rowId= index[0]
			colId= index[1]
			S= data[rowId,colId]
			ix= colId
			iy= rowId
			#iy= nRows-1-rowId
			gbin= ix + iy*nCols
			img.FillPixel(ix,iy,S,True);

		return img

	def generate_compact_sources(self):
		""" Generate list of compact sources in the map.
					- Uniform spatial distribution
					- Uniform flux distribution

				Arguments:	
					density: source density in #sources/deg^2 (e.g. 2000)
		"""
		# Compute number of sources to be generated given map area in pixels
		#area= (self.nx*self.ny)*self.pixsize/(3600.*3600.) # in deg^2
		area= ((self.nx-2*self.marginx)*(self.ny-2*self.marginy))*self.pixsize/(3600.*3600.) # in deg^2
		nsources= int(round(self.source_density*area))
		S_min= (self.zmin*self.bkg_rms) + self.bkg_level
		S_max= (self.zmax*self.bkg_rms) + self.bkg_level
		lgS_min= np.log(S_min)
		lgS_max= np.log(S_max)
		randomize_flux= False
		if self.zmin<self.zmax:
			randomize_flux= True

		print 'INFO: Generating #',nsources,' compact sources in map...'

		# Compute blob sigma pars given beam info
		sigmax= self.compute_beam_sigma(self.beam_bmaj)
		sigmay= self.compute_beam_sigma(self.beam_bmin)
		theta= self.beam_bpa + 90. # NB: BPA is the positional angle of the major axis measuring from North (up) counter clockwise, while theta is measured wrt to x axis
		source_max_scale= 2*max(self.beam_bmaj,self.beam_bmin)		

		## Start generation loop
		sources_data = Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.nx, y_width=2*self.ny)(self.gridx, self.gridy)
		mask_data = Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.nx, y_width=2*self.ny)(self.gridx, self.gridy)
		
		for index in range(0,nsources):	
			if index%100==0 :
				print ("INFO: Generating compact source no. %s/%s" % (index+1,nsources))
			
			## Generate random coordinates
			#x0= np.random.uniform(0,self.nx)
			#y0= np.random.uniform(0,self.ny)
			#x0= np.random.uniform(0,self.nx-1)
			#y0= np.random.uniform(0,self.ny-1)
			x0= np.random.uniform(self.marginx,self.nx-self.marginx-1)
			y0= np.random.uniform(self.marginy,self.ny-self.marginy-1)

			## Compute amplitude given significance level and bkg
			## Generate flux uniform in log
			if randomize_flux:
				lgS= np.random.uniform(lgS_min,lgS_max)
				S= np.exp(lgS)
				z= (S-self.bkg_level)/self.bkg_rms
			else:
				S= (self.zmin*self.bkg_rms) + self.bkg_level
				z= self.zmin
	
			## Generate blob
			#blob_data= self.generate_blob(ampl=S,x0=x0,y0=y0,sigmax=sigmax/self.pixsize,sigmay=sigmay/self.pixsize,theta=theta,trunc_thr=self.trunc_model_zmin)
			blob_data= self.generate_blob(ampl=S,x0=x0,y0=y0,sigmax=sigmax/self.pixsize,sigmay=sigmay/self.pixsize,theta=theta,trunc_thr=self.trunc_thr)
			if blob_data is None:
				print('Failed to generate blob (hint: too large trunc threshold), skip and regenerate...')
				continue

			sources_data+= blob_data

			## Set model map
			#ix= int(np.floor(x0))
			#iy= int(np.floor(y0))
			ix= int(np.round(x0))
			iy= int(np.round(y0))
			mask_data[iy,ix]+= S

			# Make Caesar source	
			source_name= 'S' + str(index+1)
			source_id= index+1
			source_type= Caesar.Source.ePointLike
			caesar_source= self.make_caesar_source(blob_data,source_name,source_id,source_type,Caesar.Source.eBlobLike,ampl=S,x0=x0,y0=y0,source_max_scale=source_max_scale)
			if caesar_source is None:
				print('Generate source has too few pixels, skip and regenerate...')
				continue
				
			self.caesar_sources.append(caesar_source)

			print ('INFO: Source %s: Pos(%s,%s), ix=%s, iy=%s, S=%s' % (source_name,str(x0),str(y0),str(ix),str(iy),str(S)))

		return [sources_data,mask_data]


	def generate_extended_sources(self):
		""" Generate list of extended sources in the map.
					- Uniform spatial distribution
					- Uniform flux distribution

				Arguments:	
					density: source density in #sources/deg^2 (e.g. 2000)
		"""
		
		# Compute number of sources to be generated given map area in pixels
		#area= (self.nx*self.ny)*self.pixsize/(3600.*3600.) # in deg^2
		area= ((self.nx-2*self.marginx)*(self.ny-2*self.marginy))*self.pixsize/(3600.*3600.) # in deg^2
		nsources= int(round(self.ext_source_density*area))
		S_min= (self.zmin_ext*self.bkg_rms) + self.bkg_level
		S_max= (self.zmax_ext*self.bkg_rms) + self.bkg_level
		lgS_min= np.log(S_min)
		lgS_max= np.log(S_max)
		randomize_flux= False
		if self.zmin_ext<self.zmax_ext:
			randomize_flux= True

		print 'INFO: Generating #',nsources,' extended sources in map...'

		print('INFO: zmin_ext=%s, zmax_ext=%s, Smin=%s, Smax=%s' % (str(self.zmin_ext),str(self.zmax_ext),str(S_min),str(S_max)) )

		
		## Start generation loop
		sources_data = Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.nx, y_width=2*self.ny)(self.gridx, self.gridy)
		ngen_sources= 0	
		if self.ext_source_type==-1:	
			nsource_types= 5
		else:
			nsource_types= 1

		#for index in range(0,nsources):	
		while (ngen_sources<nsources):
			if ngen_sources%10==0 :
				print ("INFO: Generating extended source no. %s/%s" % (ngen_sources+1,nsources))
			
			## Generate random coordinates
			#x0= random.uniform(0,self.nx)
			#y0= random.uniform(0,self.ny)
			#x0= np.random.uniform(0,self.nx-1)
			#y0= np.random.uniform(0,self.ny-1)
			x0= np.random.uniform(self.marginx,self.nx-self.marginx-1)
			y0= np.random.uniform(self.marginy,self.ny-self.marginy-1)

			## Compute amplitude given significance level and bkg
			## Generate flux uniform in log
			if randomize_flux:
				lgS= np.random.uniform(lgS_min,lgS_max)
				S= np.exp(lgS)
				z= (S-self.bkg_level)/self.bkg_rms
				#z= random.uniform(self.zmin_ext,self.zmax_ext)
				#S= (z*self.bkg_rms) + self.bkg_level
			else:
				S= (self.zmin_ext*self.bkg_rms) + self.bkg_level
				z= self.zmin_ext 			

			## Generate random type (1=ring, 2=ellipse, ...)
			if self.ext_source_type==-1:
				source_sim_type= random.randint(1, nsource_types)
			else:
				source_sim_type= self.ext_source_type

			source_max_scale= 0.	

			if source_sim_type==1: # Ring2D Sector model
				source_sim_type= Caesar.Source.eRingLike
				ring_r= random.uniform(self.ring_rmin,self.ring_rmax) 
				ring_w= random.uniform(self.ring_width_min,self.ring_width_max)
				#source_data= self.generate_ring(S,x0,y0,ring_r/self.pixsize,ring_w/self.pixsize) # convert radius/width from arcsec to pixels		
				theta1= random.uniform(-180,180)
				theta2= random.uniform(-180,180)
				theta_min= min(theta1,theta2)
				theta_max= max(theta1,theta2)
				dtheta= theta_max-theta_min
				r= ring_r
				R= ring_r + ring_w
				sector_diagonal= np.sqrt( r*r + R*R - 2*r*R*np.cos(np.deg2rad(dtheta)) )
				sector_arc= 2*R*np.pi*dtheta/360.
				source_max_scale= max(max(sector_arc,ring_w),sector_diagonal)
				source_data= self.generate_ring_sector(S,x0,y0,ring_r/self.pixsize,ring_w/self.pixsize,theta_min,theta_max) # convert radius/width from arcsec to pixels
				
			elif source_sim_type==2: # Ellipse 2D model
				source_sim_type= Caesar.Source.eEllipseLike
				ellipse_bmaj= random.uniform(self.ellipse_rmin,self.ellipse_rmax) 
				#ellipse_bmin= random.uniform(self.ellipse_rmin,self.ellipse_rmax)
				ellipse_bmin= random.uniform(max(self.ellipse_rratiomin*ellipse_bmaj,self.ellipse_rmin),self.ellipse_rmax)
				ellipse_theta= random.uniform(0,360)
				source_max_scale= max(ellipse_bmaj,ellipse_bmin)
				source_data= self.generate_ellipse(S,x0,y0,ellipse_bmaj/self.pixsize,ellipse_bmin/self.pixsize,ellipse_theta) # convert radius/width from arcsec to pixels

			elif source_sim_type==3: # bubble + shell model
				source_sim_type= Caesar.Source.eBubbleLike
				bubble_r= random.uniform(self.disk_rmin,self.disk_rmax) 
				shell_excess= random.uniform(self.shell_disk_ampl_ratio_min,self.shell_disk_ampl_ratio_max)
				shell_S= S*(1+shell_excess)
				shell_r= random.uniform(bubble_r*self.shell_disk_radius_ratio_min,bubble_r*self.shell_disk_radius_ratio_max)
				shell_width= random.uniform(0,bubble_r-shell_r)
				theta1= random.uniform(-180,180)
				theta2= random.uniform(-180,180)
				theta_min= min(theta1,theta2)
				theta_max= max(theta1,theta2)
				source_max_scale= bubble_r*2
				source_data= self.generate_bubble(S,x0,y0,bubble_r,shell_S,shell_r,shell_width,theta_min,theta_max)
				
			#elif source_sim_type==4: # Airy disk
			#	source_sim_type= Caesar.Source.eDiskLike
			#	disk_r= random.uniform(self.disk_rmin,self.disk_rmax) 
			#	source_data= self.generate_airy_disk(S,x0,y0,disk_r)

			elif source_sim_type==4: # Sersic
				source_sim_type= Caesar.Source.eDiskLike
				sersic_r= random.uniform(self.disk_rmin,self.disk_rmax)
				sersic_theta= random.uniform(0,360)
				sersic_ell= random.uniform(0,1)
				source_max_scale= 2*sersic_r
				##source_data= self.generate_sersic(S,x0,y0,sersic_r,sersic_ell,self.sersic_index,sersic_theta)
				source_data= self.generate_sersic(S,x0,y0,sersic_r,sersic_ell,self.sersic_index,sersic_theta,trunc_thr=self.trunc_thr)

			elif source_sim_type==5: # Gaussian Blob like
				source_sim_type= Caesar.Source.eBlobLike
				blob_bmaj= random.uniform(self.ellipse_rmin,self.ellipse_rmax)
				#blob_bmin= random.uniform(self.ellipse_rmin,self.ellipse_rmax)
				blob_bmin= random.uniform(max(self.ellipse_rratiomin*blob_bmaj,self.ellipse_rmin),blob_bmaj)
				blob_theta= random.uniform(0,360)
				source_max_scale= 2*max(blob_bmin,blob_bmaj)
				#source_data= self.generate_blob(ampl=S,x0=x0,y0=y0,sigmax=blob_bmaj/self.pixsize,sigmay=blob_bmin/self.pixsize,theta=blob_theta,trunc_thr=self.zmin_ext)
				source_data= self.generate_blob(ampl=S,x0=x0,y0=y0,sigmax=blob_bmaj/self.pixsize,sigmay=blob_bmin/self.pixsize,theta=blob_theta,trunc_thr=self.trunc_thr)
				if source_data is None:
					print('Failed to generate blob (hint: too large trunc threshold), skip and regenerate...')
					continue

			else:
				print('ERROR: Invalid source type given!')
				continue			
	
			## Check if source data contains all zeros (e.g. truncation removed all data)
			if np.count_nonzero(source_data)<=0:
				print('WARN: Generated extended source data contains all zeros, regenerate...')
				continue

			## Check if source pixels and its contour has been already taken before
			source_indexes= (source_data!=0) # get all source data pixels (others are 0)
			source_indexes_xright= (np.roll(source_data,1,axis=1)!=0)
			source_indexes_xleft= (np.roll(source_data,-1,axis=1)!=0)
			source_indexes_yright= (np.roll(source_data,1,axis=0)!=0)
			source_indexes_yleft= (np.roll(source_data,-1,axis=0)!=0)
			source_mask_indexes= (source_indexes | source_indexes_xright | source_indexes_xleft | source_indexes_yright | source_indexes_yleft)
	
			#source_mask= np.where(source_data!=0,1,0)
			taken_pixels= np.where(sources_data[source_mask_indexes]!=0) # get list of taken pixels in main mask corresponding to this source
			has_taken_pixels= np.any(taken_pixels)
			if has_taken_pixels:
				print 'INFO: Source pixels have been already taken by a previous generated source, regenerate...'
				continue
				
			# Add to extended source data and mask
			sources_data+= source_data
			ngen_sources+= 1

			# Set model map
			ix= int(np.round(x0))
			iy= int(np.round(y0))
			
			# Make Caesar source	
			source_name= 'Sext' + str(ngen_sources)
			source_id= ngen_sources
			source_type= Caesar.Source.eExtended
			caesar_source= self.make_caesar_source(source_data,source_name,source_id,source_type,source_sim_type,None,None,None,source_max_scale)
			if caesar_source is None:
				print('Generate source has too few pixels, skip and regenerate...')
				continue

			self.caesar_sources.append(caesar_source)

			print ('INFO: Ext Source %s: Pos(%s,%s), ix=%s, iy=%s, S=%s' % (source_name,str(x0),str(y0),str(ix),str(iy),str(S)))


		return sources_data

	
	#####################################
	###         GENERATE MAP           ##
	#####################################
	def generate_map(self):
		""" Generate sky map """
		
		## == INITIALIZE DATA ==
		print ('INFO: Initializing simulator data...')
		self.init()

		## == GENERATE EMPTY IMAGE ==
		data = Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.nx, y_width=2*self.ny)(self.gridx, self.gridy)
		mask_data = Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.nx, y_width=2*self.ny)(self.gridx, self.gridy)

		## == GENERATE BKG ==
		if self.simulate_bkg:
			print ('INFO: Generating map bkg...')
			bkg_data= self.generate_bkg()
			data+= bkg_data
	
		## == GENERATE COMPACT SOURCES ==
		if self.simulate_compact_sources:
			print ('INFO: Generating compact sources...')
			[compact_source_data,compact_source_mask_data] = self.generate_compact_sources()
			data+= compact_source_data
			mask_data+= compact_source_mask_data
	
		## == GENERATE EXTENDED SOURCES ==
		if self.simulate_extended_sources:
			print ('INFO: Generating extended sources...')
			ext_source_data = self.generate_extended_sources()
			data+= ext_source_data
			mask_data+= ext_source_data

		## == MAKE FINAL MAP ==
		print ('INFO: Creating final map with bkg + sources added...')

		## Sum data in cumulative map
		#data= bkg_data + compact_source_data + ext_source_data
		#mask_data= compact_source_mask_data + ext_source_data

		## Cast data from float64 to float32 
		data_casted = data.astype(np.float32)
		mask_data_casted = mask_data.astype(np.float32)

		## Convert data from Jy/pixel to Jy/beam
		## Jy/pixel= Jy/beam / beamArea(pixels)
		scaleFactor= self.beam_area
		data_casted*= scaleFactor

		## Create Caesar skymodel image from data (units= Jy/pixel)
		print ('INFO: Creating Caesar image from data...')
		##self.caesar_img= self.make_caesar_image(data_casted)    # set toy sim map data
		self.caesar_img= self.make_caesar_image(mask_data_casted) # set skymodel map data

		## == WRITE MAPS TO FITS FILES ==
		print ('INFO: Writing images to FITS...')
		self.write_map(data_casted,self.mapfilename)
		self.write_source_map(mask_data_casted,self.modelfilename)

		## == WRITE IMG & SOURCES TO ROOT FILE ==
		print ('INFO: Writing image & source collection to ROOT file...')
		self.save()

		## == WRITE DS9 REGION FILE ==
		print ('INFO: Writing DS9 regions...')
		self.write_ds9_regions()

		#return [data_casted,mask_data_casted]


	def write_ds9_regions(self):
		""" Write DS9 regions with sim sources """
		
		## Open file	
		fout = open(self.ds9filename, 'wb')
	
		## Write file header
		fout.write('global color=white font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n')
		fout.write('image\n')

		## Write source contour region
		for item in self.caesar_sources:
			regionInfo= item.GetDS9Region(True)
			fout.write(regionInfo)
			fout.write('\n')			

		fout.close();

	def write_casa_mask(self,boxsize=10):
		""" Write CASA mask file around simulated sources"""

		## Create a WCS structure
		w = wcs.WCS(naxis=2)
		w.wcs.crpix = [self.crpix1, self.crpix2]
		w.wcs.cdelt = np.array([-self.pixsize/3600., self.pixsize/3600.])
		w.wcs.crval = [self.crval1, self.crval2]
		w.wcs.ctype = [self.ctype1, self.ctype2]
		#w.wcs.set_pv([(2, 1, 45.0)])

		## Create mask ascii file with header
		f = open(str(self.casafilename), 'wb')
		f.write("#CRTFv0\n")
		#f.write("global coord = J2000, color=blue\n")

		## Create a CASA box around the source
		for item in self.caesar_sources:
			ix_min= item.GetIxMin()
			ix_max= item.GetIxMax()
			iy_min= item.GetIyMin()
			iy_max= item.GetIyMax()

			# Set box coordinates
			pixcrd = np.array([[max(0,ix_min-boxsize/2.), max(0,iy_min-boxsize/2.)], [min(self.nx-1,ix_max+boxsize/2.), min(self.ny-1,iy_max+boxsize/2.)]], np.float_)
			
			# Convert pixel coordinates to world coordinates
			world = w.wcs_pix2world(pixcrd, 1)
			print(world)
			f.write("box [ [{0}deg,{1}deg], [{2}deg,{3}deg] ]\n".format(min(world[0,0],world[1,0]),min(world[0,1],world[1,1]),max(world[0,0],world[1,0]),max(world[0,1],world[1,1])))

		# Close ascii file
		f.close()

	def draw_map(self,data):
		""" Draw map data """
		plt.imshow(data, origin='lower', cmap="hot")
		pylab.show()

	def write_map(self,data,outputfile):
		""" Write FITS image with sim data """

		# Define FITS header
		header= fits.Header()	
		header.set('SIMPLE','T')
		header.set('BITPIX','-32')
		header.set('NAXIS1', str(self.nx))
		header.set('NAXIS2', str(self.ny))
		#header.set('NAXIS3', 1)
		#header.set('NAXIS4', 1)
		header.set('BUNIT', 'JY/BEAM')
		header.set('BMAJ', self.beam_bmaj/3600.)
		header.set('BMIN', self.beam_bmin/3600.)
		header.set('BPA', self.beam_bpa)
		header.set('BSCALE',1.)
		header.set('BZERO',0.)
		header.set('CDELT1',-self.pixsize/3600.)
		header.set('CDELT2',self.pixsize/3600.)
		header.set('CTYPE1',self.ctype1)
		header.set('CTYPE2',self.ctype2)
		header.set('CRPIX1',self.crpix1)
		header.set('CRPIX2',self.crpix2)
		header.set('CRVAL1',self.crval1)
		header.set('CRVAL2',self.crval2)

		# Define HDU
		hdu = fits.PrimaryHDU(data=data,header=header)
		hdulist = fits.HDUList([hdu])
		hdulist.writeto(outputfile,overwrite=True)

	def write_source_map(self,data,outputfile):
		""" Write FITS image with sim mask data """
		
		# Define FITS header
		header= fits.Header()	
		header.set('SIMPLE','T')
		header.set('BITPIX','-32')
		header.set('NAXIS1', str(self.nx))
		header.set('NAXIS2', str(self.ny))
		header.set('BUNIT', 'JY/pixel')
		header.set('BMAJ', self.beam_bmaj/3600.)
		header.set('BMIN', self.beam_bmin/3600.)
		header.set('BPA', self.beam_bpa)
		header.set('BSCALE',1.)
		header.set('BZERO',0.)
		header.set('CDELT1',-self.pixsize/3600.)
		header.set('CDELT2',self.pixsize/3600.)
		header.set('CTYPE1',self.ctype1)
		header.set('CTYPE2',self.ctype2)
		header.set('CRPIX1',self.crpix1)
		header.set('CRPIX2',self.crpix2)
		header.set('CRVAL1',self.crval1)
		header.set('CRVAL2',self.crval2)

		# Define HDU
		hdu = fits.PrimaryHDU(data=data,header=header)
		hdulist = fits.HDUList([hdu])
		hdulist.writeto(outputfile,overwrite=True)

	def save(self):
		""" Write img & source collection to ROOT file """
	
		# Loop over sources
		print ('Filling #%s sources to ROOT tree...' % str(len(self.caesar_sources)) )
		for item in self.caesar_sources:
			#self.cs= item
			item.Copy(self.cs)
			self.cs.Print()
			self.outtree.Fill()

		# Write to file
		self.outfile.cd()
		self.caesar_img.Write()
		self.outtree.Write()
		self.outfile.Close()

		# Write CASA mask file
		self.write_casa_mask(boxsize=self.mask_boxsize)		

	
###########################




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

	# - Image args
	Nx= args.nx
	Ny= args.ny
	marginX= args.marginx
	marginY= args.marginy
	pixsize= args.pixsize
	ctype1= args.ctype1
	ctype2= args.ctype2
	crpix1= args.crpix1
	crpix2= args.crpix2
	crval1= args.crval1
	crval2= args.crval2

	#- Source model
	model_trunc_zmin= args.model_trunc_zmin
	trunc_thr= args.trunc_thr
	npixels_min= args.npixels_min

	# - Bkg info args
	enable_bkg= args.enable_bkg
	bkg_level= args.bkg_level
	bkg_rms= args.bkg_rms

	# - Compact source args
	enable_compactsources= args.enable_compactsources 
	Bmaj= args.bmaj
	Bmin= args.bmin
	Bpa= args.bpa
	Zmin= args.zmin
	Zmax= args.zmax
	source_density= args.source_density

	# - Extended source args
	enable_extsources= args.enable_extsources
	ext_source_type= args.ext_source_type
	Zmin_ext= args.zmin_ext
	Zmax_ext= args.zmax_ext
	ext_source_density= args.ext_source_density
	ext_scale_min= args.ext_scale_min
	ext_scale_max= args.ext_scale_max
	ring_rmin= args.ring_rmin
	ring_rmax= args.ring_rmax
	ring_wmin= args.ring_wmin
	ring_wmax= args.ring_wmax
	ellipse_rmin= args.ellipse_rmin
	ellipse_rmax= args.ellipse_rmax
	disk_shell_ampl_ratio_min= args.disk_shell_ampl_ratio_min 
	disk_shell_ampl_ratio_max= args.disk_shell_ampl_ratio_max
	disk_shell_radius_ratio_min= args.disk_shell_radius_ratio_min 
	disk_shell_radius_ratio_max= args.disk_shell_radius_ratio_max

	# - Output args
	outputfile= args.outputfile
	mask_outputfile= args.outputfile_model
	outputfile_sources= args.outputfile_sources
	outputfile_ds9region= args.outputfile_ds9region
	outputfile_casaregion= args.outputfile_casaregion
	mask_boxsize= args.mask_boxsize

	print("*** ARGS ***")
	print("Nx: %s" % Nx)
	print("Ny: %s" % Ny)
	print("Margin X: %s" % marginX)
	print("Margin Y: %s" % marginY)
	print("pixsize: %s" % pixsize)
	print("ctype: (%s %s)" % (ctype1,ctype2))
	print("crpix: (%s %s)" % (crpix1,crpix2))
	print("crval: (%s %s)" % (crval1,crval2))
	print("Beam (Bmaj/Bmin/Bpa): (%s,%s,%s)" % (Bmaj, Bmin, Bpa))
	print("Enable bkg? %s" % str(enable_bkg) )
	print("Bkg info (level,rms): (%s,%s)" % (bkg_level, bkg_rms))
	print("Enable compact sources? %s" % str(enable_compactsources) )
	print("Source significance range: (%s,%s)" % (Zmin, Zmax))
	print("Source density (deg^-2): %s" % source_density)	
	print("Enable extended sources? %s" % str(enable_extsources) )
	print("Extended source type %s" %str(ext_source_type) )
	print("Extended source significance range: (%s,%s)" % (Zmin_ext, Zmax_ext))
	print("Extended source density (deg^-2): %s" % ext_source_density)
	print("Extended source scale min/max: (%s,%s)" % (ext_scale_min, ext_scale_max))
	print("Output filename: %s " % outputfile)
	print("Model trunc thr: %s " % str(trunc_thr))
	print("Mask output filename: %s " % mask_outputfile)
	print("Mask box size: %s " % mask_boxsize)
	print("************")

	## Generate simulated sky map
	print ('INFO: Generate simulated sky map...')
	simulator= SkyMapSimulator(Nx,Ny,pixsize)
	simulator.set_margins(marginX,marginY)
	simulator.set_ref_pix(crpix1,crpix2)
	simulator.set_ref_pix_coords(crval1,crval2)
	simulator.set_coord_system_type(ctype1,ctype2)

	simulator.set_model_trunc_thr(trunc_thr)
	simulator.set_model_trunc_significance(model_trunc_zmin)
	simulator.set_npixels_min(npixels_min)
	simulator.set_map_filename(outputfile)
	simulator.set_model_filename(mask_outputfile)
	simulator.set_source_filename(outputfile_sources)
	simulator.set_ds9region_filename(outputfile_ds9region)
	simulator.set_casaregion_filename(outputfile_casaregion)
	simulator.enable_bkg(enable_bkg)
	simulator.set_bkg_pars(bkg_level,bkg_rms)
	simulator.set_beam_info(Bmaj,Bmin,Bpa)	
	simulator.enable_compact_sources(enable_compactsources)
	simulator.set_source_significance_range(Zmin,Zmax)
	simulator.set_source_density(source_density)
	simulator.enable_extended_sources(enable_extsources)
	simulator.set_ext_source_type(ext_source_type)
	simulator.set_ext_source_significance_range(Zmin_ext,Zmax_ext)
	simulator.set_ext_source_density(ext_source_density)
	#simulator.set_ring_pars(ring_rmin,ring_rmax,ring_wmin,ring_wmax)
	simulator.set_ring_pars(ext_scale_min,ext_scale_max,ring_wmin,ring_wmax)
	#simulator.set_ellipse_pars(ellipse_rmin,ellipse_rmax)
	simulator.set_ellipse_pars(ext_scale_min,ext_scale_max)
	simulator.set_disk_pars(ext_scale_min,ext_scale_max)
	simulator.set_disk_shell_pars(disk_shell_ampl_ratio_min,disk_shell_ampl_ratio_max,disk_shell_radius_ratio_min,disk_shell_radius_ratio_max)
	simulator.set_mask_box_size(mask_boxsize)
	#[data, mask_data]= simulator.generate_map()
	simulator.generate_map()

	## Write fits
	#print ('INFO: Writing images to FITS...')
	#simulator.write_map(data,outputfile)
	#simulator.write_source_map(mask_data,mask_outputfile)

	## Write sources to ROOT
	#print ('INFO: Writing source collection to ROOT TTree...')
	#simulator.write_source_tree(outputfile_sources)

	## Draw image
	#print ('INFO: Draw image...')
	#plt.imshow(data, origin='lower', cmap="hot")
	#pylab.show()

###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	#main()
	sys.exit(main())

