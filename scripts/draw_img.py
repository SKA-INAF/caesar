
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

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## ASTRO
from astropy.io import fits
from astropy.wcs import WCS
from astropy.units import Quantity
from astropy.modeling.parameters import Parameter
from astropy.modeling.core import Fittable2DModel
from astropy import wcs
from astropy.wcs.utils import wcs_to_celestial_frame
from astropy import units as u
from astropy.visualization import ZScaleInterval, SqrtStretch, LinearStretch, LogStretch, AsinhStretch, ImageNormalize, MinMaxInterval
from astropy.visualization import ContrastBiasStretch
from astropy.visualization import make_lupton_rgb
#import aplpy

import regions
from regions import DS9Parser
from regions import read_ds9

from shapely.geometry import Polygon

## Graphics modules
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

## LOGGER
import logging
import logging.config
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s - %(message)s",datefmt='%Y-%m-%d %H:%M:%S')
logger= logging.getLogger(__name__)
logger.setLevel(logging.INFO)



###########################
##     ARGS
###########################
def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")

	parser.add_argument('-img','--img', dest='img', required=True, type=str, default='', help='Input fits file') 
	parser.add_argument('-region','--region', dest='region', required=False, type=str, default='', help='Region filename to be superimposed (optional)') 
	parser.add_argument('-zmin','--zmin', dest='zmin', required=False, type=float, default=0, help='Min zscale value (default=0)') 
	parser.add_argument('-zmax','--zmax', dest='zmax', required=False, type=float, default=0, help='Max zscale value (default=0)') 
	parser.add_argument('-cmap','--cmap', dest='cmap', required=False, type=str, default='afmhot', help='Color palette') 
	parser.add_argument('-contrast','--contrast', dest='contrast', required=False, type=float, default=0.3, help='zscale contrast value (default=0.3)')
	parser.add_argument('--showgrid', dest='showgrid', action='store_true')	
	parser.set_defaults(showgrid=False)
	parser.add_argument('--showlabels', dest='showlabels', action='store_true')	
	parser.set_defaults(showlabels=False)
	parser.add_argument('--wcs', dest='wcs', action='store_true')	
	parser.set_defaults(wcs=False)
	parser.add_argument('--save', dest='save', action='store_true')	
	parser.set_defaults(save=False)
	parser.add_argument('-outfile','--outfile', dest='outfile', required=False, type=str, default='plot.png', help='Output filename where to save the plot (if enabled)') 
	
	args = parser.parse_args()	

	return args



##############
##   MAIN   ##
##############
def main():
	"""Main function"""
	
	#===========================
	#==   PARSE ARGS
	#===========================
	logger.info("Parse script args ...")
	try:
		args= get_args()
	except Exception as ex:
		logger.error("Failed to get and parse options (err=%s)",str(ex))
		return 1

	# - Input files
	imgfile= args.img
	regionfile= args.region
	zmin= args.zmin
	zmax= args.zmax
	cmap= args.cmap
	contrast= args.contrast
	plot_wcs= args.wcs
	draw_grid= args.showgrid
	draw_labels= args.showlabels
	save= args.save
	outfile= args.outfile

	#===========================
	#==   READ REGION
	#===========================
	if regionfile!="":
		logger.info("Reading region file %s ..." % regionfile)
		regs= regions.read_ds9(regionfile)

	#===========================
	#==   READ IMAGE
	#===========================
	# - Read fits
	hdu= fits.open(imgfile)
	data= hdu[0].data
	header= hdu[0].header
	#header['BMAJ']= bmaj/3600
	#header['BMIN']= bmin/3600
	#header['BPA']= pa

	# - Remove 3 and 4 channels
	nchan= len(data.shape)
	logger.info("Input image has %d channels..." % nchan)

	if nchan==4:
		logger.info("Removing fit degenerated axis 3 & 4 ...")
		data= data[0,0,:,:]
		header['NAXIS'] = 2
		if 'NAXIS3' in header:
			del header['NAXIS3']
		if 'NAXIS4' in header:
			del header['NAXIS4']
		if 'CTYPE3' in header:
			del header['CTYPE3']
		if 'CRVAL3' in header:
			del header['CRVAL3']
		if 'CDELT3' in header:
			del header['CDELT3']
		if 'CRPIX3' in header:
			del header['CRPIX3']
		if 'CROTA3' in header:
			del header['CROTA3']
		if 'CTYPE4' in header:
			del header['CTYPE4']
		if 'CRVAL4' in header:
			del header['CRVAL4']
		if 'CDELT4' in header:
			del header['CDELT4']
		if 'CRPIX4' in header:
			del header['CRPIX4']
		if 'CROTA4' in header:
			del header['CROTA4']
		if 'CUNIT3' in header:
			del header['CUNIT3']
		if 'CUNIT4' in header:
			del header['CUNIT4']

	elif nchan==3:
		logger.info("Removing fit degenerated axis 3 ...")
		data= data[0,:,:]
		header['NAXIS'] = 2
		if 'NAXIS3' in header:
			del header['NAXIS3']
		if 'CTYPE3' in header:
			del header['CTYPE3']
		if 'CRVAL3' in header:
			del header['CRVAL3']
		if 'CDELT3' in header:
			del header['CDELT3']
		if 'CRPIX3' in header:
			del header['CRPIX3']
		if 'CROTA3' in header:
			del header['CROTA3']
		if 'CUNIT3' in header:
			del header['CUNIT3']

	elif nchan==2:	

		if 'NAXIS3' in header:
			del header['NAXIS3']
		if 'CTYPE3' in header:		
			del header['CTYPE3']
		if 'CRVAL3' in header:			
			del header['CRVAL3']
		if 'CDELT3' in header:	
			del header['CDELT3']
		if 'CRPIX3' in header:	
			del header['CRPIX3']
		if 'CROTA3' in header:
			del header['CROTA3']
		if 'CUNIT3' in header:	
			del header['CUNIT3']

	hdu[0].header= header
	#hdu[0].header['BMAJ']= bmaj/3600
	#hdu[0].header['BMIN']= bmin/3600
	#hdu[0].header['BPA']= pa

	bunit= 'z'
	if 'BUNIT' in header:
		bunit= header['BUNIT']

	# - Get WCS info
	wcs = WCS(header, naxis=2)
	cs_name= "image"
	try:
		cs= wcs_to_celestial_frame(wcs)
		cs_name= cs.name
	except Exception as e:
		logger.warn("Failed to get celestial frame from wcs (err=%s), disable plot_wcs option ..." % str(e))
		plot_wcs= False

	# - Convert to mJy
	#data*= 1000
	#zmin*= 1.e+3
	#zmax*= 1.e+3
	#logger.info("data min/max=%f/%f" % (np.min(hdu[0].data),np.max(hdu[0].data)))
	
	# - Set stretch	
	stretch= LinearStretch()
	norm = ImageNormalize(data, interval=ZScaleInterval(contrast=contrast), stretch=stretch)
	
	# - Plot image	
	fig = plt.figure(figsize=(10,10))
	if plot_wcs:
		ax = fig.add_subplot(1, 1, 1, projection=wcs, label='overlays')
	else:
		ax = fig.add_subplot(1, 1, 1)
	if zmin<zmax:
		im= ax.imshow(data, origin='lower', vmin=zmin, vmax=zmax, cmap=cmap, norm=norm)
	else:
		im= ax.imshow(data, origin='lower', cmap=cmap, norm=norm)
	
	# - Set axis titles
	if plot_wcs:
		if cs_name=='galactic':
			ax.set_xlabel('Galactic Longitude (deg)',size=18, labelpad=0.7)
			ax.set_ylabel('Galactic Latitude (deg)',size=18)
		else:
			ax.set_xlabel('Right Ascension (deg)',size=18, labelpad=0.7)
			ax.set_ylabel('Declination (deg)',size=18)
	else:
		ax.set_xlabel('x',size=18, labelpad=0.7)
		ax.set_ylabel('y',size=18)

	# - Set axis ticks
	if plot_wcs:
		ax.coords[0].set_ticks(number=8)
		ax.coords[1].set_ticks(number=12)
		ax.coords[0].set_ticklabel(size=14)
		ax.coords[1].set_ticklabel(size=14)
		ax.coords[0].set_major_formatter('d.dd')
		ax.coords[1].set_major_formatter('d.dd')
		ax.coords[0].display_minor_ticks(True)
		ax.coords[1].display_minor_ticks(True)
	

	# - Draw grid
	if draw_grid and plot_wcs:
		logger.info("Drawing coordinate grid ...")
		overlay = ax.get_coords_overlay('galactic')
		overlay.grid(color='black', ls='dotted', alpha=0.8)
		overlay[0].set_axislabel('Galactic Longitude (deg)',size=14)
		overlay[1].set_axislabel('Galactic Latitude (deg)',size=14)
		overlay[0].set_major_formatter('d.ddd')
		overlay[1].set_major_formatter('d.ddd')
		overlay[0].set_ticks(color='black',number=10)
		overlay[1].set_ticks(color='black',number=10)
		overlay[0].set_ticklabel(size=12)
		overlay[1].set_ticklabel(size=12)

	# - Draw color bar
	color_bar_label= 'Brightness (' + bunit + ')'
	cb= fig.colorbar(im, orientation="vertical", pad=0.01, fraction=0.047)
	cb.set_label(color_bar_label, y=1.0, ha='right',size=12)
	cb.ax.tick_params(labelsize=12) 

	# - Set ticks
	plt.tick_params(axis='x', labelsize=14)
	plt.tick_params(axis='y', labelsize=14)


	# - Superimpose regions
	region_color= "lime"
	region_style= "dashed"
	#region_style= "solid"

	if regionfile!="":
		logger.info("Superimposing region ...")
		for r in regs:
			
			label= ''
			if draw_labels and 'text' in r.meta:
				label= r.meta['text']

			# - Check if region is in pixel coordinates
			is_pix_region= False
			try:
				r.to_pixel(wcs)
				is_pix_region= False
			except:
				is_pix_region= True

			if is_pix_region:	
				r_pix= r
				#r_sky= r.to_sky(wcs)
			else:
				r_pix= r.to_pixel(wcs)
				#r_sky= r		
			

			#if plot_wcs:
			#	if cs_name=='galactic':
			#		points_x= r_sky.vertices.galactic.l.value
			#		points_y= r_sky.vertices.galactic.b.value
			#	else:
			#		points_x= r_sky.vertices.ra.value
			#		points_y= r_sky.vertices.dec.value
			#	vertices= [tuple(x) for x in zip(points_x,points_y)]
			#	polygon = Polygon(vertices)
			#	x,y = polygon.exterior.xy
			#	plt.plot(x, y, color=region_color, linestyle=region_style, label=label)
			#else:
			#	r_pix.plot(ax=ax, color=region_color, linestyle=region_style, label=label)
			
			# - Even if draw_wcs is enabled, plot the region in pix coordinates and it will work
			r_pix.plot(ax=ax, color=region_color, linestyle=region_style, label=label)
			
			# - Draw labels
			if draw_labels:
				bbox= r.bounding_box
				label_color= 'black'
				if bbox is not None:
					xmin= bbox.ixmin
					xmax= bbox.ixmax
					ymin= bbox.iymin
					ymax= bbox.iymax
					x_label= xmin + 0.5*(xmax-xmin)
					#x_label= xmax
					y_label= ymax
					#y_label= ymin + 0.5*(ymax-ymin)
					plt.text(x_label, y_label, label, color=label_color)


	# - Save or display
	if save:
		plt.savefig(outfile, bbox_inches='tight')
	else:
		plt.show()

		

###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())
