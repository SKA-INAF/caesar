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

## KERAS MODULES
from keras import layers
from keras import models
from keras.models import Model
from keras.layers.normalization import BatchNormalization
from keras.layers.convolutional import Conv2D
from keras.layers.convolutional import MaxPooling2D
from keras.layers.core import Activation
from keras.layers.core import Dropout
from keras.layers.core import Lambda
from keras.layers.core import Dense
from keras.layers import Flatten
from keras.layers import Input
import tensorflow as tf


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
	parser.add_argument('-normdatamin', '--normdatamin', dest='normdatamin', required=False, type=float, default=0, action='store',help='Normalize data to 0->inf by subtracting this value (default=0)')	
	parser.add_argument('-nx', '--nx', dest='nx', required=False, type=int, default=101, action='store',help='Image width in pixels (default=101)')
	parser.add_argument('-ny', '--ny', dest='ny', required=False, type=int, default=101, action='store',help='Image height in pixels (default=101)')	
	parser.add_argument('-nsamples_bkg', '--nsamples_bkg', dest='nsamples_bkg', required=False, type=int, default=10, action='store',help='Number of train images for bkg extracted from input maps (default=10)')
	parser.add_argument('-nsamples_source', '--nsamples_source', dest='nsamples_source', required=False, type=int, default=-1, action='store',help='Number of train images extracted around sources from input maps (default=-1)')	
	parser.add_argument('-nmaxobjects', '--nmaxobjects', dest='nmaxobjects', required=False, type=int, default=5, action='store',help='Max number of predicted objects in target (default=5)')
	parser.add_argument('-ntargetpars', '--ntargetpars', dest='ntargetpars', required=False, type=int, default=6, action='store',help='Nmber of pars per objects in target (default=6)')	
	parser.add_argument('--saveimg', dest='saveimg', action='store_true')	
	parser.set_defaults(saveimg=False)

	args = parser.parse_args()	

	return args

###########################
##     READ INPUT DATA
###########################
# - Has pattern in string
def has_patterns(s,patterns):
	""" Return true if patterns are found in string """
	if not patterns:		
		return False

	found= False
	for pattern in patterns:
		found= pattern in s
		if found:
			break

	return found

# - Read filelist
def read_list(filename,skip_patterns=[]):
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
		if not line_fields:
			continue

		# Skip pattern
		skipline= has_patterns(line_fields[0],skip_patterns)
		if skipline:
			continue 		

		fields.append(line_fields)

	f.close()	

	return fields




# Read bkg data
def read_bkg_data(filename):
	""" Read input bkg data """

	# - Read list with files
	filelist_data= []
	try:
		filelist_data= read_list(filename,['#'])
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

	return filelist_data
	

# Read source data
def read_source_data(filename):
	""" Read input source data """
	filelist_data= []
	try:
		filelist_data= read_list(filename,['#'])
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

	return filelist_data

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
	np.nan_to_num(crop_data,False)

	return crop_data

# - Set bkg train data
def make_bkg_train_data(imglist,nsamples,imgsizex,imgsizey,nmaxobjects,ntargetpars,writeimg=False):
	""" Prepare bkg train data """

	# - Read data in list
	train_data= []
	target_size= nmaxobjects*ntargetpars
	target_data= []  
	target_label_size= nmaxobjects
	target_label_data= []
	Nchan= 1
	imgcounter= 0

	for item in imglist:
		imgcounter+= 1
		filename= item[0]
		print ('INFO: Reading file %s ...' % filename) 

		# - Read main bkg img
		try:
			data= read_img(filename)
		except Exception as ex:
			errmsg= 'Failed to read bkg image data (err=' + str(ex) + ')'
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

			# - Set train data as a tensor of size [Nsamples,Nx,Ny,Nchan] Nchan=1
			data_crop= data_crop.reshape(imgcropsize[0],imgcropsize[1],Nchan)
			train_data.append(data_crop)

			# - Set train target & labels
			target_data.append( np.zeros((1,target_size)) )
			target_label_data.append( np.zeros((1,target_label_size)) )

			index+= 1

	#- Convert list to array
	x_train= np.array(train_data)
	x_train= x_train.astype('float32')

	y_train= np.array(target_data)
	y_train= y_train.astype('float32')

	y_train_labels= np.array(target_label_data)
	y_train_labels= y_train_labels.astype('float32')

	return x_train,y_train,y_train_labels


# - Set source train data
def make_source_train_data(filelist,nsamples,imgsizex,imgsizey,nmaxobjects,ntargetpars,writeimg=False):
	""" Prepare source train data """

	# - Read data in list
	train_data= []
	target_size= nmaxobjects*ntargetpars
	target_data= []  
	target_label_size= nmaxobjects
	target_label_data= []
	Nchan= 1
	imgcounter= 0

	for item in filelist:
		imgcounter+= 1
		filename= item[0]
		filename_spar= item[1]
		print ('INFO: Reading files: %s, %s ...' % (filename,filename_spar) ) 

		# - Read main source img
		try:
			data= read_img(filename)
		except Exception as ex:
			errmsg= 'Failed to read source image data (err=' + str(ex) + ')'
			print "ERROR: " + errmsg
			raise IOError(errmsg)
	
		imgsize= np.shape(data)
		nx= imgsize[1]
		ny= imgsize[0]
		marginx= imgsizex/2
		marginy= imgsizey/2
		print ('INFO: Source image no. %d has size (%d,%d)' % (imgcounter,nx,ny) )

		# - Read source pars
		source_pars= []
		skip_patterns= ['#']
		try:
			source_pars= read_list(filename_spar,skip_patterns)
		except IOError:
			errmsg= 'Cannot read file: ' + filename_spar
			print "ERROR: " + errmsg
			raise IOError(errmsg)

		source_pars_size= np.shape(source_pars)
		nsources= source_pars_size[0]
		npars= source_pars_size[1]
		print ('DEBUG: nsources=%d, npars=%d' % (nsources,npars) )

		# - Extract sources img
		index= 0
		nsources_gen= nsources
		if nsamples!=-1 and nsamples<=nsources:
			nsources_gen= nsamples
	
		for index in range(nsources_gen):
			if index%100==0 :
				print ("INFO: Generating source image no. %s/%s from image %d ..." % (index+1,nsamples,imgcounter))
	
			# - Get source position & name
			sname= source_pars[index][0]
			source_x0= float(source_pars[index][1])
			source_y0= float(source_pars[index][2])
			x0= int(source_x0)
			y0= int(source_y0)
			xmin= int(x0-imgsizex/2)
			xmax= int(x0+imgsizex/2)
			ymin= int(y0-imgsizey/2)
			ymax= int(y0+imgsizey/2)		
			print ('DEBUG: sname=%s, (x0,y0)=(%d,%d)' % (sname,x0,y0))
			
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
			outfilename= 'train_source_' + str(index+1) + '-RUN' + str(imgcounter) + '.fits'
			if writeimg:
				write_fits(data_crop,outfilename)

			# - Set train data as a tensor of size [Nsamples,Nx,Ny,Nchan] Nchan=1
			data_crop= data_crop.reshape(imgcropsize[0],imgcropsize[1],Nchan)
			train_data.append(data_crop)

			# - Find all sources in range and sort sources in field from brightest to fainter			
			sources_in_field= []
			for sources in source_pars:	
				xx= float(sources[1])
				yy= float(sources[2])
				x= int(xx)
				y= int(yy)
				if x>=xmin and x<=xmax and y>=ymin and y<=ymax:
					sources_in_field.append(sources)

			sources_in_field.sort(key=lambda S: S[3],reverse=True)
			nsources_in_field= len(sources_in_field)
			
			# - Set train targets
			targets= np.zeros((1,target_size))
			target_labels= np.zeros((1,target_label_size))
			
			nobjs= min(nsources_in_field,nmaxobjects)
			par_counter= 0
			for k in range(nobjs):
				target_labels[0,k]= 1
				x0= float(sources_in_field[k][1])
				y0= float(sources_in_field[k][2])
				S= float(sources_in_field[k][3])
				sigmaX= float(sources_in_field[k][4])
				sigmaY= float(sources_in_field[k][5])
				theta= float(sources_in_field[k][6])	
				targets[0,par_counter+0]= x0 - xmin
				targets[0,par_counter+1]= y0 - ymin
				targets[0,par_counter+2]= S
				targets[0,par_counter+3]= sigmaX
				targets[0,par_counter+4]= sigmaY
				targets[0,par_counter+5]= theta
				par_counter+= 6

			target_data.append(targets)
			target_label_data.append(target_labels)
			
	#- Convert list to array
	x_train= np.array(train_data)
	x_train= x_train.astype('float32')

	y_train= np.array(target_data)
	y_train= y_train.astype('float32')

	y_train_labels= np.array(target_label_data)
	y_train_labels= y_train_labels.astype('float32')

	return x_train,y_train,y_train_labels


###########################
##     BUILD NETWORK
###########################
def build_network(img_height, img_width, nmaxobjects, ntargetpars, conv_nfilt_min=16, conv_nfilt_max=32, conv_kern_size=3, conv_act='relu', pool_size=2, dense_size_min=16,dense_size_max=32, dense_act='relu'):
	""" Building deep network """

	dropout= 0.25
	#model = models.Sequential()
	#model.add(layers.Conv2D(conv_nfilt_min, (conv_kern_size,conv_kern_size), activation=conv_tf, input_shape=(img_height,img_width, 1)))
	#model.add(layers.MaxPooling2D((pool_size, pool_size)))
	#model.add(layers.Conv2D(conv_nfilt_max, (conv_kern_size,conv_kern_size), activation=conv_tf))
	#model.add(layers.MaxPooling2D((pool_size, pool_size)))
	#model.add(layers.Conv2D(conv_nfilt_max, (conv_kern_size,conv_kern_size), activation=conv_tf))

	#- Input layer
	nchan= 1	
	inputShape = (img_height, img_width, nchan)
	inputs= Input(shape=inputShape,dtype='float', name='input')

	#- Convolutional layers
	x = layers.Conv2D(filters=conv_nfilt_min, kernel_size=(conv_kern_size,conv_kern_size), activation=conv_act, padding="same")(inputs)
	x = layers.MaxPooling2D(pool_size=(pool_size, pool_size),strides=None,padding='valid')(x)
	#x = layers.Dropout(dropout)(x)

	x = layers.Conv2D(filters=conv_nfilt_max, kernel_size=(conv_kern_size,conv_kern_size), activation=conv_act, padding="same")(x)
	x = layers.MaxPooling2D(pool_size=(pool_size, pool_size),strides=None,padding='valid')(x)
	#x = layers.Dropout(dropout)(x)

	x = layers.Conv2D(filters=conv_nfilt_max, kernel_size=(conv_kern_size,conv_kern_size), activation=conv_act, padding="same")(x)
	x = layers.MaxPooling2D(pool_size=(pool_size, pool_size),strides=None,padding='valid')(x)
	#x = layers.Dropout(dropout)(x)

	#- Fully connected layers
	x = layers.Flatten()(x)
	x = layers.Dense(dense_size_max, activation=dense_act)(x)
	x = layers.Dense(dense_size_min, activation=dense_act)(x)

	# - Output layers
	type_prediction = layers.Dense(nmaxobjects, activation='sigmoid', name='type')(x)
	pars_prediction = layers.Dense(nmaxobjects*ntargetpars, activation='linear', name='pars')(x)

	#- Create NN model
	model = Model(
			inputs=inputs,
			outputs=[type_prediction, pars_prediction],
			name="SourceNet"
	)

	#- Set loss function
	#model.compile(optimizer='rmsprop',loss=['mse', 'categorical_crossentropy', 'binary_crossentropy'])
	#model.compile(optimizer='rmsprop',loss={'age': 'mse','income': 'categorical_crossentropy','gender': 'binary_crossentropy'})

	return model
	
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
	nsamples_bkg= args.nsamples_bkg
	nsamples_source= args.nsamples_source
	nmaxobjects= args.nmaxobjects
	ntargetpars= args.ntargetpars
	saveimg= args.saveimg
	normdatamin= args.normdatamin
	
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
		filenames_source= read_source_data(filelist_source)
	except Exception as ex:
		print("Failed to read source data (err=%s)",str(ex))
		return 1

	#===========================
	#==   GENERATE TRAIN DATA
	#===========================
	print ('INFO: Generating bkg train data ...')
	(x_train_bkg,y_train_bkg,y_train_labels_bkg)= make_bkg_train_data(
		imglist=filenames_bkg,
		nsamples=nsamples_bkg,
		imgsizex=nx,imgsizey=ny,
		nmaxobjects=nmaxobjects,ntargetpars=ntargetpars,
		writeimg=saveimg
	)

	print 'INFO: x_train_bkg size', np.shape(x_train_bkg)
	print 'INFO: y_train_bkg size', np.shape(y_train_bkg)
	print 'INFO: y_train_labels_bkg size', np.shape(y_train_labels_bkg)
	print 'INFO: y_train_bkg',y_train_bkg
	print 'INFO: y_train_labels_bkg',y_train_labels_bkg

	print ('INFO: Generating source train data ...')
	(x_train_source,y_train_source,y_train_labels_source)= make_source_train_data(
		filelist=filenames_source,
		nsamples=nsamples_source,
		imgsizex=nx,imgsizey=ny,
		nmaxobjects=nmaxobjects,ntargetpars=ntargetpars,
		writeimg=saveimg
	)

	print 'INFO: x_train_source size', np.shape(x_train_source)
	print 'INFO: y_train_source size', np.shape(y_train_source)
	print 'INFO: y_train_labels_source size', np.shape(y_train_labels_source)
	print 'INFO: y_train_source',y_train_source
	print 'INFO: y_train_labels_source',y_train_labels_source

	#===========================
	#==   BUILD NN
	#===========================
	print ('INFO: Building network architecture ...')
	model= build_network(
		img_height=nx, img_width=ny, 
		nmaxobjects=nmaxobjects, ntargetpars=ntargetpars
	)

	model.summary()

	

###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

