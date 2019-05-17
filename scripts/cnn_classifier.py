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
import keras
from keras import layers
from keras import models
from keras import optimizers
from keras.utils import plot_model
from keras import backend
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

## ADDON ML MODULES
from sklearn.model_selection import train_test_split

## GRAPHICS MODULES
import matplotlib.pyplot as plt



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
	parser.add_argument('-normdatamin', '--normdatamin', dest='normdatamin', required=False, type=float, default=-0.0100, action='store',help='Normalization min used to scale data in (0,1) range (default=-100 mJy/beam)')	
	parser.add_argument('-normdatamax', '--normdatamax', dest='normdatamax', required=False, type=float, default=10, action='store',help='Normalization max used to scale data in (0,1) range (default=10 Jy/beam)')	
	parser.add_argument('-nx', '--nx', dest='nx', required=False, type=int, default=101, action='store',help='Image width in pixels (default=101)')
	parser.add_argument('-ny', '--ny', dest='ny', required=False, type=int, default=101, action='store',help='Image height in pixels (default=101)')	
	parser.add_argument('-nsamples_bkg', '--nsamples_bkg', dest='nsamples_bkg', required=False, type=int, default=10, action='store',help='Number of train images for bkg extracted from input maps (default=10)')
	parser.add_argument('-nsamples_source', '--nsamples_source', dest='nsamples_source', required=False, type=int, default=-1, action='store',help='Number of train images extracted around sources from input maps (default=-1)')	
	parser.add_argument('-nmaxobjects', '--nmaxobjects', dest='nmaxobjects', required=False, type=int, default=5, action='store',help='Max number of predicted objects in target (default=5)')
	parser.add_argument('-ntargetpars', '--ntargetpars', dest='ntargetpars', required=False, type=int, default=6, action='store',help='Nmber of pars per objects in target (default=6)')	
	parser.add_argument('-conv_nfilt_min', '--conv_nfilt_min', dest='conv_nfilt_min', required=False, type=int, default=16, action='store',help='Number of min convolution filters used (default=16)')	
	parser.add_argument('-conv_nfilt_max', '--conv_nfilt_max', dest='conv_nfilt_max', required=False, type=int, default=32, action='store',help='Number of max convolution filters used (default=32)')
	parser.add_argument('-dense_size_min', '--dense_size_min', dest='dense_size_min', required=False, type=int, default=16, action='store',help='Number of min neurons used in dense layer(default=16)')
	parser.add_argument('-dense_size_max', '--dense_size_max', dest='dense_size_max', required=False, type=int, default=32, action='store',help='Number of max neurons used in dense layer(default=32)')
	parser.add_argument('-test_size', '--test_size', dest='test_size', required=False, type=float, default=0.2, action='store',help='Fraction of input data used for testing the network (default=0.2)')
	parser.add_argument('-spars_loss_weight', '--spars_loss_weight', dest='spars_loss_weight', required=False, type=float, default=1, action='store',help='Loss weight to be given to source pars learning (default=1)')
	parser.add_argument('-labels_loss_weight', '--labels_loss_weight', dest='labels_loss_weight', required=False, type=float, default=1, action='store',help='Loss weight to be given to source labels learning (default=1)')
	parser.add_argument('-nepochs', '--nepochs', dest='nepochs', required=False, type=int, default=100, action='store',help='Number of epochs used in network training (default=100)')
	parser.add_argument('--saveimg', dest='saveimg', action='store_true')	
	parser.set_defaults(saveimg=False)

	parser.add_argument('-outfile_loss', '--outfile_loss', dest='outfile_loss', required=False, type=str, default='nn_loss.png', action='store',help='Name of NN loss plot file (default=nn_loss.png)')
	parser.add_argument('-outfile_accuracy', '--outfile_accuracy', dest='outfile_accuracy', required=False, type=str, default='nn_accuracy.png', action='store',help='Name of NN accuracy plot file (default=nn_accuracy.png)')
	parser.add_argument('-outfile_model', '--outfile_model', dest='outfile_model', required=False, type=str, default='nn_model.png', action='store',help='Name of NN model plot file (default=nn_model.png)')

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
	crop_data= data[ymin:ymax+1,xmin:xmax+1]
	
	#print ('DEBUG: (xmin,xmax)=(%d,%d), (ymin,ymax)=(%d,%d)' % (xmin,xmax,ymin,ymax))

	#- Replace NAN with zeros and inf with large numbers
	np.nan_to_num(crop_data,False)

	return crop_data

# - Set bkg train data
def make_bkg_train_data(imglist,nsamples,imgsizex,imgsizey,nmaxobjects,ntargetpars,normmin,normmax,writeimg=False):
	""" Prepare bkg train data """

	# - Read data in list
	#train_data= []
	input_data= []	
	#target_size= nmaxobjects*ntargetpars
	#target_data= []
	output_size= nmaxobjects*ntargetpars
	output_data= []  
	#target_label_size= nmaxobjects
	#target_label_data= []
	output_label_size= nmaxobjects
	output_label_data= []
	
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
			#print ('DEBUG: (x0,y0)=(%d,%d)' % (x0,y0))

			# - Extract crop img data
			data_crop= crop_img(data,x0,y0,imgsizex,imgsizey)
			imgcropsize= np.shape(data_crop)
			#print ('INFO: Sample image no. %s/%s from image %d has size (%d,%d)' % (index+1,nsamples,imgcounter,imgcropsize[1],imgcropsize[0]) )	
 
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
			#train_data.append(data_crop)
			input_data.append(data_crop)

			# - Set train target & labels
			#target_data.append( np.zeros((1,target_size)) )
			#target_label_data.append( np.zeros((1,target_label_size)) )
			output_data.append( np.zeros((1,output_size)) )
			output_label_data.append( np.zeros((1,output_label_size)) )

			index+= 1

	#- Convert list to array
	#x_train= np.array(train_data)
	#x_train= x_train.astype('float32')
	inputs= np.array(input_data)
	inputs= inputs.astype('float32')
	
	#y_train= np.array(target_data)
	#y_train= y_train.astype('float32')
	outputs= np.array(output_data)
	outputs= outputs.astype('float32')

	outputs_shape= outputs.shape
	N= outputs_shape[0]
	outputs= outputs.reshape((N,output_size))

	#y_train_labels= np.array(target_label_data)
	#y_train_labels= y_train_labels.astype('float32')
	outputs_labels= np.array(output_label_data)
	outputs_labels= outputs_labels.astype('float32')
	outputs_labels= outputs_labels.reshape((N,output_label_size))

	# - Normalize to [0,1]
	inputs= (inputs - normmin)/(normmax-normmin)

	return inputs,outputs,outputs_labels


# - Set source train data
def make_source_train_data(filelist,nsamples,imgsizex,imgsizey,nmaxobjects,ntargetpars,normmin,normmax,writeimg=False):
	""" Prepare source train data """

	# - Read data in list
	#train_data= []
	input_data= []
	
	#target_size= nmaxobjects*ntargetpars
	#target_data= []  
	output_size= nmaxobjects*ntargetpars
	output_data= []  
	
	#target_label_size= nmaxobjects
	#target_label_data= []
	output_label_size= nmaxobjects
	output_label_data= []
	
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
		#print ('DEBUG: nsources=%d, npars=%d' % (nsources,npars) )

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
			#print ('DEBUG: sname=%s, (x0,y0)=(%d,%d)' % (sname,x0,y0))
			
			# - Extract crop img data
			data_crop= crop_img(data,x0,y0,imgsizex,imgsizey)
			imgcropsize= np.shape(data_crop)
			#print ('DEBUG: Sample image no. %s/%s from image %d has size (%d,%d)' % (index+1,nsamples,imgcounter,imgcropsize[1],imgcropsize[0]) )	
 
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
			#train_data.append(data_crop)
			input_data.append(data_crop)

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
			#targets= np.zeros((1,target_size))
			#target_labels= np.zeros((1,target_label_size))
			targets= np.zeros((1,output_size))
			target_labels= np.zeros((1,output_label_size))
			
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

			#target_data.append(targets)
			#target_label_data.append(target_labels)
			output_data.append(targets)
			output_label_data.append(target_labels)
			
	#- Convert list to array
	#x_train= np.array(train_data)
	#x_train= x_train.astype('float32')
	inputs= np.array(input_data)
	inputs= inputs.astype('float32')

	#y_train= np.array(target_data)
	#y_train= y_train.astype('float32')
	outputs= np.array(output_data)
	outputs= outputs.astype('float32')

	outputs_shape= outputs.shape
	N= outputs_shape[0]
	outputs= outputs.reshape((N,output_size))

	#y_train_labels= np.array(target_label_data)
	#y_train_labels= y_train_labels.astype('float32')
	outputs_labels= np.array(output_label_data)
	outputs_labels= outputs_labels.astype('float32')
	outputs_labels= outputs_labels.reshape((N,output_label_size))

	# - Normalize to [0,1]
	#x_train= (x_train - normmin)/(normmax-normmin)
	inputs= (inputs - normmin)/(normmax-normmin)

	#return x_train,y_train,y_train_labels
	return inputs,outputs,outputs_labels


###########################
##     BUILD NETWORK
###########################
def build_network(img_height, img_width, nmaxobjects, ntargetpars, conv_nfilt_min=16, conv_nfilt_max=32, conv_kern_size_min=3, conv_kern_size_max=3,conv_act='relu', pool_size=2, dense_size_min=16,dense_size_max=32, dense_act='relu'):
	""" Building deep network """

	dropout= 0.25
	dropout_dense= 0.5
	padding= "same"
	#padding= "valid"

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
	x = layers.Conv2D(filters=conv_nfilt_min, kernel_size=(conv_kern_size_min,conv_kern_size_min), activation=conv_act, padding=padding)(inputs)
	x = layers.MaxPooling2D(pool_size=(pool_size, pool_size),strides=None,padding='valid')(x)
	x = layers.Dropout(dropout)(x)

	x = layers.Conv2D(filters=conv_nfilt_max, kernel_size=(conv_kern_size_min,conv_kern_size_min), activation=conv_act, padding=padding)(x)
	x = layers.MaxPooling2D(pool_size=(pool_size, pool_size),strides=None,padding='valid')(x)
	x = layers.Dropout(dropout)(x)

	x = layers.Conv2D(filters=conv_nfilt_max, kernel_size=(conv_kern_size_max,conv_kern_size_max), activation=conv_act, padding=padding)(x)
	x = layers.MaxPooling2D(pool_size=(pool_size, pool_size),strides=None,padding='valid')(x)
	x = layers.Dropout(dropout)(x)

	#- Fully connected layers
	x = layers.Flatten()(x)
	x = layers.Dense(dense_size_max, activation=dense_act)(x)
	x = layers.Dropout(dropout_dense)(x)
	#x = layers.Dense(dense_size_min, activation=dense_act)(x)
	#x = layers.Dropout(dropout_dense)(x)

	# - Output layers
	type_prediction = layers.Dense(nmaxobjects, activation='sigmoid', name='type')(x)
	pars_prediction = layers.Dense(nmaxobjects*ntargetpars, activation='linear', name='pars')(x)

	#- Create NN model
	model = Model(
			inputs=inputs,
			outputs=[type_prediction, pars_prediction],
			name="SourceNet"
	)

	return model

#####################################
##     DEFINE NETWORK ADDON METRICS
#####################################
#- These metrics were removed in Keras
#  See https://gist.github.com/pishangujeniya/ca8dd46a5d5cf0b0391b712c1a03b9b6

def precision(y_true, y_pred):	
	"""Precision metric.	
	Only computes a batch-wise average of precision. Computes the precision, a
	metric for multi-label classification of how many selected items are
	relevant.
	"""	
	true_positives = keras.backend.sum(keras.backend.round(keras.backend.clip(y_true * y_pred, 0, 1)))	
	predicted_positives = keras.backend.sum(keras.backend.round(keras.backend.clip(y_pred, 0, 1)))	
	precision = true_positives / (predicted_positives + keras.backend.epsilon())	
	return precision

def recall(y_true, y_pred):	
	"""Recall metric.	
	Only computes a batch-wise average of recall. Computes the recall, a metric
	for multi-label classification of how many relevant items are selected.	
	"""	
	true_positives = keras.backend.sum(keras.backend.round(keras.backend.clip(y_true * y_pred, 0, 1)))	
	possible_positives = keras.backend.sum(keras.backend.round(keras.backend.clip(y_true, 0, 1)))	
	recall = true_positives / (possible_positives + keras.backend.epsilon())	
	return recall

def f1_score(y_true, y_pred):
	"""Computes the F1 Score
	Only computes a batch-wise average of recall. Computes the recall, a metric
	for multi-label classification of how many relevant items are selected.	
	"""
	p = precision(y_true, y_pred)
	r = recall(y_true, y_pred)
	return (2 * p * r) / (p + r + keras.backend.epsilon())

	
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
	normdatamax= args.normdatamax
	conv_nfilt_min= args.conv_nfilt_min
	conv_nfilt_max= args.conv_nfilt_max
	dense_size_min= args.dense_size_min
	dense_size_max= args.dense_size_max
	test_size= args.test_size
	spars_loss_weight= args.spars_loss_weight
	labels_loss_weight= args.labels_loss_weight
	nepochs= args.nepochs
	outfile_loss= args.outfile_loss
	outfile_accuracy= args.outfile_accuracy
	outfile_model= args.outfile_model

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
	# - Extract train data for bkg from images
	print ('INFO: Generating bkg train data ...')
	#(x_train_bkg,y_train_bkg,y_train_labels_bkg)= make_bkg_train_data(
	(inputs_bkg,outputs_bkg,outputs_labels_bkg)= make_bkg_train_data(
		imglist=filenames_bkg,
		nsamples=nsamples_bkg,
		imgsizex=nx,imgsizey=ny,
		nmaxobjects=nmaxobjects,ntargetpars=ntargetpars,
		normmin=normdatamin,normmax=normdatamax,
		writeimg=saveimg
	)

	print 'DEBUG: inputs_bkg size=', np.shape(inputs_bkg)
	print 'DEBUG: outputs_bkg size=', np.shape(outputs_bkg)
	print 'DEBUG: outputs_labels_bkg size', np.shape(outputs_labels_bkg)
	print 'DEBUG: outputs_bkg=',outputs_bkg
	print 'DEBUG: outputs_labels_bkg=',outputs_labels_bkg

	# - Extract train data for sources from images
	print ('INFO: Generating source train data ...')
	#(x_train_source,y_train_source,y_train_labels_source)= make_source_train_data(
	(inputs_source,outputs_source,outputs_labels_source)= make_source_train_data(
		filelist=filenames_source,
		nsamples=nsamples_source,
		imgsizex=nx,imgsizey=ny,
		nmaxobjects=nmaxobjects,ntargetpars=ntargetpars,
		normmin=normdatamin,normmax=normdatamax,
		writeimg=saveimg
	)

	print 'DEBUG: inputs_source size=', np.shape(inputs_source)
	print 'DEBUG: outputs_source size=', np.shape(outputs_source)
	print 'DEBUG: outputs_labels_source size=', np.shape(outputs_labels_source)
	print 'DEBUG: outputs_source=',outputs_source
	print 'DEBUG: outputs_labels_source=',outputs_labels_source


	# - Merge data for bkg & sources
	print 'INFO: Merging train data for bkg & sources ...'
	#inputs= []
	#inputs.append(inputs_bkg)
	#inputs.append(inputs_source)
	inputs= np.concatenate((inputs_bkg,inputs_source))

	#outputs= []	
	#outputs.append(outputs_bkg)
	#outputs.append(outputs_source)
	outputs= np.concatenate((outputs_bkg,outputs_source))

	#outputs_labels= []	
	#outputs_labels.append(outputs_labels_bkg)
	#outputs_labels.append(outputs_labels_source)
	outputs_labels= np.concatenate((outputs_labels_bkg,outputs_labels_source))

	# - Shuffle data before splitting in test & validation sample
	print 'INFO: Shuffling train data ...'
	indices= np.arange(inputs.shape[0])
	np.random.shuffle(indices)
	inputs= inputs[indices]
	outputs= outputs[indices]
	outputs_labels= outputs_labels[indices]
	
	print 'DEBUG: inputs size=', np.shape(inputs)
	print 'DEBUG: outputs size=', np.shape(outputs)
	print 'DEBUG: outputs_labels size=', np.shape(outputs_labels)

	# - Partition the data into training and cross-validation splits
	print 'INFO: Splitting data into train & test samples ...'
	split= train_test_split(
		inputs,outputs,outputs_labels, 
		test_size=test_size, 
		random_state=None
	)
	(inputs_train, inputs_test, outputs_train, outputs_test, outputs_labels_train, outputs_labels_test) = split

	print 'DEBUG: inputs_train size=', np.shape(inputs_train)
	print 'DEBUG: inputs_test size=', np.shape(inputs_test)
	print 'DEBUG: outputs_train size=', np.shape(outputs_train)
	print 'DEBUG: outputs_test size=', np.shape(outputs_test)
	print 'DEBUG: outputs_labels_train size=', np.shape(outputs_labels_train)
	print 'DEBUG: outputs_labels_test size=', np.shape(outputs_labels_test)

	#===========================
	#==   BUILD NN
	#===========================
	#- Create the network
	print ('INFO: Building network architecture ...')
	model= build_network(
		img_height=nx, img_width=ny, 
		nmaxobjects=nmaxobjects, ntargetpars=ntargetpars,
		conv_nfilt_min=conv_nfilt_min,conv_nfilt_max=conv_nfilt_max,
		dense_size_min=dense_size_min,dense_size_max=dense_size_max
	)

	# - Print network architecture
	model.summary()

	#- Set optimizer & loss function per each output
	print ('INFO: Compiling network...')
	opt= optimizers.RMSprop(lr=1.e-4)
	#opt= Adam(lr=INIT_LR, decay=INIT_LR / nepochs)
	
	losses = {
		"type": "binary_crossentropy",
		"pars": "mse"
	}
	lossWeights = {
		"type": labels_loss_weight,
		"pars": spars_loss_weight
	}

	model.compile(optimizer=opt,loss=losses, loss_weights=lossWeights, metrics=['accuracy',precision,recall,f1_score])

	#===========================
	#==   TRAIN NN
	#===========================
	print ('INFO: Training network...')
	fitout= model.fit(
		x=inputs_train, 
		y={"type": outputs_labels_train,"pars": outputs_train},
		validation_data=(inputs_test,{"type": outputs_labels_test,"pars": outputs_test}),
		#batch_size=64
		epochs=nepochs,
		verbose=1
	)
	
	#===========================
	#==   PLOT NN RESULTS
	#===========================
	# - Plot the network	
	plot_model(model, to_file=outfile_model)

	# - Plot the total loss, type loss, spars loss
	lossNames = ["loss", "type_loss", "pars_loss"]
	plt.style.use("ggplot")
	(fig, ax) = plt.subplots(3, 1, figsize=(13, 13))

	for (i, l) in enumerate(lossNames):
		# Plot the loss for both the training and validation data
		title = "Loss for {}".format(l) if l != "loss" else "Total loss"
		ax[i].set_title(title)
		ax[i].set_xlabel("Epoch #")
		ax[i].set_ylabel("Loss")
		ax[i].plot(np.arange(0, nepochs), fitout.history[l], label="TRAIN SAMPLE - " + l)
		ax[i].plot(np.arange(0, nepochs), fitout.history["val_" + l], label="TEST SAMPLE - " + l)
		ax[i].legend()

	plt.tight_layout()
	plt.savefig(outfile_loss)
	plt.close()

	# - Plot the accuracy
	accuracyNames = ["type_acc", "pars_acc"]
	plt.style.use("ggplot")
	(fig, ax) = plt.subplots(2, 1, figsize=(8, 8))

	for (i, l) in enumerate(accuracyNames):
		# Plot the loss for both the training and validation data
		ax[i].set_title("Accuracy for {}".format(l))
		ax[i].set_xlabel("Epoch #")
		ax[i].set_ylabel("Accuracy")
		ax[i].plot(np.arange(0, nepochs), fitout.history[l], label="TRAIN SAMPLE - " + l)
		ax[i].plot(np.arange(0, nepochs), fitout.history["val_" + l], label="TEST SAMPLE - " + l)
		ax[i].legend()

	plt.tight_layout()
	plt.savefig(outfile_accuracy)
	plt.close()

###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

