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
from astropy.modeling.models import Box2D, Gaussian2D


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
from keras.models import load_model
from keras import backend
from keras.models import Model
from keras.preprocessing import image
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

#from keract import get_activations

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

	parser.add_argument('-inputimg', '--inputimg', dest='inputimg', required=True, type=str,action='store',help='Input image for which we want to draw nn response')
	parser.add_argument('-model', '--model', dest='model', required=True, type=str,action='store',help='NN model file to be loaded')
	parser.add_argument('-normdatamin', '--normdatamin', dest='normdatamin', required=False, type=float, default=-0.0100, action='store',help='Normalization min used to scale data in (0,1) range (default=-100 mJy/beam)')	
	parser.add_argument('-normdatamax', '--normdatamax', dest='normdatamax', required=False, type=float, default=10, action='store',help='Normalization max used to scale data in (0,1) range (default=10 Jy/beam)')	
	parser.add_argument('--normalize_inputs', dest='normalize_inputs', action='store_true')	
	parser.set_defaults(normalize_inputs=False)	

	args = parser.parse_args()	

	return args

####################
##   READ IMAGE   ##
####################
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


def get_activations(model,inputs):
	""" """

	# - Create activation model	
	layer_name= None
	#layer_outputs = [layer.output for layer in model.layers]
	layer_outputs = [layer.output for layer in model.layers if layer.name == layer_name or layer_name is None][1:]

	activation_model= models.Model(inputs=model.input, outputs=layer_outputs)

	# - Get layer names
	layer_names = []
	for layer in model.layers:
		layer_names.append(layer.name)

	print("== LAYER NAMES ==")
	print layer_names

	# - Get activations
	activations= activation_model.predict(inputs)
	images_per_row = 16

	#for layer_name, layer_activation in zip(layer_names, activations):
	#	n_features = layer_activation.shape[-1]
	#	size = layer_activation.shape[1]
	#	print("INFO: layer_name=%s, n_features=%d, size=%d" % (layer_name,n_features,size))
	#	print layer_activation.shape

	for layer_name, layer_activation in zip(layer_names, activations):

		# Skip input layer
		if layer_name=='input':
			continue

		# Draw only conv layers
		if layer_activation.ndim!=4:
			continue

		n_features = layer_activation.shape[-1]
		size = layer_activation.shape[1]
		n_cols = n_features # images_per_row

		print layer_activation.shape
		print("INFO: layer_name=%s, n_features=%d, size=%d, n_cols=%d" % (layer_name,n_features,size,n_cols))
		

		# - Fill channel image
		channel_images= []
		for index in range(n_features):
			channel_image = layer_activation[0,:,:, index]
			#channel_image -= channel_image.mean()
			#channel_image /= channel_image.std()
			#channel_image *= 64
			#channel_image += 128
			#channel_image = np.clip(channel_image, 0, 255).astype('uint8')
			channel_images.append(channel_image)

		# - Organize images in a grid
		nimg_row= int(np.sqrt(n_features))
		nimg_col= int(np.ceil(n_features/float(nimg_row)))
		print("INFO: Organize %d feature maps in pretty image of size (%d,%d)" % (n_features,nimg_row,nimg_col))

		display_grid = np.zeros((size*nimg_col, size*nimg_row))
		
		for col in range(nimg_col):
			for row in range(nimg_row):
				index= col * nimg_row + row
				if index>=n_features:
					continue
				print("INFO: Filling image at index=%d (col=%d, row=%d)" % (index,col,row))
				channel_image= channel_images[index]
				display_grid[col * size : (col + 1) * size, row * size : (row + 1) * size] = channel_image
			
		#display_grid = np.zeros((size * n_cols, images_per_row * size))

		#for col in range(n_cols):
		#	for row in range(images_per_row):
		#		index= col * images_per_row + row
		#		print("INFO: Filling image at index=%d (col=%d, row=%d)" % (index,col,row))
				
		#		if index>=n_features:
		#			continue

		#		channel_image = layer_activation[0,:,:, index]
		#		channel_image -= channel_image.mean()
		#		channel_image /= channel_image.std()
		#		channel_image *= 64
		#		channel_image += 128
		#		channel_image = np.clip(channel_image, 0, 255).astype('uint8')
		#		display_grid[col * size : (col + 1) * size, row * size : (row + 1) * size] = channel_image
				
		
		# - Draw plots
		scale = 1. / size
		plt.figure(figsize=(scale * display_grid.shape[1],scale * display_grid.shape[0]))
		plt.title(layer_name)
		plt.grid(False)
		plt.imshow(display_grid, aspect='auto', cmap='viridis')
		plt.show()

##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	#===========================
	#==   PARSE ARGS
	#===========================
	print("INFO: Get script args ...")
	try:
		args= get_args()
	except Exception as ex:
		print("ERROR: Failed to get and parse options (err=%s)",str(ex))
		return 1

	# - Input file
	inputimg= args.inputimg
	model_filename= args.model
	normmin= args.normdatamin
	normmax= args.normdatamax
	normalize_inputs= args.normalize_inputs

	labels_loss_weight= 1
	spars_loss_weight= 1

	#===========================
	#==   LOAD MODEL
	#===========================
	print("INFO: Loading nn model from file %s ..." % model_filename)
	model= load_model(model_filename)

	model.summary()

	#print("INFO: Compiling network ...")
	#opt= optimizers.RMSprop(lr=1.e-4)
	#losses = {
	#	"type": "binary_crossentropy",
	#	"pars": "mse"
	#}
	#lossWeights = {
	#	"type": labels_loss_weight,
	#	"pars": spars_loss_weight
	#}
	#model.compile(optimizer=opt,loss=losses, loss_weights=lossWeights, metrics=['accuracy'])

	#===========================
	#==   READ INPUTS
	#===========================
	print("INFO: Reading input image %s ..." % inputimg)
	data= read_img(inputimg)

	if normalize_inputs:
		print("INFO: Normalizing inputs in range [%s,%s] ..." % (normmin,normmax))
		data= (data - normmin)/(normmax-normmin)

	#nchannels= 1
	#data= data.reshape(data.shape[0],data.shape[1],nchannels)

	#- Set inputs
	#input_data= []
	#input_data.append(data)
	
	#inputs= np.array(input_data)
	#inputs= inputs.astype('float32')
	#if normalize_inputs:
	#	print("INFO: Normalizing inputs in range [%s,%s] ..." % (normmin,normmax))
	#	inputs= (inputs - normmin)/(normmax-normmin)

	# - Set input image tensor
	img_tensor = image.img_to_array(data)
	img_tensor = np.expand_dims(img_tensor, axis=0)
	print(img_tensor.shape)

	# - Draw input image
	plt.imshow(data)
	plt.show()

	#===========================
	#==   GET ACTIVATIONS
	#===========================
	print("INFO: Retrieving model activations for all layers ...")
	#model_activations= get_activations(model, inputs)

	#print model_activations

	get_activations(model,img_tensor)


###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

