#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2006-2012 Filip Wasilewski <http://en.ig.ma/>
# See COPYING for license details.

import optparse
import os
import sys
import operator
import wave
import math
#import pdb
from functools import partial
if os.name == 'nt':
    from time import clock  # noqa
else:
    from time import time as clock  # noqa

import Image  # PIL
import numpy  # http://www.scipy.org
from numpy import array

def image2array(image):
    """PIL Image to NumPy array"""
    assert image.mode in ('L', 'RGB', 'CMYK')
    return numpy.asarray(image)
    #arr = numpy.fromstring(image.tostring(), numpy.uint8)
    #arr.shape = (image.size[1], image.size[0], len(image.getbands()))
    #return arr.swapaxes(0, 2).swapaxes(1, 2).astype(numpy.floatnumframes)


def array2image(arr, mode):
    """NumPy array to PIL Image"""
    #arr = arr.swapaxes(1, 2).swapaxes(0, 2)
    #arr[arr < 0] = 0
    #arr[arr > 255] = 255
    #arr = numpy.fix(arr).astype(numpy.uint8)
    #return Image.fromstring(mode, arr.shape[1::-1], arr.tostring())
    return Image.fromarray(numpy.uint8(arr), mode = mode)


def load_image(path, mode=None, size=None):
    """Load image"""
    im = Image.open(path)

    if im.mode not in ('L', 'P', 'RGB', 'CMYK'):
        raise TypeError("Image mode must be 'L', 'P', 'RGB' or 'CMYK'")

    if mode is not None:
        if mode == 'P':
            raise ValueError("Mode must be 'L', 'RGB' or 'CMYK'")
        im = im.convert(mode)
    elif im.mode == 'P':
        im = im.convert('RGB')

    if size is not None and im.size != size:
        im = im.resize(size, Image.ANTIALIAS)
    return im

def savearray(array, name, target="RGB"):
    wcf = Image.fromarray((array).clip(0.0,255.0).astype('uint8'))
    wcf.convert(target).save(name, "PNG")

#images  = array([numpy.array(Image.open("{0}.png".format(str(i).zfill(2))).convert('RGB')) for i in range(0,81)])
#images  = array([numpy.array(Image.open("{0}.png".format(str(i).zfill(2))).convert('RGB')) for i in range(14,22)])
#images  = array([numpy.array(Image.open("{0}.png".format(str(i).zfill(2))).convert('RGB')) for i in range(0,80)])
#images  = array([numpy.array(Image.open("psy/{0}.png".format(str(i).zfill(8))).convert('RGB')) for i in range(1,50)])
#images  = array([numpy.array(image2array(Image.open("ASDSpin/{0}.png".format(str(i).zfill(8))).convert('RGB'))) for i in range(1,33)])
#images  = array([numpy.array(image2array(Image.open("mito/{0}.png".format(str(i).zfill(8))).convert('RGB'))) for i in range(1,33)])
#images  = array([numpy.array(image2array(Image.open("../analoguecompression/py/{0}.png".format(str(i).zfill(2))).convert('RGB'))) for i in range(5, 37)])
images  = array([numpy.array(image2array(Image.open("psy/{0}.png".format(str(i).zfill(8))).convert('RGB'))) for i in range(1,33)]) # split: parking, lift
#images  = array([numpy.array(image2array(Image.open("psy/{0}.png".format(str(i).zfill(8))).convert('RGB'))) for i in range(1,17)]) #before lift # NEEDS SPECIAL CODE
#images  = array([numpy.array(image2array(Image.open("psy/{0}.png".format(str(i).zfill(8))).convert('RGB'))) for i in range(15,47)]) # lift scene
#images  = array([numpy.array(image2array(Image.open("wom/{0}.png".format(str(i).zfill(2))).convert('RGB'))) for i in range(1,33)]) # wom: Clip
##images  = array([numpy.array(image2array(Image.open("psy/{0}.png".format(str(i).zfill(8))).convert('RGB'))) for i in range(1,3)])
print "Loaded", len(images)
if(len(images)!=32):
    print "WARNING: Image length not 32, will cause problems in compressor"
#images  = array([numpy.array(Image.open("{0}.png".format(str(i).zfill(2))).convert('L')) for i in range(12,20)])
#images  = array([numpy.array(Image.open("{0}.png".format(str(i).zfill(2))).convert('RGB')) for i in range(12,28)])
#images  = array([numpy.array(Image.open("{0}.png".format(str(i).zfill(2))).convert('RGB')) for i in range(12,14)])
#From somewhere on stackoverflow
#rgbtoyuv = array([[ 0.18266261,  0.61447292,  0.06197099],                                                                                                                                                                      
#                  [-0.10067201, -0.33865837,  0.43933038],                                                                                                                                                                      
#                  [ 0.4391415 , -0.39891048, -0.04023103]])  
#Inverse
#yuvtorgb = array([[1.164,     0.,  1.793],
#           [1.164, -0.213, -0.533],
#           [1.164,  2.112,     0.]])
#From wikipedia
rgbtoyuv = array([[ 0.0126    , 0.7152     , 0.7222],                                                                                                                                                                      
                  [ -0.09991  , -0.33609   , 0.436 ],
                  [ 0.615     , -0.55861   , -0.05639]])
#rgbtoyuv = array([[ 0.21259997,  0.71520041,  0.07219962],
#                  [-0.09990694, -0.33609358,  0.43600052],
#                  [ 0.61499772, -0.5586063 , -0.05639142]])
yuvtorgb = array([[ 1   , 0         ,  1.28033 ],
                  [ 1   , -0.21482  , -0.38059 ],
                  [ 1   ,  2.12798  ,      0.0 ]])
#http://stackoverflow.com/questions/7041172/pils-colour-space-conversion-ycbcr-rgb
#yuvtorgb[:,0]  *= 255./219.
#yuvtorgb[:,1:] *= 255./112.
    
yuvimgs = array([127.0+(127.0*numpy.dot(image / 375.0, rgbtoyuv.T)) for image in images])
#yuv = [numpy.zeros(yuvimgs.shape[0:3])]*3
#for i in xrange(3):
#    yuv[i] = yuvimgs[:,:,:,i]

yuv = [ yuvimgs[:,:,:, i].astype('uint8') for i in range(3)]
yuvc = reduce(lambda t,e: numpy.concatenate((t, e)), yuv)
yuvc.tofile("input.raw")
#revmerged = numpy.zeros(list(yuv[0].shape)+[3])
#for i in range(3):
#    revmerged[:,:,:,i] = (yuv[i].clip(0.0, 255.0)-127.0)/127.0    #for i in range(len(revmerged)):
#for i in range(len(revmerged)):
#    savearray(numpy.dot(revmerged[i], yuvtorgb.T)*350.0, "wcReversed{0}.png".format(str(i).zfill(4)), "RGB");
