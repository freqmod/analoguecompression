#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Copyright 2013 Frederik M.J. Vestre AGplv3

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
#BW
numframes = 32
if False:
    rev = numpy.fromfile("output.raw", dtype="int8").reshape([1,numframes,128,256])
    #rev = numpy.fromfile("input.raw", dtype="uint8").reshape([1,numframes,128,256])
    revmerged = numpy.zeros([numframes,128,256,1])
    for i in range(1):
        revmerged[:,:,:,i] = (rev[i].clip(-127.0, 127.0)+127.0)    #int16_t
        #revmerged[:,:,:,i] = (rev[i].clip(0.0, 255.0)-127.0)/127.0    #uint8_t
    for i in range(len(revmerged)):
        savearray(numpy.squeeze(revmerged[i]), "wcFromC{0}.png".format(str(i).zfill(4)), "L");
else:
    #YUV
    rev = numpy.fromfile("output.raw", dtype="int16").reshape([3,numframes,128,256])
    #rev = numpy.fromfile("input.raw", dtype="uint8").reshape([3,numframes,128,256])
    revmerged = numpy.zeros([numframes,128,256,3])
    for i in range(3):
        revmerged[:,:,:,i] = ((rev[i]).clip(-127.0, 127.0))/127.0    #int16_t
        #revmerged[:,:,:,i] = (rev[i].clip(0.0, 255.0)-127.0)/127.0    #uint8_t
    for i in range(len(revmerged)):
        savearray(numpy.dot(revmerged[i], yuvtorgb.T)*375.0, "wcFromC{0}.png".format(str(i).zfill(4)), "RGB");
