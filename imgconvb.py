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
#from numpy import fft

#import pywt
from pywt import idwt2, dwt2, dwtn, dwt, idwt, Wavelet, MODES, dwt_max_level
from pywt.numerix import transpose, array, as_float_array, default_dtype, apply_along_axis
import pywttest_helper as wthlp

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


#import ImageMath
#from numpy import array
#from math import sqrt
from math import ceil
from pywt import upcoef
def all_permutations(symbols, length):
    if length < 1:
        return []
    result = symbols
    for length in range(1, length):
        result = [permutation + symbol for permutation in result for symbol in symbols]
    return result

def _upcoef(data, wavelet, type):
    """Adapts pywt.downcoef call for apply_along_axis"""
    return upcoef(type, data, wavelet, level=1)

def idwtn(data, wavelet):
    """
    Single-level n-dimensional Discrete Wavelet Transform.

    data     - n-dimensional array
    wavelet - wavelet to use (Wavelet object or name string)
    mode    - signal extension mode, see MODES

    Results are arranged in a dictionary, where key specifies
    the transform type on each dimension and value is a n-dimensional
    coefficients array.

    (a=L d=H)

    For example, for a 2D case the result will look something like this:
        {
            'aa': <coeffs>  # A(LL) - approx. on 1st dim, approx. on 2nd dim
            'ad': <coeffs>  # H(LH) - approx. on 1st dim, det. on 2nd dim
            'da': <coeffs>  # V(HL) - det. on 1st dim, approx. on 2nd dim
            'dd': <coeffs>  # D(HH) - det. on 1st dim, det. on 2nd dim
        }
    """
    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    #print "IDWTN:", data
    shape = len(data.keys()[0])
    #print data
    #print "Shape:", shape, data.keys()
    prev = data
    new = [{}]
    #for axis in range(shape):
    #    print "Axis:", axis
    for dimension in range(shape, 0, -1):
        new = {}
        permutations = all_permutations(['a','d'], dimension-1) if dimension > 1 else ['']
        for outer in permutations:
            #print "Outer:", dimension, outer, prev.keys()
#            new_coeffs = []
            #approx = data[''] if outer=='a'*(shape-1) else prev[outer+'a']
            new[outer] = (
                    apply_along_axis(_upcoef, dimension-1, prev[outer+'a'], wavelet, 'a')+
                    apply_along_axis(_upcoef, dimension-1, prev[outer+'d'], wavelet, 'd')
                )
            #print "Shapes:", [array.shape for array in new_coeffs]
        prev = new
    return prev['']

def wavedecn(data, wavelet, mode='sym', level=None, shape=2):
    """
    Multilevel ND Discrete Wavelet Transform.

    data    - 2D input data
    wavelet - wavelet to use (Wavelet object or name string)
    mode    - signal extension mode, see MODES
    level   - decomposition level. If level is None then it will be
              calculated using `dwt_max_level` function .

    Returns coefficients list - [cAn, (cHn, cVn, cDn), ... (cH1, cV1, cD1)]
    """

    data = as_float_array(data)

    if len(data.shape) != shape:
        raise ValueError("Expected {0}D input data.".format(shape))

    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    if level is None:
        size = min(data.shape)
        level = dwt_max_level(size, wavelet.dec_len) # Not working properly for n
    elif level < 0:
        raise ValueError(
            "Level value of %d is too low . Minimum level is 0." % level)

    coeffs_list = []

    a = data
    for i in range(level):
        res = dwtn(a, wavelet, mode)
        #print "Keys:", array(a).shape, i, res.keys()
        a = res['a'*shape]
        coeff_permutations = filter(lambda i: i != 'a'*shape, all_permutations(['d','a'], shape))
        coeffs_list.append({ key : res[key] for key in coeff_permutations})

    coeffs_list.append({'':a})
    coeffs_list.reverse()

    return coeffs_list


def waverecn(coeffs, wavelet, shape=2):
    """
    Multilevel ND Inverse Discrete Wavelet Transform.

    coeffs  - coefficients list [cAn, (cHn, cVn, cDn), ... (cH1, cV1, cD1)]
    wavelet - wavelet to use (Wavelet object or name string)

    Returns 2D array of reconstructed data.
    """

    if not isinstance(coeffs, (list, tuple)):
        raise ValueError("Expected sequence of coefficient arrays.")

    if len(coeffs) < 2:
        raise ValueError(
            "Coefficient list too short (minimum 2 arrays required).")

    a = coeffs[0]['']
    ds = coeffs[1:]

    for i, d in enumerate(ds):
        #print "DNextlevel", i
        # Crop out center from padded 
        cropped = (a if not d or a.shape == d['d'*shape].shape else 
                   a[[slice(int(axe/2.0 - dxe/2.0),
                      int(axe/2.0 + dxe/2.0)) 
                      for axe, dxe in zip(a.shape, d['d'*shape].shape)]])
        td = {'a'*shape:cropped}
        td.update(d)
        #print "Send TD:", [key+str([k+":"+str(elm.shape) for k,elm in val.iteritems()]) for key, val in td.iteritems()]
        #print "Send TD:", [key+":"+str(val.shape) for key, val in td.iteritems()]
        a = idwtn(td, wavelet)
        #TODO: Needs to be cropped back to unpadded versions to match the other blocks in parent transform
    return a

def save2d(prefix, idx, cidx, part, wc, bandmod):
    #wcf = Image.fromarray(wc*1000)
    #print "Save2d:", wc.shape
    wcf = Image.fromarray(bandmod(wc, idx))
    wci = wcf.convert("I")
    wci.convert("RGB").save(prefix.format(idx, cidx, part))

def saveiwtn(prefix, idx, iwt, mod=lambda x:x, bandmod = lambda x,y:x):
    for name, wc in iwt.iteritems():
        if len(wc.shape) < 2:
            raise ValueError("Too small shape")
        if len(wc.shape) == 2:
            save2d(prefix, idx, '', mod(wc), bandmod)
        #TODO: WORKS FOR 3D, DOES THIS WORK FOR 4D? Sceptic
        stack = [0]
        while len(wc[tuple(stack)].shape)>2:
            stack.append(0)
        #print "SWS:", stack,";", wc.shape
        while True:
            #print "Save2d:", name, stack
            save2d(prefix, idx, stack, name, mod(wc[tuple(stack)]), bandmod)
            #print "range(",len(stack),",-2,0,-1)"
            for idx in range(len(stack)-1,-1,-1):
                #print "FLIDX:", idx
                if stack[idx] < len(wc[tuple(stack[0:idx])])-1:
                    stack[idx] +=1
                    break #for loop
                else:
                    stack[idx] = 0
            else:
                #all elements traversed, finished
                break #while loop

def wavesaven(prefix, coeffs, shape=2, maxlevel=-1, mod=lambda x:x, bandmod = lambda x,y:x):
    sliced_coeffs = coeffs if maxlevel == -1 else coeffs[0:maxlevel+1]
    for i, wc in enumerate(sliced_coeffs):
        saveiwtn(prefix, i, wc, mod, bandmod)
 
#def iswt2(data, wavelet, level, start_level=0):
def savearray(array, name, target="RGB"):
    wcf = Image.fromarray((array).clip(0.0,255.0).astype('uint8'))
    wcf.convert(target).save(name, "PNG")

def savemultilevel(wavedec, idx = 0, level = 0, indices = []):
    if not hasattr(wavedec[0][0], '__iter__'):
        savearray(wavedec, "wc{0}Leaf.png".format(indices), "L")
        return
    approx = wavedec[0]
    details = wavedec[1:]
    savearray(approx, "wc{0}.png".format(indices), "L")
    for idx, detail in enumerate(details):
        savemultilevel(detail, idx, level+1, indices+[idx]) 

def extractplane(plane, data, depth):
    bits = data.astype('uint16') & (numpy.ones(data.shape, dtype=numpy.uint8)*plane)
    return numpy.piecewise(data, (bits == plane, bits == 0), (plane, 0))
    #return data

def diffenc(data):
    shaped = numpy.reshape(data,-1)
    convolved = numpy.convolve(shaped, [1,-1])[0:len(shaped)]
    diffenced = numpy.reshape(convolved, data.shape)
    return diffenced

def diffdec(data):
    shaped = numpy.reshape(data,-1)
    convolved = numpy.add.accumulate(shaped)
    diffenced = numpy.reshape(convolved, data.shape)
    return diffenced

def findlimit(data, part = 0.5):
    #find limit for the lowest part so that 
    #len(data < limit) = len(data) * part
    shaped = numpy.reshape(data,-1)
    l = len(shaped)
    return numpy.sort(shaped)[::-1][(l*part)-1 if l*part < l-1 else l-1]
def findlimitf(data, part = 0.5):
    #find limit for the lowest part so that 
    #len(data < limit) = len(data) * part
    shaped = numpy.reshape(data,-1)
    l = len(shaped)
    return numpy.sort(shaped)[(l*part)-1 if l*part < l-1 else l-1]


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
    
#yuvimgs = array([127.0+(127.0*numpy.dot(image / 350.0, rgbtoyuv.T)) for image in images])
#yuv = [numpy.zeros(yuvimgs.shape[0:3])]*3
#for i in xrange(3):
#    yuv[i] = yuvimgs[:,:,:,i]

#yuv = [ yuvimgs[:,:,:, i].astype('uint8') for i in range(3)]
#yuvc = reduce(lambda t,e: numpy.concatenate((t, e)), yuv)
#yuvc.tofile("input.raw")
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
