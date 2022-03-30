#!/usr/bin/env python

#################################################################################$$
# Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$

# Usage PerMoleculeIPDRatio.py <bamIn> <bamOut>
# model used is hardcoded in the script
# input bam is an unalined HiFi kenetics bam processed with ccs-kinetics-bystrandify https://ccs.how/faq/kinetics.html
# output bam is an unalinged bam with f array tag indicating ipd ratio
# code is single threaded, split the unaligned bam for parallel execution


# Conda environment is the best way to deal with dependencies. 
# Everything should be available through conda / bioconda, in a python 3 environment
import numpy as np
from pbcore.io import BamReader, IndexedBamReader, IndexedFastaReader, AlignmentSet
import pickle
import sys
import itertools
import pandas as pd
import gzip
import pysam
from array import array 
from tqdm import tqdm, trange

# Kinetics tools can be found https://github.com/PacificBiosciences/kineticsTools 
# and can be installed in to a conda environment using setup.py
#from kineticsTools.KineticWorker import KineticWorkerProcess
#from kineticsTools.ResultWriter import KineticsWriter
from kineticsTools.ipdModel import IpdModel, GbmContextModel
#from kineticsTools import ReferenceUtils, loader
from kineticsTools.sharedArray import SharedArray

bamIn = pysam.AlignmentFile(sys.argv[1], check_sq=False)
bamOut = pysam.AlignmentFile(sys.argv[2], "wb", header=bamIn.header)

def snippetFunc(refId, pre, post):
    """
    Return a function that returns a snippet of the reference sequence around a given position
    """

    refArray = refId.getNumpyWrapper()

    def f(tplPos, tplStrand):
        """Closure for returning a reference snippet. The reference is padded with N's for bases falling outside the extents of the reference"""
        # skip over the padding
        tplPos += pad

        # Forward strand
        if tplStrand == 0:
            slc = refArray[(tplPos - pre):(tplPos + 1 + post)]
            slc = np.right_shift(slc, 4)
            return "".join(c for c in seqMapNp[slc])

        # Reverse strand
        else:
            slc = refArray[(tplPos + pre):(tplPos - post - 1):-1]
            slc = np.right_shift(slc, 4)
            return "".join(c for c in seqMapComplementNp[slc])

    return f

def _makeFramepoints():
    B = 2
    t = 6
    T = 2**t

    framepoints = []
    next = 0
    for i in range(256//T):
        grain = B**i
        nextOnes = next + grain * np.arange(0, T)
        next = nextOnes[-1] + grain
        framepoints = framepoints + list(nextOnes)
    return np.array(framepoints, dtype=np.uint16)

def _makeLookup(framepoints):
    # (frame -> code) involves some kind of rounding
    # basic round-to-nearest
    frameToCode = np.empty(shape=max(framepoints)+1, dtype=int)
    for i, (fl, fu) in enumerate(zip(framepoints, framepoints[1:])):
        if (fu > fl + 1):
            m = (fl + fu)//2
            for f in range(fl, m):
                frameToCode[f] = i
                frameToCode[f] = i + 1
        else:
            frameToCode[fl] = i
    # Extra entry for last:
    frameToCode[fu] = i + 1
    return frameToCode, fu

_framepoints = _makeFramepoints()
_frameToCode, _maxFramepoint = _makeLookup(_framepoints)


def framesToCode(nframes):
    nframes = np.minimum(_maxFramepoint, nframes)
    return _frameToCode[nframes]

def codeToFrames(code):
    return _framepoints[code]

byte = np.dtype('byte')
uint8 = np.dtype('uint8')

# Map for ascii encoded bases to integers 0-3 -- will be used to define a 24-bit lookup code
# for fetching predicted IPDs from the kinetic LUT.

# We start everything at 0, so anything will map to 'A' unless it appears
# in this table
lutCodeMap = np.zeros(256, dtype=uint8)
maps = {'a': 0, 'A': 0, 'c': 1, 'C': 1, 'g': 2, 'G': 2, 't': 3, 'T': 3}
for k in maps:
    lutCodeMap[ord(k)] = maps[k]
lutReverseMap = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

seqCodeMap = np.ones(256, dtype=uint8) * 4
for k in maps:
    seqCodeMap[ord(k)] = maps[k]
seqMap = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}
seqMapNp = np.array(['A', 'C', 'G', 'T', 'N'])

seqMapComplement = {0: 'T', 1: 'G', 2: 'C', 3: 'A', 4: 'N'}
seqMapComplementNp = np.array(['T', 'G', 'C', 'A', 'N'])

# Base letters for modification calling
# 'H' : m6A, 'I' : m5C, 'J' : m4C, 'K' : m5C/TET
baseToCode = {'N': 0, 'A': 0, 'C': 1, 'G': 2,
              'T': 3, 'H': 4, 'I': 5, 'J': 6, 'K': 7}
baseToCanonicalCode = {'N': 0, 'A': 0, 'C': 1,
                       'G': 2, 'T': 3, 'H': 0, 'I': 1, 'J': 1, 'K': 1}

codeToBase = dict([(y, x) for (x, y) in baseToCode.items()])

pre = 10
post = 4
pad = 30
base4 = 4 ** np.array(range(pre + post + 1))
refDict = {}
refLengthDict = {}

refid = 0

lutPath = "./kineticsTools/resources/SP3-C3.npz.gz"

with gzip.open(lutPath, "rb") as npz_in:
    gbmModelData = np.load(npz_in, allow_pickle=True)
    gbmModel = GbmContextModel(gbmModelData, -1)

for read in tqdm(bamIn, desc='Reads'):
    rawSeq = read.seq
    refSeq = np.frombuffer(rawSeq.encode("utf-8"), dtype=byte)
  
    # Store the reference length
    
    length = len(rawSeq)

    # Make a shared array
    sa = SharedArray(dtype='B', shape=len(rawSeq) + pad * 2)
    saWrap = sa.getNumpyWrapper()

    # Lut Codes convert Ns to As so that we don't put Ns into the Gbm Model
    # Seq Codes leaves Ns as Ns for getting reference snippets out
    innerLutCodes = lutCodeMap[refSeq]
    innerSeqCodes = seqCodeMap[refSeq]
    innerCodes = np.bitwise_or(innerLutCodes, np.left_shift(innerSeqCodes, 4))

    saWrap[pad:(len(rawSeq) + pad)] = innerCodes

    # Padding codes -- the lut array is padded with 0s the sequence
    # array is padded with N's (4)
    outerCodes = np.left_shift(np.ones(pad, dtype=uint8) * 4, 4)
    saWrap[0:pad] = outerCodes
    saWrap[(len(rawSeq) + pad):(len(rawSeq) + 2 * pad)] = outerCodes

    snipFunction = snippetFunc(sa, post, pre)
    sites = range(0,length)
    contexts = [snipFunction(sites[x], 1) for x in sites]
    
    control = gbmModel.getPredictions(contexts)
    ip = read.get_tag('ip')
    ipFrames = codeToFrames(read.get_tag('ip'))
    ipFramesNorm = ipFrames / np.mean(ipFrames)
    ipr = ipFramesNorm / control
    
    read.set_tag('in', array('f', ipFramesNorm))
    read.set_tag('ir', array('f', ipr))
    bamOut.write(read)
    
bamOut.close()