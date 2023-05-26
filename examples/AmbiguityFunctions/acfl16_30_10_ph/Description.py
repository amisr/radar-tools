#!/usr/bin/env python

# install https://github.com/amisr/radar-tools on python3
from radartools import amb_func_tools
import numpy as np

a16code,signs=amb_func_tools.a16rand()
amb_func_tools.compute_lamb(ltype=4, hfile='blackman_10.00usec_051010.fco', MCIC=[5,10],
        in1=np.transpose(signs), in2=[30,3,'r20101112_008b_480.txt'], outdir="", lags=[])
