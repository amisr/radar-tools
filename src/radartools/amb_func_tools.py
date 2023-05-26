#! /usr/bin/env python

"""
xxxxx

~M. Nicolls
last revised: xx/xx/2007

"""

import sys
import os.path
import h5py # importing to avoid weird matplotlib interference with tables
import tables
import numpy as np
import scipy, scipy.fftpack
import scipy.io
import glob
#import matplotlib
#matplotlib.use('TkAgg'); matplotlib.interactive(True)
#matplotlib.use('Agg');# matplotlib.interactive(True)
import matplotlib.pyplot as plt
#import pylab
import string

##############################

v_lightspeed=299792458

def conv2d(a, b):
    x = np.zeros((a.shape[0], a.shape[1] + b.shape[1]-1), dtype=a.dtype)
    a = np.fliplr(a)
 
    for ii in range(a.shape[0]):
        x[ii,:] = np.convolve(a[ii,:], b[ii,:])

    return x

def baudexpand(x, baud):

    if baud.size==1:
        px = np.reshape(np.repeat(x[:,np.newaxis], int(baud), 0), (1,x.size*int(baud)))
    else:
        px = np.reshape(baud*np.transpose(x), (1,x.size*baud.size))
        xxxx
        
    return px

def a4():

    codeset=[[1,1,1,1],
        [1,-1,1,1],
        [1,1,-1,1],
        [1,-1,-1,1],
        [1,1,1,-1],
        [1,-1,1,-1],
        [1,1,-1,-1],
        [1,-1,-1,-1]]
        
    return np.array(codeset)
    
def code88():
    codec=[1,-2,1,-5,3,-1,2,-2,1,-1,2,-3,1,-2,2,-1,1,-1,4,-1,1,-1,1,-1,3,-6,1,-2,1,-1,1,-6,1,-1,2,-2,3,-3,1,-1,2,-2,3,-3,1,-1]
    
    code=[]
    for ii in range(len(codec)):
        for jj in range(np.absolute(codec[ii])):
            code.append(np.sign(codec[ii]))
    
    return np.array(code)


def a2():
    # Calculates randomized set of A16 codes, strong condition
    #
    # Octal
    aoct = ['000','000','001']
    a=[]
    for ii in range(len(aoct)):
        a.append(string.atoi(aoct[ii],8))
    
    a2code=np.zeros(2,dtype='float64')
    signs=np.zeros((2,3),dtype='float64')
    for jpul in range(2):
        codetmp = 0
        signstmp = []
        for ibaud in range(len(a)):
            if walsh(a[ibaud],jpul)==1:
                signstmp.append(1)
                codetmp = codetmp + np.power(2,(len(aoct)-1-ibaud))
            else:
                signstmp.append(-1)
        a2code[jpul] = codetmp
        print (signstmp)
        signs[jpul,:] = signstmp

    return a2code,signs

def a16rand():
    # Calculates randomized set of A16 codes, strong condition
    #
    # Octal
    codes_oct = ['000','001','002','004','010','020','017','037','021','014','031','035','024','006','015','032']
    codes=[]
    for code in codes_oct:
        codes.append(int(code,8))

    randomizer = [1,1,1,1,-1,-1,-1,1,-1,-1,1,-1,-1,1,-1,1]

    a16code=np.zeros(32,dtype='float64')
    signs=np.zeros((32,16),dtype='float64')
    for jpul in range(32):
        codetmp = 0
        signstmp = []
        for ibaud,(randomizer_i, code_i) in enumerate(zip(randomizer,codes)):
            if randomizer_i*walsh(code_i,jpul)==1:
                signstmp.append(1)
                codetmp = codetmp + np.power(2,(len(codes) - 1 - ibaud))
            else:
                signstmp.append(-1)
        a16code[jpul] = codetmp
        signs[jpul,:] = signstmp

    return a16code,signs
    
def a8rand():
    # Calculates randomized set of A8 codes, strong condition
    #
    # Octal
    aoct = ['000','001','002','004','010','003','007','016']
    a=[]
    for ii in range(len(aoct)):
        a.append(string.atoi(aoct[ii],8))

    randomizer = [1,-1,1,-1,1,-1,-1,1]
    
    a16code=np.zeros(16,dtype='float64')
    signs=np.zeros((16,8),dtype='float64')
    for jpul in range(16):
        codetmp = 0
        signstmp = []
        for ibaud in range(8):
            if randomizer[ibaud]*walsh(a[ibaud],jpul)==1:
                signstmp.append(1)
                codetmp = codetmp + np.power(2,(7-ibaud))
            else:
                signstmp.append(-1)
        a16code[jpul] = codetmp
        signs[jpul,:] = signstmp

    return a16code,signs

# walsh
def walsh(i,j):
    # calculates elements of a walsh sign matrix for alternating codes
    bita = np.bitwise_and(i,j)
    s = 1
    for ibit in range(16):
        if np.bitwise_and(bita,np.power(2,ibit))>0:
            s = -s
    return s

# generate_barkercode
def generate_barkercode(leng,opt=0):
    # returns barker code pulse sequence

    if leng==1:
        c=[1]
    elif (leng==2) and (opt==0):
        c=[1,-1]
    elif (leng==2) and (opt==1):
        c=[1,1]
    elif leng==3:
        c=[1,1,-1]
    elif (leng==4) and (opt==0):
        c=[1,-1,1,1]
    elif leng==4 and (opt==1):
        c=[1,-1,-1,-1]
    elif leng==5:
        c=[1,1,1,-1,1]
    elif leng==7:
        c=[1,1,1,-1,-1,1,-1]
    elif leng==11:
        c=[1,1,1,-1,-1,-1,1,-1,-1,1,-1]
    elif leng==13:
        c=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1]
    elif leng==28:
#        c=[1,1,-1,1,1,-1,1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,-1,-1]
        c=[1,-1,-1,-1,1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,1]
#		1,-1,1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,-1,-1]

    return c

def multipulse(gaps, pattern=[]):

    mp=[]
    
    if len(pattern)==0:
        mp.append(1.0)
        for igap in range(len(gaps)):
            gap=gaps[igap]
            mp.extend(np.zeros(gap).tolist())
            mp.append(1.0)
        
    else: # pattern is repeated without extra breaks
        mp.extend(pattern)
        leng = len(pattern)
        for igap in range(len(gaps)):
            gap=gaps[igap]
            mp.extend(np.zeros(gap).tolist())
            mp.extend(pattern)
   
    return np.array(mp)

# mpulse
def dpulse(codeset, clen, baud, os, h, lags):
    #
    # computes a barker code range-lag ambiguity function
    #
    # codeset - coding of pulse
    # baud - baud length (microseconds)
    # os - oversampling factor
    # h - impulse response function

    baud=np.array(baud)
    
    print (codeset)
    xxxxx
    
    hbc=np.squeeze(np.fliplr(baudexpand(codeset, baud))) # this is the total impulse response of a Barker decoder with delta function impulse response
    env=np.squeeze(baudexpand(multipulse(pattern,codeset),baud)) # pulse transmission envelope 
    
    # sampling of barker code
    sc=np.zeros(hbc.shape,hbc.dtype)
    sc[::baud/os]=1
    hbc=hbc*sc

    # delta function
    hdelt=np.zeros(h.shape)
    hdelt[np.floor(hdelt.shape[0]/2.0)]=1.0
        
    # to get the total impulse response, convolve filter impulse response with decoder
    h=np.convolve(h,hbc)
    hDelt=np.convolve(hdelt,hbc)
    
    wta=np.zeros((h.shape[0],env.shape[0]+h.shape[0]),dtype='float32')
    wtaDelt=np.zeros((hDelt.shape[0],env.shape[0]+h.shape[0]),dtype='float32')

    for tau in range(h.shape[0]):
        wta[tau,tau:tau+env.shape[0]]=h[tau]*env
        wtaDelt[tau,tau:tau+env.shape[0]]=hDelt[tau]*env

    mlag=np.array(lags).max()

    wtt=np.zeros((baud*mlag+2*h.shape[0],env.shape[0]+h.shape[0]),dtype='float32')
    wlag=np.zeros((len(lags),baud*mlag+2*h.shape[0]),dtype='float32')
    wrng=np.zeros((len(lags),env.shape[0]+h.shape[0]),dtype='float32')
    try:
        wttall=np.zeros((len(lags),baud*mlag+2*h.shape[0],env.shape[0]+h.shape[0]),dtype='float32')
        comp_all=1
    except:
        print ('Unable to allocate memory for full, 3D ambiguity function')
        comp_all=0
        wttall=[]
        
    offs = np.arange(2*h.shape[0]-1)

    for ii in range(len(lags)):
        lgind=lags[ii]
        lag = baud*lgind
        print ('Lag-%d' % (lag))
        if lgind==0:
            tmp=np.transpose(conv2d(np.transpose(wta[:,lag:]), np.transpose(wta[:,0:])))
        else:
            tmp=np.transpose(conv2d(np.transpose(wta[:,lag:]), np.transpose(wta[:,0:-lag])))
        if comp_all:
            wttall[ii,lag+offs,lag:]=tmp
        else:
            wlag[ii,lag+offs]=wlag[ii,lag+offs]+np.sum(tmp,axis=1)
            wrng[ii,lag:]=wrng[ii,lag:]+np.sum(tmp,axis=0)
        tmp2=wtt[lag+offs,lag:]
        I=np.where(tmp>tmp2)
        tmp2[I]=tmp[I]
        wtt[lag+offs,lag:]=tmp2

    wttDelt=np.transpose(conv2d(np.transpose(wtaDelt), np.transpose(wtaDelt)))
    fact=np.sum(np.sum(wttDelt))
    
    tau=np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r=np.arange(wtt.shape[1],dtype='float64')

    wtt=wtt/fact #/(leng*leng*baud*os*os)
    wttall=wttall/fact #/(leng*leng*baud*os*os)

    # lag ambiguity function
    wlag=np.sum(wttall,axis=2)
    wlagsum=np.sum(wlag,axis=1)
    
    # range ambiguity function
    wrng=np.sum(wttall,axis=1)
    wrngsum=np.sum(wrng,axis=1)

    comp_all=1

    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all

# mpulse
def mpulse(codeset, pattern, baud, os, h, lags):
    #
    # computes a barker code range-lag ambiguity function
    #
    # codeset - coding of pulse
    # pattern - multipulse pattern
    # baud - baud length (microseconds)
    # os - oversampling factor
    # h - impulse response function

    baud=np.array(baud)
    
    print (multipulse(pattern,codeset))
    print (os)
    
    hbc=np.squeeze(np.fliplr(baudexpand(codeset, baud))) # this is the total impulse response of a Barker decoder with delta function impulse response
    env=np.squeeze(baudexpand(multipulse(pattern,codeset),baud)) # pulse transmission envelope 
    
    # sampling of barker code
    sc=np.zeros(hbc.shape,hbc.dtype)
    sc[::baud/os]=1
    hbc=hbc*sc

    # delta function
    hdelt=np.zeros(h.shape)
    hdelt[np.floor(hdelt.shape[0]/2.0)]=1.0
        
    # to get the total impulse response, convolve filter impulse response with decoder
    h=np.convolve(h,hbc)
    hDelt=np.convolve(hdelt,hbc)
    
    wta=np.zeros((h.shape[0],env.shape[0]+h.shape[0]),dtype='float32')
    wtaDelt=np.zeros((hDelt.shape[0],env.shape[0]+h.shape[0]),dtype='float32')

    for tau in range(h.shape[0]):
        wta[tau,tau:tau+env.shape[0]]=h[tau]*env
        wtaDelt[tau,tau:tau+env.shape[0]]=hDelt[tau]*env

    mlag=np.array(lags).max()

    wtt=np.zeros((baud*mlag+2*h.shape[0],env.shape[0]+h.shape[0]),dtype='float32')
    wlag=np.zeros((len(lags),baud*mlag+2*h.shape[0]),dtype='float32')
    wrng=np.zeros((len(lags),env.shape[0]+h.shape[0]),dtype='float32')
    try:
        wttall=np.zeros((len(lags),baud*mlag+2*h.shape[0],env.shape[0]+h.shape[0]),dtype='float32')
        comp_all=1
    except:
        print ('Unable to allocate memory for full, 3D ambiguity function')
        comp_all=0
        wttall=[]
        
    offs = np.arange(2*h.shape[0]-1)

    for ii in range(len(lags)):
        lgind=lags[ii]
        lag = baud*lgind
        print ('Lag-%d' % (lag))
        if lgind==0:
            tmp=np.transpose(conv2d(np.transpose(wta[:,lag:]), np.transpose(wta[:,0:])))
        else:
            tmp=np.transpose(conv2d(np.transpose(wta[:,lag:]), np.transpose(wta[:,0:-lag])))
        if comp_all:
            wttall[ii,lag+offs,lag:]=tmp
        else:
            wlag[ii,lag+offs]=wlag[ii,lag+offs]+np.sum(tmp,axis=1)
            wrng[ii,lag:]=wrng[ii,lag:]+np.sum(tmp,axis=0)
        tmp2=wtt[lag+offs,lag:]
        I=np.where(tmp>tmp2)
        tmp2[I]=tmp[I]
        wtt[lag+offs,lag:]=tmp2

    wttDelt=np.transpose(conv2d(np.transpose(wtaDelt), np.transpose(wtaDelt)))
    fact=np.sum(np.sum(wttDelt))
    
    tau=np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r=np.arange(wtt.shape[1],dtype='float64')

    wtt=wtt/fact #/(leng*leng*baud*os*os)
    wttall=wttall/fact #/(leng*leng*baud*os*os)

    # lag ambiguity function
    wlag=np.sum(wttall,axis=2)
    wlagsum=np.sum(wlag,axis=1)
    
    # range ambiguity function
    wrng=np.sum(wttall,axis=1)
    wrngsum=np.sum(wrng,axis=1)

    comp_all=1

    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all


# barker_code
def coh_code(codeset, baud, os, h):
    #
    # computes a barker code range-lag ambiguity function
    #
    # leng - barker code lenggth
    # baud - baud length (microseconds)
    # os - oversampling factor
    # h - impulse response function

    baud=np.array(baud)

#    codeset=np.array(generate_barkercode(leng,opt=type)) # generate barker code
    
    hbc=np.squeeze(np.fliplr(baudexpand(codeset, baud))) # this is the total impulse response of a Barker decoder with delta function impulse response
    env=np.flipud(hbc) # pulse transmission envelope 
    
    # sampling of barker code
    sc=np.zeros(hbc.shape,hbc.dtype)
    sc[::baud//os]=1
    hbc=hbc*sc

    # delta function
    hdelt=np.zeros(h.shape)
    hdelt[int(np.floor(hdelt.shape[0]/2.0))]=1.0
        
    # to get the total impulse response, convolve filter impulse response with decoder
    h=np.convolve(h,hbc)
    hDelt=np.convolve(hdelt,hbc)

    wta=np.zeros((h.shape[0],env.shape[0]+h.shape[0]),dtype='float64')
    wtaDelt=np.zeros((hDelt.shape[0],env.shape[0]+h.shape[0]),dtype='float64')

    for tau in range(h.shape[0]):
        wta[tau,tau:tau+env.shape[0]]=h[tau]*env
        wtaDelt[tau,tau:tau+env.shape[0]]=hDelt[tau]*env
    
    wtt=np.transpose(conv2d(np.transpose(wta), np.transpose(wta)))
    wttall=wtt[np.newaxis,:,:]

    wttDelt=np.transpose(conv2d(np.transpose(wtaDelt), np.transpose(wtaDelt)))
    fact=np.sum(np.sum(wttDelt))

    tau=np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r=np.arange(wtt.shape[1],dtype='float64')

    wtt=wtt/fact #/(leng*leng*baud*os*os)
    wttall=wttall/fact #/(leng*leng*baud*os*os)

    # lag ambiguity function
    wlag=np.sum(wttall,axis=2)
    wlagsum=np.sum(wlag,axis=1)
    
    # range ambiguity function
    wrng=np.sum(wttall,axis=1)
    wrngsum=np.sum(wrng,axis=1)

    comp_all=1

    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all, h


# barker_code
def barker_code(leng, baud, os, h, Fs, btype=0):
    #
    # computes a barker code range-lag ambiguity function
    #
    # leng - barker code lenggth
    # baud - baud length (microseconds)
    # os - oversampling factor
    # h - impulse response function
    # Fs - nyquist
    # btype - in cases where there are two BCs

    baud=np.array(baud)

    codeset=np.array(generate_barkercode(leng,opt=btype)) # generate barker code
    
    hbc=np.squeeze(np.fliplr(baudexpand(codeset, baud))) # this is the total impulse response of a Barker decoder with delta function impulse response
    env=np.flipud(hbc) # pulse transmission envelope 

    """
    # messed up pulse (additional 30us on end of transmit pulse)
    tcodeset=np.concatenate((codeset,np.ones(3)))  
    env=np.squeeze(baudexpand(tcodeset, baud))
    """
    
    # sampling of barker code
    sc=np.zeros(hbc.shape,hbc.dtype)
    sc[::baud//os]=1
    hbc=hbc*sc/10.0
    print (np.sum(hbc))

    # delta function
    hdelt=np.zeros(h.shape)
    hdelt[int(np.floor(hdelt.shape[0]/2.0))]=1.0
    
    BWfilt=Fs*np.sum(h**2.0)/np.sum(h)**2.0
    print (BWfilt)
    
    tpp=np.absolute(scipy.fftpack.fftshift(scipy.fftpack.fft(h)))**2.0
    tf=np.linspace(-Fs,Fs,tpp.shape[0])
    df=np.median(np.diff(tf))
    
#    BWfilt=0.5*df*np.sum(tpp)**2.0/np.sum(tpp**2.0)
#    print BWfilt
    
#    pylab.plot(tf,10.0*np.log10(tpp/tpp.max()))        
                  
    # to get the total impulse response, convolve filter impulse response with decoder
    h=np.convolve(h,hbc)
    hDelt=np.convolve(hdelt,hbc)
    
    tp=np.absolute(scipy.fftpack.fftshift(scipy.fftpack.fft(h-np.mean(h))))**2.0
    tf=np.linspace(-Fs,Fs,tp.shape[0])
    df=np.median(np.diff(tf))
    
    #pylab.plot(tf,10.0*np.log10(tp/tp.max()))
    #adssdaas
       
    BW=Fs*np.sum(h**2.0)/np.sum(h)**2.0
    print ("BW",BW)

#    BWfilt=0.5*df*np.sum(tp)**2.0/np.sum(tp**2.0)
#    print BWfilt

#    xxxx
    
    wta=np.zeros((h.shape[0],env.shape[0]+h.shape[0]),dtype='float64')
    wtaDelt=np.zeros((hDelt.shape[0],env.shape[0]+h.shape[0]),dtype='float64')

    for tau in range(h.shape[0]):
        wta[tau,tau:tau+env.shape[0]]=h[tau]*env
        wtaDelt[tau,tau:tau+env.shape[0]]=hDelt[tau]*env
    
    wtt=np.transpose(conv2d(np.transpose(wta), np.transpose(wta)))
    wttall=wtt[np.newaxis,:,:]

    wttDelt=np.transpose(conv2d(np.transpose(wtaDelt), np.transpose(wtaDelt)))
    fact=np.sum(np.sum(wttDelt))

    tau=np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r=np.arange(wtt.shape[1],dtype='float64')

    wtt=wtt/fact #/(leng*leng*baud*os*os)
    wttall=wttall/fact #/(leng*leng*baud*os*os)

    # lag ambiguity function
    wlag=np.sum(wttall,axis=2)
    wlagsum=np.sum(wlag,axis=1)
    
    # range ambiguity function
    wrng=np.sum(wttall,axis=1)
    wrngsum=np.sum(wrng,axis=1)

    comp_all=1

    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all, BWfilt, BW

# alt_code_phaseerr
def alt_code_phaserr(codeset, baud, h, lags):
    #
    # computes an alternating code range-lag ambiguity function for the case of a coding error
    #
    # codeset - sign sequence for the alternating code (should be nbaud x nscans size)
    # baud - baud length (microseconds)
    # h - impulse response function
    # lags - lags that will be computed, in units of baud

    [nbaud, scancount] = codeset.shape
    baud=np.array(baud)

    cwta=[]
    for ii in range(scancount):
        env=baudexpand(codeset[:,ii], baud)
        
        # introduce the errors
        env2=env.copy()
        III=np.where(codeset[:,ii]==-1)[0]
#       print III
        for aa in range(len(III)):
            II=III[aa]
            env2[0,II*baud+2./3.*baud:(II+1)*baud]=1
#       pylab.plot(np.squeeze(env2),'-r')
#       pylab.hold(1)
#       pylab.plot(np.squeeze(env),'-k')
#       xxx
        env=np.squeeze(np.fliplr(env2))
        
        cwta.append(np.zeros((h.shape[0],env.shape[0]+h.shape[0]),dtype='float64'))
    
        for tau in range(h.shape[0]):
            cwta[ii][tau,tau:tau+env.shape[0]] = h[tau]*env

    wtt=np.zeros((baud*len(lags)+2*h.shape[0],env.shape[0]+h.shape[0]),dtype='float64')
    wttall=np.zeros((len(lags),baud*len(lags)+2*h.shape[0],env.shape[0]+h.shape[0],nbaud),dtype='float32')

    offs = np.arange(2*h.shape[0]-1)

    for iii in range(len(lags)):
        lgind=lags[iii]
        cwtt=[]
        for ii in range(scancount):
            lag = baud*lgind;
            print ('Lag-%d' % (lag))
            if lgind==0:
                tmp=np.transpose(conv2d(np.transpose(cwta[ii][:,lag:]), np.transpose(cwta[ii][:,0:])))
                wtt[offs,:]+=tmp
                wttall[iii,offs,:,0] = wtt[offs,:]+tmp
            else:
                tmp=np.transpose(conv2d(np.transpose(cwta[ii][:,lag:]), np.transpose(cwta[ii][:,0:-lag])))
            cwtt.append(tmp)

        if lgind != 0:
            # decode alt. codes
            for jj in range(nbaud-lgind):
                decode = codeset[jj,:]*codeset[jj+lgind,:]
                wtmp = np.zeros(cwtt[0].shape)
                for ii in range(scancount):
                    wtmp+=decode[ii]*cwtt[ii]

                # wtt contains the maximum value
                wttall[iii,lag+offs,lag:,jj]=wtmp.copy()
                tmp2=wtt[lag+offs,lag:]
                I=np.where(wtmp>tmp2)
                tmp2[I]=wtmp[I]
                wtt[lag+offs,lag:] = tmp2           

    tau=np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r=np.arange(wtt.shape[1],dtype='float64')
    
    wtt=wtt/(nbaud*scancount*baud)
    wttall=wttall/(nbaud*scancount*baud)
    
    # lag ambiguity function
    wlag=np.sum(np.sum(wttall,axis=3),axis=2)
    wlagsum=np.sum(wlag,axis=1)
        
    tmp=np.zeros((wttall.shape[0],wttall.shape[1],wttall.shape[2]+baud*nbaud),dtype='float32')
    for ii in range(wttall.shape[3]):
        tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],ii*baud),dtype='float32'),wttall[:,:,:,ii],np.zeros((wttall.shape[0],wttall.shape[1],baud*nbaud-ii*baud),dtype='float32')),axis=2)  
#       if ii==0:
#           tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*nbaud),dtype='float32'),wttall[:,:,:,ii]),axis=2)
#       else:
#           tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),dtype='float32'),wttall[:,:,:-baud*ii,ii]),axis=2)
    wttall=tmp
    
    wtt=np.concatenate((wtt,np.zeros((wtt.shape[0],baud*nbaud),dtype='float32')),axis=1)
    r=np.arange(wtt.shape[1],dtype='float64')

#   tmp=np.zeros((wttall.shape[0],wttall.shape[1],wttall.shape[2]),dtype='float32')
#   for ii in range(wttall.shape[3]):
#       if ii==0:
#           tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),dtype='float32'),wttall[:,:,:,ii]),axis=2)
#       else:
#           tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),dtype='float32'),wttall[:,:,:-baud*ii,ii]),axis=2)
#   wttall=tmp

    # range ambiguity function
    wrng=np.sum(wttall,axis=1)
    wrngsum=np.sum(wrng,axis=1)
    
    comp_all=1
    
    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all

# alt_code_fraclag
def alt_code_fraclag_old(codeset, baud, fraction, h, lags):
    """
    computes an alternating code range-lag ambiguity function

    codeset - sign sequence for the alternating code (should be nbaud x nscans size)
    baud - baud length (microseconds)
    fraction - fractional lag
    h - impulse response function
    lags - lags that will be computed, in units of baud
    """

    [nbaud, scancount] = codeset.shape

    newcodeset=np.zeros((codeset.shape[0]*fraction,codeset.shape[1]))
    for aa in range(codeset.shape[0]):
        newcodeset[aa*fraction:aa*fraction+fraction,:]=codeset[aa,:]

    nx=2*fraction-1
    decodeset=np.zeros(((nbaud-1)*nx,(nbaud-1)*nx+1,newcodeset.shape[1]))
    lagmat=np.zeros((nbaud-1)*nx+1)
    for aa in range(decodeset.shape[2]):
        for bb in range(1,codeset.shape[0]):
            sc=codeset[bb:,aa]*codeset[:-bb,aa]
            sc=np.repeat(sc[:,np.newaxis],nx,axis=1)
            decodeset[range(0,newcodeset.shape[0]-bb*fraction,fraction),1+(bb-1)*nx:1+(bb)*nx,aa]=sc
            lagmat[1+(bb-1)*nx:1+(bb)*nx]=np.arange((bb-1)*fraction+(fraction-(fraction-1)),(bb-1)*fraction+(fraction+(fraction+1))-1)
    
    lags=lagmat
    
    codeset=newcodeset
    baud=baud/fraction
    
    [nbaud, scancount] = codeset.shape
    baud=np.array(baud)
    
    cwta=[]
    for ii in range(scancount):
        env=np.squeeze(np.fliplr(baudexpand(codeset[:,ii], baud)))
        cwta.append(np.zeros((h.shape[0],env.shape[0]+h.shape[0]),dtype='float64'))
    
        for tau in range(h.shape[0]):
            cwta[ii][tau,tau:tau+env.shape[0]] = h[tau]*env

    wtt=np.zeros((baud*len(lags)+2*h.shape[0],env.shape[0]+h.shape[0]),dtype='float64')
    wlag=np.zeros((len(lags),baud*len(lags)+2*h.shape[0]),dtype='float32')
    wrng=np.zeros((len(lags),env.shape[0]+h.shape[0]),dtype='float32')    
    try:
        wttall=np.zeros((len(lags),baud*len(lags)+2*h.shape[0],env.shape[0]+h.shape[0],nbaud),dtype='float32')
        comp_all=1
    except:
        print ('Unable to allocate memory for full, 3D ambiguity function')
        comp_all=0
        wttall=[]
    comp_all=0

    offs = np.arange(2*h.shape[0]-1)

    for iii in range(lagmat.shape[0]):
        lgind=int(lagmat[iii])
        ind2do=iii
        cwtt=[]
        for ii in range(scancount):
            lag = baud*lgind;
            print ('Lag-%d' % (lag))
            if lgind==0:
                tmp=np.transpose(conv2d(np.transpose(cwta[ii][:,lag:]), np.transpose(cwta[ii][:,0:])))
                wtt[offs,:]+=tmp
                if comp_all:
                    wttall[ind2do,offs,:,0] = wtt[offs,:]+tmp
                else:
                    wlag[ind2do,lag+offs]=wlag[ind2do,lag+offs]+np.sum(tmp,axis=1)
                    wrng[ind2do,lag:]=wrng[ind2do,lag:]+np.sum(tmp,axis=0)                    
            else:
                tmp=np.transpose(conv2d(np.transpose(cwta[ii][:,lag:]), np.transpose(cwta[ii][:,0:-lag])))
            cwtt.append(tmp)

        if lgind != 0:
            # decode alt. codes
            for jj in range(nbaud-lgind):
            
                decode=decodeset[jj,iii,:]
                            
                wtmp = np.zeros(cwtt[0].shape)
                for ii in range(scancount):
                    wtmp+=decode[ii]*cwtt[ii]

                # wtt contains the maximum value
                if comp_all:
                    wttall[ind2do,lag+offs,lag:,jj]=wttall[ind2do,lag+offs,lag:,jj]+wtmp.copy()
                else:
                    wlag[ind2do,lag+offs]=wlag[ind2do,lag+offs]+np.sum(wtmp,axis=1)
                    wrng[ind2do,lag:]=wrng[ind2do,lag:]+np.sum(wtmp,axis=0)
                tmp2=wtt[lag+offs,lag:]
                I=np.where(wtmp>tmp2)
                tmp2[I]=wtmp[I]
                wtt[lag+offs,lag:] = tmp2           
#           raw_input("Enter something: ")

    tau=np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r=np.arange(wtt.shape[1],dtype='float64')
    
    wtt=wtt/(nbaud*scancount*baud)
        
    if comp_all:
        wttall=wttall/(nbaud*scancount*baud)
        tmp=np.zeros((wttall.shape[0],wttall.shape[1],wttall.shape[2]),dtype='float32')
        for ii in range(wttall.shape[3]):
            if ii==0:
                tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),dtype='float32'),wttall[:,:,:,ii]),axis=2)
            else:
                tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),dtype='float32'),wttall[:,:,:-baud*ii,ii]),axis=2)
        wttall=tmp
    
    # range and lag ambiguity function
    if comp_all:
        wlag=np.sum(np.sum(wttall,axis=3),axis=2)
        wrng=np.sum(wttall,axis=1)
    else:
        wlag=wlag/(nbaud*scancount*baud)
        wrng=wrng/(nbaud*scancount*baud)
    wlagsum=np.sum(wlag,axis=1)
    wrngsum=np.sum(wrng,axis=1)
    
    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all, lags


# alt_code_fraclag
def alt_code_fraclag(codeset, baud, fraction, h, lags, ph=''):
    """
    computes an alternating code range-lag ambiguity function

    codeset - sign sequence for the alternating code (should be nbaud x nscans size)
    baud - baud length (microseconds)
    fraction - fractional lag
    h - impulse response function
    lags - lags that will be computed, in units of baud
    ph - file with phase corrections
    """

    [nbaud, scancount] = codeset.shape

    newcodeset = np.zeros((codeset.shape[0] * fraction, codeset.shape[1]))
    for aa in range(codeset.shape[0]):
        newcodeset[aa * fraction:aa * fraction + fraction,:] = codeset[aa,:]

    nx = 2*fraction-1
    decodeset = np.zeros(((nbaud-1)*nx,(nbaud-1)*nx+1,newcodeset.shape[1]))
    lagmat = np.zeros((nbaud-1) * nx + 1)
    for aa in range(decodeset.shape[2]):
        for bb in range(1,codeset.shape[0]):
            sc = codeset[bb:,aa] * codeset[:-bb,aa]
            sc = np.repeat(sc[:,np.newaxis], nx, axis=1)
            decodeset[range(0,newcodeset.shape[0] - bb * fraction, fraction),
                           1 + (bb - 1) * nx : 1 + (bb) * nx, aa] = sc
            lagmat[1+(bb-1)*nx:1+(bb)*nx] = np.arange((bb-1) * fraction
                    + (fraction - (fraction - 1)), (bb-1) * fraction
                    + (fraction + (fraction + 1)) - 1)

    codeset = newcodeset
    baud = baud // fraction

    [nbaud, scancount] = codeset.shape
    baud = np.array(baud)

    # pulse envelope
    if ph=='':
        envsc = np.ones((nbaud*baud),dtype='float32')
    else:
        ph = np.array(np.loadtxt(ph))
        envsc = np.ones((nbaud*baud), dtype='complex64')
        envsc = np.absolute(envsc) * np.exp(1.0j*ph)

    cwta=[]
    for ii in range(scancount):
        env = envsc * np.squeeze(np.fliplr(baudexpand(codeset[:,ii], baud)))
        cwta.append(np.zeros((h.shape[0], env.shape[0] + h.shape[0]), dtype=env.dtype))

        for tau in range(h.shape[0]):
            cwta[ii][tau,tau:tau+env.shape[0]] = h[tau]*env

    wtt = np.zeros(( baud * len(lags) + 2 * h.shape[0], env.shape[0] + h.shape[0]),
            dtype=env.dtype)
    wlag = np.zeros((len(lags), baud*len(lags)+2*h.shape[0]),dtype=env.dtype)
    wrng = np.zeros((len(lags),env.shape[0]+h.shape[0]),dtype=env.dtype)
    try:
        forgetit   # this is to break the code
        wttall= np.zeros((len(lags),baud*len(lags)+2*h.shape[0],
            env.shape[0]+h.shape[0],nbaud),dtype=env.dtype)
        comp_all=1
    except:
        print ('Unable to allocate memory for full, 3D ambiguity function')
        comp_all=0
        wttall=[]	

    offs = np.arange(2*h.shape[0]-1)

    for iii in range(lagmat.shape[0]):
        print(iii, "out of", lagmat.shape[0])
        lgind=int(lagmat[iii])
        cwtt=[]
        for ii in range(scancount):
            lag = baud*lgind;
            if ii % 10 ==0:
                print ('Lag-%d' % (lag), ii, "out of scancount:", scancount)
            if lgind==0:
                tmp = np.transpose(conv2d(np.transpose(cwta[ii][:,lag:]),
                    np.transpose(cwta[ii][:,0:])))
                wtt[offs,:]+=tmp
                if comp_all:
                    wttall[lgind,offs,:,0] = wtt[offs,:]+tmp
                else:
                    wlag[lgind,lag+offs] = wlag[lgind, lag+offs] + np.sum(tmp,axis=1)
                    wrng[lgind,lag:] = wrng[lgind, lag:] + np.sum(tmp, axis=0)
            else:
                tmp = np.transpose(conv2d(np.transpose(cwta[ii][:,lag:]),
                    np.transpose(cwta[ii][:,0:-lag])))
            cwtt.append(tmp)

        if lgind != 0:
            # decode alt. codes
            for jj in range(nbaud-lgind):

                decode=decodeset[jj,iii,:]

                wtmp = np.zeros(cwtt[0].shape, dtype=env.dtype)
                for ii in range(scancount):
                    wtmp+=decode[ii]*cwtt[ii]

                # wtt contains the maximum value
                if comp_all:
                    wttall[lgind,lag+offs,lag:,jj]=wttall[lgind,lag+offs,lag:,jj]+wtmp.copy()
                else:
                    wlag[lgind,lag+offs] = wlag[lgind,lag+offs]+np.sum(wtmp,axis=1)
                    wrng[lgind,lag:] = wrng[lgind,lag:]+np.sum(wtmp,axis=0)
                tmp2=wtt[lag+offs,lag:]
                I = np.where(wtmp > tmp2)
                tmp2[I] = wtmp[I]
                wtt[lag+offs,lag:] = tmp2
#           raw_input("Enter something: ")

    tau = np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r = np.arange(wtt.shape[1],dtype='float64')

    wtt = wtt / (nbaud*scancount*baud)

    if comp_all:
        wttall = wttall // (nbaud*scancount*baud)
        tmp = np.zeros((wttall.shape[0],wttall.shape[1],wttall.shape[2]),dtype=env.dtype)
        for ii in range(wttall.shape[3]):
            if ii==0:
                tmp += np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),
                    dtype=env.dtype),wttall[:,:,:,ii]),axis=2)
            else:
                tmp += np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),
                    dtype=env.dtype),wttall[:,:,:-baud*ii,ii]),axis=2)
        wttall = tmp
    # range and lag ambiguity function
    if comp_all:
        wlag = np.sum(np.sum(wttall, axis=3), axis=2)
        wrng = np.sum(wttall, axis=1)
    else:
        wlag = wlag / (nbaud * scancount * baud)
        wrng = wrng / (nbaud * scancount * baud)
    wlagsum = np.sum(wlag, axis=1)
    wrngsum = np.sum(wrng, axis=1)
    
    lags = np.arange(lagmat.max()+1)

    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all, lags


# alt_code
def alt_code(codeset, baud, h, lags, ph=''):
    #
    # computes an alternating code range-lag ambiguity function
    #
    # codeset - sign sequence for the alternating code (should be nbaud x nscans size)
    # baud - baud length (microseconds)
    # h - impulse response function
    # lags - lags that will be computed, in units of baud

    [nbaud, scancount] = codeset.shape
    baud=np.array(baud)

    # pulse envelope
    if ph=='':
        envsc=np.ones((nbaud*baud),dtype='float32')
    else:
        ph=np.array(np.loadtxt(ph))
        envsc = np.ones((nbaud*baud),dtype='complex64')
        envsc = np.absolute(envsc)*np.exp(1.0j*ph)

    cwta=[]
    for ii in range(scancount):
        env=envsc*np.squeeze(np.fliplr(baudexpand(codeset[:,ii], baud)))
        cwta.append(np.zeros((h.shape[0],env.shape[0]+h.shape[0]),dtype=env.dtype))

        for tau in range(h.shape[0]):
            cwta[ii][tau,tau:tau+env.shape[0]] = h[tau]*env

    wtt = np.zeros((baud*len(lags)+2*h.shape[0],env.shape[0]+h.shape[0]),
            dtype=env.dtype)
    wttall = np.zeros((len(lags),baud*len(lags)+2*h.shape[0],env.shape[0]+h.shape[0],nbaud),
            dtype=env.dtype)

    offs = np.arange(2*h.shape[0]-1)

    for iii in range(len(lags)):
        lgind=lags[iii]
        cwtt=[]
        for ii in range(scancount):
            lag = baud*lgind;
            print ('Lag-%d' % (lag))
            if lgind==0:
                tmp=np.transpose(conv2d(np.transpose(cwta[ii][:,lag:]), np.transpose(np.conj(cwta[ii][:,0:]))))
                wtt[offs,:]+=tmp
                wttall[iii,offs,:,0] = wtt[offs,:]+tmp
            else:
                tmp=np.transpose(conv2d(np.transpose(cwta[ii][:,lag:]), np.transpose(np.conj(cwta[ii][:,0:-lag]))))
            cwtt.append(tmp)

        if lgind != 0:
            # decode alt. codes
            for jj in range(nbaud-lgind):
                decode = codeset[jj,:]*codeset[jj+lgind,:]
                wtmp = np.zeros(cwtt[0].shape)
                for ii in range(scancount):
                    wtmp+=decode[ii]*cwtt[ii]

                # wtt contains the maximum value
                wttall[iii,lag+offs,lag:,jj]=wtmp.copy()
                tmp2=wtt[lag+offs,lag:]
                I=np.where(np.absolute(wtmp)>np.absolute(tmp2))
                tmp2[I]=wtmp[I]
                wtt[lag+offs,lag:] = tmp2           

    tau=np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r=np.arange(wtt.shape[1],dtype='float64')

    wtt=wtt/(nbaud*scancount*baud)
    wttall=wttall/(nbaud*scancount*baud)

    # lag ambiguity function
    wlag=np.sum(np.sum(wttall,axis=3),axis=2)
    wlagsum=np.sum(wlag,axis=1)

    tmp=np.zeros((wttall.shape[0],wttall.shape[1],wttall.shape[2]),dtype=env.dtype)
    for ii in range(wttall.shape[3]):
        if ii==0:
            tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),dtype=env.dtype),wttall[:,:,:,ii]),axis=2)
        else:
            tmp+=np.concatenate((np.zeros((wttall.shape[0],wttall.shape[1],baud*ii),dtype=env.dtype),wttall[:,:,:-baud*ii,ii]),axis=2)
    wttall=tmp
    
    # range ambiguity function
    wrng=np.sum(wttall,axis=1)
    wrngsum=np.sum(wrng,axis=1)
    
    comp_all=1

    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all

# long_pulse
def long_pulse(nbaud, baud, h, lags, ph=''):
    #
    # computes a long pulse range-lag ambiguity function
    #
    # nbaud - number of bauds
    # baud - baud length (microseconds)
    # h - impulse response function
    # lags - lags that will be computed, in units of baud

    # pulse envelope
    if ph=='':
        env=np.ones((nbaud*baud),dtype='float32')
    else:
        ph=np.array(np.loadtxt(ph))
        env=np.ones((nbaud*baud),dtype='complex64')
        #env = np.absolute(env)*np.exp(1.0j*np.arange(480)*0.001)
        env = np.absolute(env)*np.exp(1.0j*ph)
    
    wta=np.zeros((h.shape[0],env.shape[0]+h.shape[0]),dtype=env.dtype)

    for tau in range(h.shape[0]):
        wta[tau,tau:tau+env.shape[0]]=h[tau]*env
    
    wtt=np.zeros((baud*len(lags)+2*h.shape[0],env.shape[0]+h.shape[0]),dtype=env.dtype)
    wlag=np.zeros((len(lags),baud*len(lags)+2*h.shape[0]),dtype=env.dtype)
    wrng=np.zeros((len(lags),env.shape[0]+h.shape[0]),dtype=env.dtype)
    try:
        wttall=np.zeros((len(lags),baud*len(lags)+2*h.shape[0],env.shape[0]+h.shape[0]),dtype=env.dtype)
        comp_all=1
    except:
        print ('Unable to allocate memory for full, 3D ambiguity function')
        comp_all=0
        wttall=[]

    offs = np.arange(2*h.shape[0]-1)

    for ii in range(len(lags)):
        lgind=lags[ii]
        lag = baud*lgind
        print ('Lag-%d' % (lag))
        if lgind==0:
            tmp=np.transpose(conv2d(np.transpose(wta[:,lag:]), np.transpose(np.conj(wta[:,0:]))))
        else:
            tmp=np.transpose(conv2d(np.transpose(wta[:,lag:]), np.transpose(np.conj(wta[:,0:-lag]))))
        if comp_all:
            wttall[ii,lag+offs,lag:]=tmp
        else:
            wlag[ii,lag+offs]=wlag[ii,lag+offs]+np.sum(tmp,axis=1)
            wrng[ii,lag:]=wrng[ii,lag:]+np.sum(tmp,axis=0)
        tmp2=wtt[lag+offs,lag:]
        I=np.where(np.absolute(tmp)>np.absolute(tmp2))
        tmp2[I]=tmp[I]
        wtt[lag+offs,lag:]=tmp2

    tau=np.arange(wtt.shape[0],dtype='float64')-h.shape[0]+1
    r=np.arange(wtt.shape[1],dtype='float64')
    
    wtt=wtt/(nbaud*baud)
    if comp_all:
        wttall=wttall/(nbaud*baud)
        wlag=np.sum(wttall,axis=2)
        wrng=np.sum(wttall,axis=1)
    else:
        wlag=wlag/(nbaud*baud)
        wrng=wrng/(nbaud*baud)
    
    wlagsum=np.sum(wlag,axis=1)
    wrngsum=np.sum(wrng,axis=1)
    
    
    return tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all


# plot_lamb
def plot_lamb(x,y,wtt,wttall,wlag,wrng,plotall=1,plotson=0,plotdir='',lags=[]):

    if plotson:
        if not os.path.exists(plotdir):
            try: 
                os.mkdir(plotdir)
            except:
                plotson=0

    labsize=12
    figsize=(12, 9)
 
    #print wlag

    iscomplex = np.iscomplex(wlag).any()
    
    #xxx

    fig = plt.figure(num=1,figsize=figsize)
    if not iscomplex:
        ax = fig.add_axes(111)
        ax.plot(x,np.transpose(wlag))
        ax.set_ylabel('Scaling', fontsize=labsize)
        ax.set_xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
        ax.set_title('Lag Ambiguity Function', fontsize=labsize, horizontalalignment='center')
    else:
        ax = fig.add_axes(221)
        ax.plot(x,np.transpose(wlag.real))
        ax.set_ylabel("Real", fontsize=labsize)
        ax.set_title('Lag Ambiguity Function', fontsize=labsize, horizontalalignment='center')
        ax = fig.add_axes(222)
        ax.plot(x,np.transpose(wlag.imag))
        ax.set_title('Lag Ambiguity Function', fontsize=labsize, horizontalalignment='center')
        ax.set_ylabel("Imag", fontsize=labsize)
        ax = fig.add_axes(223)
        ax.plot(x,np.transpose(np.absolute(wlag)))
        ax.set_xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
        ax.set_ylabel("Magitude", fontsize=labsize)
        ax = fig.add_axes(224)
        ax.plot(x,np.transpose((np.arctan2(wlag.imag,wlag.real))))
        ax.set_xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
        ax.set_ylabel("Phase", fontsize=labsize)
    if plotson:
        oname='wlag.png'
        fullname =  os.path.join(plotdir,oname)
        print(f"Saving fig: {fullname}")
        fig.savefig(fullname)

    fig = plt.figure(num=2, figsize=figsize)
    if not iscomplex:
        ax = fig.add_axes(111)
        ax.plot(y,np.transpose(wrng))
        ax.set_ylabel('Scaling', fontsize=labsize)
        ax.set_xlabel(r"$Range \ (\mu s)$", fontsize=labsize)
        ax.set_title('Range Ambiguity Function', fontsize=labsize, horizontalalignment='center')
    else:
        ax = fig.add_axes(221)
        ax.plot(y,np.transpose(wrng.real))
        ax.set_ylabel("Real", fontsize=labsize)
        ax.set_title('Range Ambiguity Function', fontsize=labsize, horizontalalignment='center')
        ax = fig.add_axes(222)
        ax.plot(y,np.transpose(wrng.imag))
        ax.set_title('Range Ambiguity Function', fontsize=labsize, horizontalalignment='center')
        ax.set_ylabel("Imag", fontsize=labsize)
        ax = fig.add_axes(223)
        ax.plot(y,np.transpose(np.absolute(wrng)))
        ax.set_xlabel(r"$Range \ (\mu s)$", fontsize=labsize)
        ax.set_ylabel("Magitude", fontsize=labsize)
        ax = fig.add_axes(224)
        ax.plot(y,np.transpose((np.arctan2(wrng.imag,wrng.real))))
        ax.set_xlabel(r"$Range \ (\mu s)$", fontsize=labsize)
        ax.set_ylabel("Phase", fontsize=labsize)
    if plotson:
        oname='wrng.png'
        fig.savefig(os.path.join(plotdir,oname))
    
    fig = plt.figure(num=3, figsize=figsize)
    if not iscomplex:
        ax = fig.add_axes(111)
        pcm = ax.imshow(np.transpose(wtt),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()])
        ax.set_ylabel(r"$Range \ (\mu s)$", fontsize=labsize)
        ax.set_xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
        ax.set_title('Full 2d Ambiguity Function', fontsize=labsize, horizontalalignment='center')
        fig.colorbar(pcm)
    else:
        ax = fig.add_axes(221)
        pcm = ax.imshow(np.transpose(wtt.real),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()],aspect='auto')
        fig.colorbar(pcm)
        ax.set_title('Real', fontsize=labsize, horizontalalignment='center')
        ax.set_ylabel(r"$Range \ (\mu s)$", fontsize=labsize)                
        ax = fig.add_axes(222)
        pcm = ax.imshow(np.transpose(wtt.imag),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()],aspect='auto')
        fig.colorbar(pcm)
        ax.set_title('Imag', fontsize=labsize, horizontalalignment='center')
        ax = fig.add_axes(223)
        pcm = ax.imshow(np.transpose(np.absolute(wtt)),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()],aspect='auto')
        fig.colorbar(pcm)
        ax.set_title('Magnitude', fontsize=labsize, horizontalalignment='center')        
        ax.set_xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
        ax.set_ylabel(r"$Range \ (\mu s)$", fontsize=labsize)        
        ax = fig.add_axes(224)
        pcm = ax.imshow(np.transpose(np.arctan2(wtt.imag,wtt.real)),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()],aspect='auto')
        fig.colorbar(pcm)
        ax.set_title('Phase', fontsize=labsize, horizontalalignment='center')  
        ax.set_xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
              
    if plotson:
        oname='wtt.png'
        fig.savefig(os.path.join(plotdir,oname))
    
    #pylab.show()
    
    if plotall:
        for aa in range(len(lags)):
            pylab.close('all')
            figg=pylab.figure(figsize=figsize)
            if len(lags)==0:
                lag=aa
            else:
                lag=lags[aa]
            if not iscomplex:
                pylab.imshow(np.transpose(wttall[aa,:,:]),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()])
                ax=pylab.gca()
                ax.set_ylabel(r"$Range \ (\mu s)$", fontsize=labsize)
                ax.set_xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
                ax.set_title('2d Ambiguity Function, Lag %d' % (lag), fontsize=labsize, horizontalalignment='center')
                pylab.colorbar()
            else:
                pylab.subplot(221)
                pylab.imshow(np.transpose(wttall[aa,:,:].real),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()],aspect='auto')
                pylab.colorbar()
                pylab.title('Real, Lag %d' % (lag), fontsize=labsize, horizontalalignment='center')
                pylab.ylabel(r"$Range \ (\mu s)$", fontsize=labsize)                    
                pylab.subplot(222)
                pylab.imshow(np.transpose(wttall[aa,:,:].imag),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()],aspect='auto')
                pylab.colorbar()
                pylab.title('Imag, Lag %d' % (lag), fontsize=labsize, horizontalalignment='center')
                pylab.subplot(223)
                pylab.imshow(np.transpose(np.absolute(wttall[aa,:,:])),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()],aspect='auto')
                pylab.colorbar()
                pylab.xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
                pylab.ylabel(r"$Range \ (\mu s)$", fontsize=labsize)                    
                pylab.title('Magnitude, Lag %d' % (lag), fontsize=labsize, horizontalalignment='center')
                pylab.subplot(224)        
                pylab.imshow(np.transpose(np.arctan2(wttall[aa,:,:].imag,wttall[aa,:,:].real)),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()],aspect='auto')
                pylab.colorbar() 
                pylab.title('Phase, Lag %d' % (lag), fontsize=labsize, horizontalalignment='center')
                pylab.xlabel(r"$\tau \ (\mu s)$", fontsize=labsize)
            if plotson:
                oname='wtt-lag%d.png' % (lag)
                figg.savefig(os.path.join(plotdir,oname))
            pylab.show()
    
    return
    
def write_outputfile(fhandle,dict2do,keys2do=[],groupname='',name=''):

    if groupname == '':
        group=fhandle.root
    else:
        if fhandle.groups.__contains__('/'+groupname):
            group='/'+groupname
        else:
            group=fhandle.create_group(fhandle.root, groupname, 'Dataset')

    if len(keys2do)==0:
        try:
            fhandle.remove_node(group,name)
        except:
            ''
        fhandle.create_array(group,name, dict2do, "Dataset")
    else:
        for key in keys2do:
            fhandle.create_array(group, key, dict2do[key], "Dataset")
    
    return      
    
# compute_lamb
def compute_lamb(ltype,hfile,MCIC,in1,in2,outdir,lags=[]):
    """
    desc: Computes the impulse response function

    Input:
    ltype - pulse type (0=Long Pulse, 1=Alternating Code,
            2=Alternating Code with phase errs, 3=Barker Code,
            4=Fractional Lag Alternating Code)
    hfile - file name containing the impulse response function
           (time spacing should be one microsecond)
    lags - lags that will be computed, in units of baud (default is all)
    MCIC - length 2 list containing the decimation factors

    ltype-0:
      in1 - number of samples per pulse
      in2 - [sample time (microseconds), OPTIONAL phase file]
    ltype-1:
      in1 - codeset
      in2 - [baud length (microseconds)]
    ltype-2:
      in1 - codeset
      in2 - [baud length (microseconds)]
    ltype-3:
      in1 - code length
      in2 - [baud length (microseconds) oversampling factor]
    ltype-4:
      in1 - code length
      in2 - [baud length(microseconds) fraction]


   """

    with open(hfile, 'r') as fp:
        h = np.loadtxt(fp)

    if len(MCIC)!=2:
        print ('MCIC must be length 2')
        sys.exit(-1)
    else:
        dt = float(MCIC[0]) * float(MCIC[1]) / 50.0

    out1=[]
    # Long Pulse
    if ltype == 0:
        name1 = 'nSamples'
        name2 = ['SampleTime', 'PhaseFile']
        typedesc = 'Long Pulse'
        if len(lags) == 0:
            lags = list(range(in1))
        if len(in2) == 1:
            phFile = ''
        else: 
            phFile = in2[1]
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all \
                = long_pulse(in1, int(float(in2[0])/dt), h, lags, ph=phFile)
        PulseLength = float(in1) * float(in2[0])
        LagSp = int(float(in2[0]) / dt)
    # Alternating Code
    elif ltype==1:
        name1 = 'Codeset'
        name2 = ['BaudLength']
        typedesc='Alternating Code'
        if len(lags)==0:
            lags = list(range(in1.shape[0]))
        if len(in2) == 1:
            phFile = ''
        else:
            phFile = in2[1]
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum,comp_all \
                = alt_code(in1, int(float(in2[0])/dt), h, lags, ph=phFile)
        PulseLength = float(in1.shape[0]) * float(in2[0])
        LagSp = int(float(in2[0]) / dt)
    # Alternating Code with Phase Errors
    elif ltype==2:
        name1 = 'Codeset'
        name2 = ['BaudLength']
        typedesc = 'Alternating Code with Coding Error'
        if len(lags) == 0:
            lags = list(range(in1.shape[0]))
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all \
                = alt_code_phaserr(in1, int(float(in2[0])/dt), h, lags)
        PulseLength=float(in1.shape[0])*float(in2[0])
        LagSp=int(float(in2[0])/dt)
    # Barker Code
    elif ltype==3:
        name1 = 'Length'
        if len(in2) == 2:
            name2 = ['BaudLength', 'OversamplingFactor']
            ctype=0
        elif len(in2) == 3:
            name2 = ['BaudLength', 'OversamplingFactor', 'CodeType']
            ctype = in2[2]
        typedesc='Barker Code'
        lags = [0]
        Fs = 0.5 / dt*1.0e6
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all, BWfilt, BW \
                = barker_code(in1, int(float(in2[0])/dt), int(in2[1]), h, Fs, btype=ctype)
        PulseLength=float(in1)*float(in2[0])
        LagSp=int(float(in2[0])/dt)
        out1=[BWfilt,BW]
        nameOut1=['FilterBandwidth', 'Bandwidth']
    # Fractional Lag Alternating Code
    elif ltype == 4:
        name1 = 'Codeset'
        name2 = ['BaudLength','Fraction']
        typedesc = 'Fractional Lag Alternating Code'
        if len(lags) == 0:
            lags = list(range(in1.shape[0] * in2[1]))
        if len(in2) == 2:
            phFile = ''
        elif len(in2) == 3:
            phFile = in2[2]
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum,comp_all,lags \
                = alt_code_fraclag(in1, int(float(in2[0])/dt), in2[1], h, lags, ph=phFile)
        PulseLength=float(in1.shape[0])*float(in2[0])
        LagSp=int(float(in2[0])/in2[1])
    # Fractional Lag Alternating Code
    elif ltype==44:
        name1 = 'Codeset'
        name2 = ['BaudLength','Fraction']
        typedesc = 'Fractional Lag Alternating Code'
        if len(lags) == 0:
            lags = list(range(in1.shape[0]*in2[1]))
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum,comp_all,lags \
                = alt_code_fraclag_old(in1, int(float(in2[0])/dt), in2[1], h, lags)
        PulseLength = float(in1.shape[0]) * float(in2[0])
        LagSp = int(float(in2[0]) / dt/in2[1])
        print (LagSp)
    # Coherent Code
    elif ltype==5:    
        name1='Codeset'
        if len(in2)==2:
            name2=['BaudLength', 'OversamplingFactor']
            ctype=0
        elif len(in2)==3:
            name2=['BaudLength', 'OversamplingFactor', 'CodeType']
            ctype=in2[2]
        typedesc='Barker Code'
        lags=0
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all, htotal \
                = coh_code(in1, int(float(in2[0])/dt), int(in2[1]), h)
        PulseLength=float(in1.shape[0])*float(in2[0])
        LagSp=int(float(in2[0])/dt)
    # Multi pulse
    elif ltype==6:
        name1='Codeset'
        name2=['Baudlength', 'OversamplingFactor', 'Pattern']
        typedesc='MultiPulse'
        if len(lags)==0:
            lags = list(range(in1.shape[0]))
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all \
                = mpulse(in1, in2[2], int(float(in2[0])/dt), int(in2[1]), h, lags)
        PulseLength=float(in1.shape[0])*float(in2[0])
        LagSp=int(float(in2[0]))
        print ("LagSp",LagSp)
    # Double pulse
    elif ltype==7:
        name1='Codeset'
        name2=['Baudlength', 'OversamplingFactor', 'CodeLength']
        typedesc='MultiPulse'
        if len(lags)==0:
            lags=list(range(in1.shape[0]))
        tau, r, wttall, wtt, wlag, wlagsum, wrng, wrngsum, comp_all \
                = dpulse(in1, in2[2], int(float(in2[0])/dt), int(in2[1]), h, lags)
        PulseLength=float(in1.shape[0])*float(in2[0])
        LagSp=int(float(in2[0]))
        print ("LagSp",LagSp)
    
    in2[0]=in2[0]*1e-6

    tau=tau*dt
    r=r*dt
          
    saveOutput=1
    # write output files
    if saveOutput:
        try:
            outname=os.path.join(outdir,'AmbFuncFull.h5')
            with tables.open_file(outname, mode = "w",
                    title = "Ambiguity Function File") as fp:
                write_outputfile(fp,tau.astype('float64')*1e-6,groupname='',name='Delay')
                write_outputfile(fp,r.astype('float64')*v_lightspeed/2.0*1e-6,groupname='',name='Range')
                if comp_all:
                    write_outputfile(fp,wttall,groupname='',name='WttAll')
                write_outputfile(fp,wtt,groupname='',name='Wtt')
                write_outputfile(fp,wlag,groupname='',name='Wlag')
                write_outputfile(fp,wlagsum,groupname='',name='WlagSum') 
                write_outputfile(fp,wrng,groupname='',name='Wrange')
                write_outputfile(fp,wrngsum,groupname='',name='WrangeSum')   
                write_outputfile(fp,in1,groupname='',name=name1)
                for aa in range(len(in2)):
                    write_outputfile(fp,in2[aa],groupname='',name=name2[aa])
                for aa in range(len(out1)):
                    write_outputfile(fp,out1[aa],groupname='',name=nameOut1[aa])
                write_outputfile(fp,ltype,groupname='',name='Type')
                write_outputfile(fp,typedesc,groupname='',name='TypeDesc')
                write_outputfile(fp,hfile,groupname='',name='FilterFile')
                write_outputfile(fp,np.asarray(MCIC),groupname='',name='MCIC')
                write_outputfile(fp,dt*1e-6,groupname='',name='dt')
                write_outputfile(fp,np.array(lags).astype('float64')*float(LagSp)*1e-6,groupname='',name='Lags')
                write_outputfile(fp,PulseLength*1e-6,groupname='',name='PulseLength')
        except:
            print ('unable to create ambfuncfull')

        outname=os.path.join(outdir,'AmbFunc.h5')
        with tables.open_file(outname, mode = "w",
                title = "Ambiguity Function File") as fp:

            write_outputfile(fp,tau.astype('float64')*1e-6,groupname='',name='Delay')
            write_outputfile(fp,r.astype('float64')*v_lightspeed/2.0*1e-6,groupname='',name='Range')
            write_outputfile(fp,wlag,groupname='',name='Wlag')
            write_outputfile(fp,wlagsum,groupname='',name='WlagSum') 
            write_outputfile(fp,wrng,groupname='',name='Wrange')
            write_outputfile(fp,wrngsum,groupname='',name='WrangeSum')   
        #    write_outputfile(fp,in1,groupname='',name=name1)
        #    for aa in range(len(in2)):
        #        write_outputfile(fp,in2[aa],groupname='',name=name2[aa])
            for aa in range(len(out1)):
                write_outputfile(fp,out1[aa],groupname='',name=nameOut1[aa])
        #    write_outputfile(fp,ltype,groupname='',name='Type')
        #    write_outputfile(fp,ltypedesc,groupname='',name='TypeDesc')
        #    write_outputfile(fp,hfile,groupname='',name='FilterFile')
        #    write_outputfile(fp,MCIC,groupname='',name='MCIC')
        #    write_outputfile(fp,dt*1e-6,groupname='',name='dt')
            write_outputfile(fp,np.array(lags).astype('float64')*float(LagSp)*1e-6,groupname='',name='Lags')
        #    write_outputfile(fp,PulseLength*1e-6,groupname='',name='PulseLength')

    
    # make plots    
    plot_lamb(tau,r,wtt,wttall,wlag,wrng,plotson=1,plotdir=os.path.join(outdir,'plots'),plotall=comp_all,lags=lags)
    
    return tau,r,wttall,wtt,wlag,wrng
