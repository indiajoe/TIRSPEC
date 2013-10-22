#!/usr/bin/python
# This script can be used to quickly generate bad pixel mask
# Use flats  ~5 times different in brightness
# Enjoy!!!---------------------- indiajoe@gmail.com

# Set the upper and lower cut of threshold to detect as bad pixel
LsigCut=10
HsigCut=10

import os
import sys
import numpy as np

if len(sys.argv) != 3 or not os.path.isfile(sys.argv[1]) or not os.path.isfile(sys.argv[2]):
    print("-"*30)
    print("Usage : "+sys.argv[0]+" BrightFlat.fits NotBrightFlat.fits")
    print("Eg: "+sys.argv[0]+" FlatHigh.fits FlatLow.fits ")
    print("-"*30)
    exit(1)


import pyfits
Createpl=False
print("Using Lower Sigma Cut %f "%(LsigCut))
print("Using Upper Sigma Cut %f "%(HsigCut))
print("Give suffix .pl only if you want to generate a binary mask file using iraf.text2mask also. \n Text file will also be generated with .txt suffix.")
Outfile=raw_input("Output mask name : ")
if Outfile[-3:] == '.pl' :
    from pyraf import iraf
    Outfname=Outfile[-3:]+".txt"
    Createpl=True
else:
    Outfname=Outfile

ImgH=sys.argv[1]
ImgL=sys.argv[2]

hdulistH=pyfits.open(ImgH)
hdulistL=pyfits.open(ImgL)

Hdata=hdulistH[0].data
Ldata=hdulistL[0].data
hdulistH.close()
hdulistL.close()
#Taking Ratio of two flats
RatioHL=Hdata/Ldata

#Replacing all inf to 0
RatioHL[RatioHL==np.inf]=0
RatioHL[RatioHL==-np.inf]=0
#Replacing all nan to 0
RatioHL=np.nan_to_num(RatioHL)

MedianR=np.median(RatioHL)
StdR=np.std(RatioHL)
MadR=np.median(np.abs(RatioHL-MedianR))
print("Median Ratio = %f , Std = %f Mad =%f "%(MedianR,StdR,MadR))

#(XLbp,YLbp)=np.where(RatioHL < MedianR-LsigCut*StdR)
#(XHbp,YHbp)=np.where(RatioHL > MedianR+HsigCut*StdR)
(XLbp,YLbp)=np.where(RatioHL < MedianR-LsigCut*MadR*1.4826)
(XHbp,YHbp)=np.where(RatioHL > MedianR+HsigCut*MadR*1.4826)


OutFILE=open(Outfname,'w')
for X,Y in zip(XLbp,YLbp): OutFILE.write(str(Y+1)+' '+str(X+1)+'\n')
for X,Y in zip(XHbp,YHbp): OutFILE.write(str(Y+1)+' '+str(X+1)+'\n')
OutFILE.close()
print("Text file of list of Bad pixels created in "+Outfname)

if Createpl :
    iraf.proto(_doprint=0)
    iraf.text2mask.unlearn()
    iraf.text2mask(text=Outfname,mask=Outfile,ncols=1024,nlines=1024)
    print("Binary .pl file of Bad pixels created in "+Outfile)
    print("To view mask in IRAF : display "+ImgH+" 1 overlay="+Outfile)



    

    
