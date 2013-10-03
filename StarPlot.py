#!/usr/bin/env python
#This script plots the surface and contour plots of the region in image on TIRSPEC machine
import pyfits
import readline
import sys
import os
import os.path

Xsizeby2=30  #The length of square box to plot will be double this size +1.
Ysizeby2=30 #The width of square box to plot will be double this size +1.
Gnuplotfile='./StarPlot.gnuplot'
if len(sys.argv) < 2 or not os.path.isfile(sys.argv[1]):
    print("-"*30)
    print("Usage  : "+sys.argv[0]+" ImageName")
    print("Eg: "+sys.argv[0]+" /data/20130621/test-32.fits")
    print("Eg: "+sys.argv[0]+" /data/20130621/test-{32..38}.fits")
    print("-"*30)
    exit(1)

Fullimage=pyfits.getdata(sys.argv[1])
imgcount=1
print(sys.argv[1])
if len(sys.argv) > 2 : #Then Average all the input frames
    for img in sys.argv[2:] :
        if os.path.isfile(img):
            Fullimage+=pyfits.getdata(img)
            imgcount+=1
            print(img)
        else: print ('Image '+img+' not found and not used for averaging.')
    Fullimage/=float(imgcount)

while 1 :
    try:
        XYcoords=raw_input('Center of Location to plot (Eg: X Y) :').strip(' ')
    except (KeyboardInterrupt, EOFError):
        print("\n Exiting.. \n ")
        exit(0)
    Xcenter=int(eval(XYcoords.split()[1])) -1  #Python has opposite X and Y than ds9 or DV and starts at 0
    Ycenter=int(eval(XYcoords.split()[0])) -1
    outfile=open('/tmp/gridtoplot.table','w')
    for i in range(Xcenter-Xsizeby2,Xcenter+Xsizeby2 +1):
        for j in range(Ycenter-Ysizeby2,Ycenter+Ysizeby2 +1):
            outfile.write('%d %d %f \n'%(i,j,Fullimage[i,j]))
        outfile.write('\n')
    outfile.close()
    os.system('gnuplot '+Gnuplotfile)

        

