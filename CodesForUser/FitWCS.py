#!/usr/bin/env python
""" Fit the WCS components by least square fit of a set of stars """
import numpy as np
import scipy.optimize as optimize

from astropy import wcs
from astropy.io import fits

def fitWCSobject(XYRaDecfilename,updatefits=None):
    """ Fit the WCS components by least square fit of a set of stars in input XYRaDecfilename
    Input:
         XYRaDecfilename : File name of text file containing 4 columns of : Xcoord  Ycoord Ra(degree) Dec(degree)
         updatefits (optional): Fits file into which the WCS header kewords have to be written after fitting.
    Output:
         w : the WCS astropy coordinate object, which can be used for further convertions.
    """    
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [510, 588]     # Reference pixel to enter header
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]    #Projection type

    #Read the pixel and ra dec coordinates.
    XYRaDec=np.loadtxt(XYRaDecfilename)#'XYRaDec.txt')

    def func(ip):
        w.wcs.crval = [ip[0], ip[1]]
        w.wcs.cdelt = np.array([ip[2], ip[3]])
        RaDec=w.all_pix2world(zip(XYRaDec[:,0],XYRaDec[:,1]),1)
        return XYRaDec[:,2:].flatten()-RaDec[:,0:].flatten()

    ip=[XYRaDec[0,2],XYRaDec[0,3],8.56E-5,8.56E-5]
#    print ip
    newip,ierr=optimize.leastsq(func, ip)
    print('Ierr='+str(ierr))
    if ier not in [1, 2, 3, 4]:
        msg = "Fitting WCS failed, optimal parameters not found!! \n Verify "+XYRaDecfilename+" contents are correct and rerun the code."
        raise RuntimeError(msg)

    print('Fitted coefficents in header')
    print(w.wcs.to_header())
    print('-'*20)
    print('Input coordinates :')
    print(XYRaDec)
    print('Fitted coordinates :')
    print(w.all_pix2world(zip(XYRaDec[:,0],XYRaDec[:,1]),1))
    print('-'*20)

    if updatefits is not None:   #We have to add the WCS coefficents to fits header
        print('Updating fits header of {0} with new WCS'.format(updatefits))
        hdulist=fits.open(updatefits,mode='update')
        hdulist[0].header.extend(w.to_header(),update=True)
        hdulist[0].header.add_history('Added WCS header keys')
        hdulist.close()

    return w

def main():
    """ If run from termial """
    import sys
    import os
    
    if len(sys.argv)<3 :
        print('Usage : {0} X_Y_Ra_DecFile FitsImage'.format(sys.argv[0]))
        print('where,')
        print('X_Y_Ra_DecFile = Should contain 4 columns Xcoord Ycoord Ra(deg) Dec(deg) of atleast more than 4 stars')
        print('FitsImage = Is the fits image, whose header you want to update with WCS coordinate coefficents')
        exit(1)
    if not os.path.isfile(sys.argv[1]):
        print('Error: Text File {0} not found '.format(sys.argv[1]))
        exit(1)
    if not os.path.isfile(sys.argv[2]):
        print('Error: Fits image {0} not found '.format(sys.argv[2]))
        exit(1)
    #--End of Sanity check

    try :
        w=fitWCSobject(sys.argv[1],updatefits=sys.argv[2])
    except :
        raise
    else :
        print('Sucessfully updated WCS in {0}'.format(sys.argv[2]))


if __name__ == "__main__":                                                                                                 
    main()

        
        
