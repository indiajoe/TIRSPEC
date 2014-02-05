#!/usr/bin/env python
# The faster version of the Generate4mRamp script uses the module numexpr
# So you need to install this module before using this script.
#-----------------------------
# This code is to generate a clean slope image from up-the-ramp readout images.
# Feel free to modify and share. Code is under GNU GPL v3 or higher versions.
# Author: J. P. Ninan (indiajoe@gmail.com)
#----------------------------------------------------------

from __future__ import division
#import matplotlib.pyplot as plt  #Debugggg
import numpy as np
import numpy.ma
import pyfits
import glob
import gc
import numexpr as ne

# By default it will use minnimum of 8 or number of cores you have on your machine for evaluations
# Or you can change to any number of cores by uncommenting the line below
# ne.set_num_threads(16)  

RampNDRsuffix='-debug*.fits'   # If changed, change in the function LoadDataCube() (line #38) also where we extract integer.
FullTile=(0,1024,0,1024)
LinearThresh=10500    #Threshold at which ramp looses linearity
#SaturThresh=12500     #Threshold at which signal is almost hard saturated.

def filelist(inputstring):
    """ Returns the filenames as a python list of strings.
       If the input string is prefixed with '@' it will load the files names from inside that text file name. 
       Otherwise the input string is expected to be comma separated input of filenames"""
    if inputstring[0] == "@" : 
        try :
            infile=open(inputstring[1:],'r')
            output=[fls.rstrip() for fls in infile]
            infile.close()
        except IOError :
            print("Error: Cannot open the file "+inputstring[1:]) 
            return []
    else : output=inputstring.split(',')
    return output

def CheckNDRexist(imgname,RampNDRsuffix,NoofNDRs=None):
    """ Checks the Non distructive Readouts (NDR) exists in the current directory for given image.
        Only if they exist we should proceed with Slope calculation from them """
    NDRfiles=glob.glob(imgname[:-5]+RampNDRsuffix)
    if NoofNDRs==None : NoofNDRs=pyfits.convenience.getval(imgname,'NDRS')
    if len(NDRfiles) == NoofNDRs :
        return True
    else:
        return False

def LoadDataCube(imgname,RampNDRsuffix,NoofNDRs=None,tile=FullTile):
    """ Loads the up-the-ramp readout of the input image tile to a single ndarray cube 
        Cube's coordinates are (time,X,Y)
        NoofNDRs is the number of NDRs to load into cube (type =int). By default everything is loaded."""
    NDRfilesT=glob.glob(imgname[:-5]+RampNDRsuffix)
    NDRfiles=sorted(NDRfilesT, key=lambda k: int(k[len(imgname)-5+7:-5]))
    if NoofNDRs==None : NoofNDRs=len(NDRfiles)

    datacube=np.empty((NoofNDRs,tile[1]-tile[0],tile[3]-tile[2]),dtype=np.float32)
    for i,img in enumerate((pyfits.getdata(NDR).astype(np.float32)[tile[0]:tile[1],tile[2]:tile[3]] for NDR in NDRfiles[:NoofNDRs])) : datacube[i,:,:]=img   # This iterator method is ~12 times faster than the equivalent commented line below !!
#    datacube=np.array([pyfits.getdata(NDR).astype(np.float32)[tile[0]:tile[1],tile[2]:tile[3]] for NDR in NDRfiles[:NoofNDRs]])

    return datacube

def AverageDataCube(imglist,RampNDRsuffix,NoofNDRs=None,tile=FullTile):
    """ Loads the data cubes of the files in the input imglist (a python list of string file names) 
    and return the ndarray cube with averaged data 
    NoofNDRs is the number of NDRs to load into cube (type =int). By default everything is loaded."""
    imgcube=LoadDataCube(imglist[0],RampNDRsuffix,NoofNDRs=NoofNDRs,tile=tile)
    if len(imglist) > 1 :  #Then Average the cubes..
        for img in imglist[1:]:  imgcube+=LoadDataCube(img,RampNDRsuffix,NoofNDRs=NoofNDRs,tile=tile)
        imgcube=imgcube/len(imglist)
    return imgcube

def CRhitslocation(imgcube,thresh=10):
    """ Returns the array of positions at which Cosmic Ray hits occurred.
    Input the  image data cube, and optional threshold (thresh= ) to detect the CR hit convolution image.
    The output of the locations of CR hit are in the tuple of 3 arrays format (Time,X,Y).
    """
    term1=imgcube[:-3,:,:]
    term2=imgcube[1:-2,:,:]
    term3=imgcube[2:-1,:,:]
    term4=imgcube[3:,:,:]
    convimg=ne.evaluate("term1 - 3*term2 + 3*term3 - term4") #Convolution of image with [1,-3,3,-1] along time axis
    RN=4.3   #ReadNoise
    # if np.ma.isMaskedArray(imgcube) :
    #     convimg=imgcube.data[:-3,:,:] -3*imgcube.data[1:-2,:,:] +3*imgcube.data[2:-1,:,:] -imgcube.data[3:,:,:]
    # else:
    #     convimg=imgcube[:-3,:,:] -3*imgcube[1:-2,:,:] +3*imgcube[2:-1,:,:] -imgcube[3:,:,:] #Convolution of image with [1,-3,3,-1] along time axis
#    stdconv=np.std(convimg,axis=0)
    stdconv=np.median(convimg,axis=0)  #Nameing the variable now itself to conserve memory
    np.median(np.abs(convimg-stdconv),axis=0,out=stdconv)  # MAD for robustnus
#    stdconv=np.ma.median(np.abs(convimg-np.ma.median(convimg,axis=0)),axis=0)  # MAD for robustnus
    thresh*=1.4826
    #We should remove all the points where number of data points were less than 5
    stdconv[np.where(np.ma.count(imgcube,axis=0) < 5)]=99999
    MinNoise=3*np.sqrt(RN**2 +(3*RN)**2 +(3*RN)**2 + RN**2) # 3 Sigma from ReadNoise
    stdconv[np.where(stdconv < MinNoise/thresh)]= MinNoise/thresh   #This is 3 sigma limit of RN
    #Find the edges of ramp jumps.
    if np.ma.isMaskedArray(imgcube) :
        (T,X,Y)=np.ma.where(np.ma.array(convimg,mask=np.ma.getmaskarray(imgcube[3:,:,:])) > thresh*stdconv) 
    else: (T,X,Y)=np.ma.where(convimg > thresh*stdconv)
    T=T+2  #Converting to the original time coordinate of imgcube by correcting for the shifts. This T is the first pixel after CR hit
    return (T,X,Y)

def ReplaceCRhits(imagelist,NoofNDRs=None,tile=FullTile):
    """ Outputs an image cube with of all pixels hit by Cosmic Rays in first image replaced with those from the next clean image"""
    if type(imagelist)==str : imglist=filelist(imagelist)   #List of images
    else : imglist=imagelist
    imgcube=LoadDataCube(imglist[0],RampNDRsuffix,NoofNDRs=NoofNDRs,tile=tile)
    T,X,Y=CRhitslocation(imgcube,thresh=10) #Locations of CR hit
    Tlength=imgcube.shape[0]
    XYlist=np.unique(zip(X,Y))
    for img in imglist[1:]:  
        nextcube=LoadDataCube(img,RampNDRsuffix,NoofNDRs=NoofNDRs,tile=tile)
        foundk=[]
        for k in range(len(XYlist)) :
            i,j=XYlist[k]
            t,junkx,junky=CRhitslocation(nextcube[:,i,j].reshape((Tlength,1,1)),thresh=10) #Locations of CR hit
            if len(t)==0 :  # no CR hits in this image at this i,j pixel
                imgcube[:,i,j]=nextcube[:,i,j]
                foundk.append(k)
        #Removing all the pixels we found clean readouts from the XYlist to search in next image
        XYlist=np.delete(XYlist,foundk,axis=0) 
        if len(XYlist)==0 : break  # Break if all CR hits are replaced

    for i,j in XYlist : print("CR hit couldn't be replaced for the pixel : (%d,%d)" %(i,j))   #For debugging 
    
    return imgcube
        
def AverageWithCRrej(imagelist,NoofNDRs=None,tile=FullTile):
    """ Outputs the average of all the data cubes after rejecting the pixels hit by Cosmic Rays"""
    if type(imagelist)==str : imglist=filelist(imagelist)   #List of images
    else : imglist=imagelist
    sumcube=0
    N=np.zeros((tile[1],tile[3]),dtype=np.int8)
    for img in imglist:
        imgcube=LoadDataCube(img,RampNDRsuffix,NoofNDRs=NoofNDRs,tile=tile)
        imgcube=np.ma.array(imgcube,fill_value=0) #Keeping fill value to zero for removing the CR hit pixels
        T,X,Y=CRhitslocation(imgcube,thresh=10) #Locations of CR hit
        for i,j in np.unique(zip(X,Y)): imgcube[:,i,j]=np.ma.masked
        sumcube+=imgcube.filled()
        N[~np.ma.getmaskarray(imgcube[0,:,:])]+=1
    #If there are still some pixels which were hit by CR in every single frame we shall replace it with the last frame's values
    X,Y=np.where(N==0)
    for i,j in np.unique(zip(X,Y)):
        print("CR hit couldn't be removed for the pixel : (%d,%d)" %(i,j))   #For debugging 
        sumcube[:,i,j]=imgcube.data[:,i,j]
        N[i,j]=1
    
    return sumcube/N


def FitSlope(imgcube,time, CRcorr=False):
    """ Fits a slope by linear regression along the axis=0, with respect to time and returns the slope and constant for 
    each pixel as a 2d Matrix
    The linear fitting function is  y = alpha + beta *x
    Equations based on https://en.wikipedia.org/wiki/Simple_linear_regression#Numerical_example


    Parameters:
    -----------
    imgcube   : Masked numpy 3d array.
               Time axis should be axis=0
               The points which shouldn't be used for straight line fitting should be masked out by numpy's ma module.

    time      : 1d numpy array
               The time corresponding to each slice along the time axis of the data. The slope calculated will be in 
               units of this time.

    CRcorr    : Boolean True / False  (optional)
               Default is False. If CRcorr is True, then Cosmic ray hits are detected and removed while calculating the slope.
               Weighted average of the slopes of the regions unaffected by Cosmic ray hits are used to fill beta
               The constant alpha is taken to be that of the first fitted good region's.

    Returns
    -----------
    (beta,alpha) : masked (2d numpy array, 2d numpy array)
                  The straight line fit is of the equation  y = alpha + beta *x
                  The slope beta and constant alpha is for each pixel is returned as two 2d numpy arrays.

    """
    #The linear regression fit.  y = alpha + beta *x
    #Equation notations based on https://en.wikipedia.org/wiki/Simple_linear_regression#Numerical_example
    #Variables are 2d matrices corresponding to 2d array of pixels.
    tshape=tuple(np.roll(imgcube.shape,-1))  # Creating the tuple to resize the time array to 3d cube. We will have to do a (2,0,1) permutation later. 

    # Sx=np.ma.array(np.ones(imgcube.shape,dtype=np.int)*time[:,np.newaxis,np.newaxis],mask=np.ma.getmaskarray(imgcube)).sum(axis=0)
    # Sxx=np.ma.array(np.ones(imgcube.shape,dtype=np.int)*(time**2)[:,np.newaxis,np.newaxis],mask=np.ma.getmaskarray(imgcube)).sum(axis=0)
    Sx=np.ma.array(np.transpose(np.resize(time,tshape),(2,0,1)),mask=np.ma.getmaskarray(imgcube)).sum(axis=0,dtype=np.float64)
    Sxx=np.ma.array(np.transpose(np.resize(np.square(time),tshape),(2,0,1)),mask=np.ma.getmaskarray(imgcube)).sum(axis=0,dtype=np.float64)

    Sy=imgcube.sum(axis=0,dtype=np.float64)
    Syy=(np.square(imgcube)).sum(axis=0,dtype=np.float64)
    Sxy=(imgcube*time[:,np.newaxis,np.newaxis]).sum(axis=0,dtype=np.float64)
    n=np.ma.count(imgcube,axis=0)   #number of points used in fitting slope of a pixel
    
    beta= ne.evaluate("(n*Sxy - Sx*Sy)/ (n*Sxx - Sx*Sx)")
    alpha= ne.evaluate("Sy/n - beta*Sx/n ")

    #mask beta and alpha where n < 2
    beta=np.ma.masked_where(n<2,beta)
    alpha=np.ma.array(alpha,mask=np.ma.getmaskarray(beta))

    #If cosmic ray hit correction is asked to do and total exposure of image was more than 4 sec
    if CRcorr == True and imgcube.shape[0] > 5 :
        T,X,Y=CRhitslocation(imgcube,thresh=10)  #Locations of CR hit
        #Masking the first pixel after a CR hit.
        imgcube[(T,X,Y)]=np.ma.masked
        #Loop over each pixel (X,Y) with cosmic ray hits
        print("Number of Cosmic Ray hits : %d " %(len(T)))
        if len(T) > 1000 :  #Sanity check
            print("More than 1000 CR hits are not possible and worth removing. \n So better discard this image. Skipping CR hit removal.")
            return (beta,alpha)   #Exiting this function skipping CR hit removal...
        for i,j in np.unique(zip(X,Y)):
            localbeta=[]
            localalpha=[]
            localn=[]
            localsigma=[]
            #loop over each sections of the ramp.
            slices=np.ma.notmasked_contiguous(imgcube[:,i,j])
            if slices is None :  #When no unmasked pixels exist
                print('Couldnot remove CR (insufficent data) at %d %d '%(i,j))
                continue 
            for k in range(len(slices)) :
                n=len(imgcube[:,i,j][slices[k]])
                if  n > 3 : #At least 4 points are there to calculate slope
                    time=np.arange(slices[k].start,slices[k].stop)
                    Sx=time.sum(dtype=np.float64)
                    Sxx=(np.square(time)).sum(dtype=np.float64)
                    Sy=imgcube[:,i,j][slices[k]].sum(dtype=np.float64)
                    Syy=(np.square(imgcube[:,i,j][slices[k]])).sum(dtype=np.float64)
                    Sxy=(imgcube[:,i,j][slices[k]]*time).sum(dtype=np.float64)
                    #append localbeta, localalpha, localn and localsigma
                    localbeta.append((n*Sxy - Sx*Sy)/ (n*Sxx - Sx**2))
                    localalpha.append(Sy/n - localbeta[-1]*Sx/n)
                    localn.append(n)
                    Se2=(n*Syy -Sy**2- localbeta[-1]**2 *(n*Sxx- Sx**2))/(n*(n-2))
                    localsigma.append(np.sqrt(n*Se2/(n*Sxx-Sx**2)))   #Std dev Error of slope beta
            #calculate the average beta with weights 1/localsigma 
            if len(localsigma) > 0 : 
                beta[i,j]=np.average(localbeta,weights=1.0/np.asarray(localsigma))
                alpha[i,j]=localalpha[0]
#            except ZeroDivisionError :
            else :
                print('Couldnot remove CR hit (insufficent data) at %d %d '%(i,j))
#                print("LocalSigma:",localsigma)
#                print("Localbeta:",localbeta)

    return (beta,alpha)


def Generate_DarkTemplate(longexp,shortexp,extent=35,LThresh=LinearThresh,outputfile=None,tile=FullTile) :
    """ Generates a numpy array cube containing templates of the initial nonlinear fall 
    part of the dark current during up-the-ramp readout.

    Parameters:
    -----------
    longexp   : String 
               It is a  list of long exposure dark frames in a comma separated format.
               Or it could be a filename prefixed with '@' symbol containing the filesnames in each line of a textfile.
               Example: "LongDark1.fits,LongDar2.fits"      OR   "@Darklist.txt"

    shortexp  : String
               It is a  list of short (but longer than 'extent') exposure dark frames in a comma separated format.
               Or it could be a filename prefixed with '@' symbol containing the filesnames in each line of a textfile.
               Example: "ShortDark1.fits,ShortDar2.fits"      OR   "@Darklist.txt"
    
    extent    : int
               Integer giving the time in units of Non-Destructive readouts, up to which the non-linear part of dark 
               current template has to be made.
               Example: 35

    LThresh   : 2d numpy array    (optional)
               Lower Threshold value of the pixel count for each pixel above which we have to fit the correction for non-linearity
    
    outputfile: String   (optional)
               It is the name of the output .npy file which should be made to save the numpy ndarray cube of templates.
               If no input is given, then the function will return the ndarray template.
    
    tile      : Tuple   (optional)
               The tile of the image we should generate the dark current template.
               To be used for large number of NDRs and if you don't have enough RAM memory to load everything together.
               Format of tuple is (Xbegin,Xend,Ybegin,Yend)  Eg: (0,512,0,512) for a quadrant. 

    Returns
    -----------
    If outputfile input is provided then the template as an ndarray cube will be  written in current directory as a python pickle file 
    with that name.
    Else the ndarray cube will be returned as output.
           
               """
    
    longfiles=filelist(longexp)   #List of long exposure dark files


    imgcube=AverageDataCube(longfiles,tile=FullTile) #Average the data cubes

    #Now we have to fit a straight line to region beyond input time "extent" up-the-ramp to get the linear slope.
    imgcube=np.ma.masked_greater(imgcube[extent:],LThresh)  #Trimming the initial part and masking away the saturation
    time=np.arange(extent,imgcube.shape[0]+extent,dtype=np.int)  #1-D time axis in units of NDR number.
    
    #Now Calculating the linear regression fit.  y = alpha + beta *x
    beta,alpha = FitSlope(imgcube,time)
   
    #Now proceed to load the short exposure dark frames and subtract the slope.
    shortfiles=filelist(shortexp)  #List of short exposure dark files
    imgcube=AverageDataCube(shortfiles,tile=FullTile)

    time=np.arange(imgcube.shape[0],dtype=np.int)  #1-D time axis in units of NDR number.
    #Subtracting the y = alpha + beta *x component from the dark cube.
    #FUTURE WORK: To be split into -= and += lines if memory shortage is noticed.
    imgcube=(( imgcube-alpha[np.newaxis,:] ) - beta[np.newaxis,:]*time[:,np.newaxis,np.newaxis]).astype(np.float32)
    if outputfile == None: return imgcube
    else: 
        np.save(outputfile,imgcube.data)
        return
        
    
def ApplyNLdarkcorr(image,darkcube,tile=FullTile):
    """ Subtracts only the Non-Linear part of the initial few seconds of Dark current given in darkcube from the image's data cube

    Parameters:
    -----------
    image     : numpy 3d array cube / String 
               It is the file name of the fits image which has to be dark subtracted.
               Example: Objectframe.fits

    darkcube  : numpy ndarray
               It is the template of dark current cube which has to be subtracted from the image cube. (output of Generate_DarkTemplate)
               It is important that dark current cube was generated from dark frames taken in same condition as Object frame.  
               It is expected from user to make sure that the origin of dark template in X,Y plane is the same as the 
               origin of the tile requested in the tile input.


    tile      : Tuple    (optional)
               The tile of the image we should generate the image cube and subtract dark current template.
               To be used for large number of NDRs and if you don't have enough RAM memory to load everything together.
               Format of tuple is (Xbegin,Xend,Ybegin,Yend)  Eg: (0,512,0,512) for a quadrant. 
               It is expected from user to make sure that the origin of dark template in X,Y plane is the same as the 
               origin of the tile requested in the tile input.
             
    Returns
    -----------
    imgcube: numpy ndarray
            It returns the data cube for the image after subtracting the dark current Non-linear template.
            """
    
    if type(image)==str : imgcube=LoadDataCube(image,RampNDRsuffix,tile=tile)
    else : imgcube=image
    #We want to apply the NL dark correction only where there is an overlap of dark current template over the image.
    # It is expected from user side to make sure that the origin of dark template in X,Y plane is the same as the origin of the tile
    # requested in the tile input.
    # Following X,Y variables are in the big full image pixel coordinates. (Not useful, but may be easier to generalise later)
    Xbegin=tile[0]  
    Xend=min(tile[1],Xbegin+darkcube.shape[1])
    Ybegin=tile[2]
    Yend=min(tile[3],Ybegin+darkcube.shape[2])
    
    timebegin=0
    timeend=min(imgcube.shape[0],darkcube.shape[0])
    print('Applying NL-dark correction in region : ('+str(Xbegin)+':'+str(Xend)+','+str(Ybegin)+':'+str(Yend)+')')
    #Now converting back in to the extracted cube coordinates Xbegin-tile[0]=0, Ybegin-tile[2]=0
    tb=timebegin
    te=timeend
    xb=Xbegin-tile[0]
    xe=Xend-tile[0]
    yb=Ybegin-tile[2]
    ye=Yend-tile[2]

    imgcube[tb:te,xb:xe,yb:ye]-=darkcube[tb:te,xb:xe,yb:ye]
    
    return imgcube

def NonlinearityCorrCoeff(image,dark,order=3,LThreshFactor=0.7,UThresh=None,LThresh=None,outputfile=None,tile=FullTile):
    """ Outputs the coefficients of polynomial as 2d array for pixels. These polynomials map from observed count to non-linearity corrected actual count 

    Parameters:
    -----------
    image     : numpy 3d array cube / String 
               If it is string then it should be file name of the saturated fits image which has all the pixels saturated, and has to 
               be used for estimating the nonlinearity correction polynomial coefficients
               Example: Saturatedframe.fits
               Or If it is numpy 3d array cube. It should be the loaded data cube to use. Ex: Cube after Cosmic Ray Removal.
               IMP: It is very important that there should not be any Cosmic Ray hits in this image.

    dark      : numpy 3d array cube / String 
               If it is string then it should be file name of the dark fits image which has same number of NDR frames as the object 
               image.  Example: Darkframe.fits
               Or If it is numpy 3d array cube. It should be the loaded data cube to use. 
               Ex: Average Cube after Cosmic Ray removal in Dark frames.
               IMP: It is very important that there should not be any Cosmic Ray hits in this dark.

    order     : Int   (optional)
               It is the order of polynomial to fit for non-linear counts to  actual count function. 
               Don't make it very large, < 4 is good enough. Default is 3.
        
    LThreshFactor : float   (optional)
               It is the fraction of Uthresh, up to which detector is expected to be linear. The actual slope after dark correction 
               is calculated in the region up to Uthresh*LthreshFactor. Default value of this factor is 0.7 . 
               i.e. 0.7 factor of the Up-the-Ramp up to UThresh is expected to be linear.

    UThresh   : 2d numpy array   (optional)
               Upper Threshold values of the pixel count for each pixel up to which we have to fit the correction for non-linearity.
               It is better to keep this slightly below hard saturation.
               If None is given, UThresh is estimated from the global minima in second derivative of the up the ramp.

    LThresh   : 2d numpy array   (optional)
               Lower Threshold values of the pixel count for each pixel above which we have to fit the correction for non-linearity.
               Keep this around 2nd NDR readout value to remove the highly NonLinear Dark current region out of polynomial's fit.
               If None is given, LThresh is taken to be 2nd NDR readout of the up the ramp.

    outputfile: String   (optional)
               It is the name of the output .npy file which should be made to save the numpy ndarray cube of coefficients.
               If no input is given, then the function will return the ndarray 3d cube of coefficients.

    tile      : Tuple   (optional)
               The tile of the image we should estimated the coefficients of each pixel.
               To be used for large number of NDRs and if you don't have enough RAM memory to load everything together.
               Format of tuple is (Xbegin,Xend,Ybegin,Yend)  Eg: (0,512,0,512) for a quadrant. 

    Returns
    -----------
    coeffs : 3d coefficient arrays, where the axis=0 contains coefficients corresponding to polynomial (length= order of polynomial + 1)
             The dimension of the output 3d cube will be (order+1,X,Y)
             If Along the axis=0 the coeffs are coeffs[0,X,Y],coeffs[1,X,Y],coeffs[2,X,Y].. 
             then the polynomial is coeffs[0,X,Y]*x**order +coeffs[1,X,Y]* x**(order-1) +...+coeffs[order+1,X,Y]

    LThresh: 2d numpy array
             Lower Threshold values of the pixel count for each pixel above which we have fitted the correction for non-linearity.
             IMP: Never apply the correction polynomial coefficients below this value.

    UThresh: 2d numpy array
             Upper Threshold values of the pixel count for each pixel up to which we have fitted the correction for non-linearity.
             IMP: Never apply the correction polynomial coefficients above this value.
                         
             """
    if type(image)==str : imgcube=LoadDataCube(image,RampNDRsuffix,tile=tile)
    else : imgcube=image
    if type(dark)==str : darkcube=LoadDataCube(dark,RampNDRsuffix,tile=tile)
    else : darkcube=dark
    if imgcube.shape != darkcube.shape :
        print("Image cube and Dark cube doesn't match in dimensions: Image is "+str(imgcube.shape)+" and Dark is "+str(darkcube.shape))
        return

    if UThresh == None : #No input Upper threshold given. So use the global minima of second derivative
        ne.evaluate("imgcube-darkcube",out=imgcube)  #Removing the dark current temporarily for now
        term1=imgcube[:-2,:,:]
        term2=imgcube[1:-1,:,:]
        term3=imgcube[2:,:,:]
        imgcube2ndD=ne.evaluate("term1 -2*term2 + term3")  #Convolving with [1,-2,1] for second derivative
        imgcube2ndD=imgcube2ndD[2:,:,:]  #removing the first two entries before finding minima
        Cuttoff=np.argmin(imgcube2ndD,axis=0)+2+2   #taking global minima
        del imgcube2ndD      #Delete to free memory
        ne.evaluate("imgcube+darkcube",out=imgcube)  #Putting back the dark current for finding upper threshold
        i,j=np.ogrid[0:imgcube.shape[1],0:imgcube.shape[2]]
        UThresh=imgcube[Cuttoff,i,j]       #Extracting the threshold values by fancy indexing
        
        
    #Now we have to fit a straight line to region below the Linear Threshold in up-the-ramp to get the linear slope.
    imgcube=np.ma.masked_greater(imgcube,UThresh*LThreshFactor)  # Masking away the Non-linear part near saturation

    timelength=imgcube.shape[0]
    time=np.arange(timelength,dtype=np.int)  #1-D time axis in units of NDR number.
    
    #Now Calculating the linear regression fit for image-dark data.  y = alpha + beta *x
    beta,alpha = FitSlope(imgcube-darkcube,time)
    
    #Removing the previous mask, and masking out linear and hard saturation region.
    imgcube.mask = np.ma.nomask

#    imgcube=np.ma.masked_greater(imgcube,UThresh)  #For No lower bound case

    if LThresh==None: #Lets define the lower bound of the region where polynomial fit is valid. 
        LThresh=imgcube[2,:,:]  # Change the index here to choose any slice below which to be masked away before polynomial fit.

    imgcube=np.ma.masked_less(np.ma.masked_greater(imgcube,UThresh),LThresh)
    
    #fitting of polynomial between original data and straight line.
    coeffs=np.zeros((order+1,imgcube.shape[1],imgcube.shape[2]))  #An array to store coeffs
    coeffs[-2,:,:]=1   #Setting by default all correction polynomial  to be simply y=x . i.e. no correction.
    #The following for loops are very slope. (To be replaced by Nadia's (STSCI) astropy pull request later)
    
    for i in range(imgcube.shape[1]) :
        for j in range(imgcube.shape[2]) :
            # Fit only if there are order+1 points to fit, and the pixel is super saturated at end.
            if np.ma.count(imgcube[:,i,j]) > order +1 and imgcube[:,i,j].mask[-1] == True  :   
                y=np.ma.array(time,mask=np.ma.getmaskarray(imgcube[:,i,j]))
                y=(alpha[i,j] + beta[i,j]*y )+darkcube[:,i,j]  #Fitted slope + dark current
                coeffoffit,residuals,rank,singular_values,rcond=np.ma.polyfit(imgcube[:,i,j],y,order,full=True)
                if rank != order+1 : print("Polyfit may be poorly conditioned at pixel : (%d,%d)" %(i,j))
                else : coeffs[:,i,j]=coeffoffit
                # plt.plot(y,imgcube[:,i,j])  # Debuggg
                # plt.show()                   #Debuggg 
            else: print("Non linear correction not applicable for pixel : (%d,%d)" %(i,j))   #For debugging 
    
    if outputfile == None: 
        return (coeffs,LThresh,UThresh)
    else: 
        np.save(outputfile,coeffs)
        return
    
    
def ApplyNonLinearityCorrection(image,NLcoeffs,LThresh,UThresh,tile=FullTile):
    """ Applies the non linearity correction above LThresh upto UThresh of a image (or datacube)
         and returns a numpy ma masked cube, with all points above UThresh masked

    Parameters:
    -----------
    image     : numpy 3d array / String 
                If it is a string then it should be the file name of the fits image to load and process.
                Ex: StarImage.fits
                Or it could be the 3d data cube of up-the-ramp readout. Time axis should be axis=0
               
    NLcoeffs  : String or numpy 3d array
                Either the file name of the npy file which has the save non linearity correction coefficents.
                Or the 3d cube of coefficients from the output of NonlinearityCorrCoeff()

    LThresh   : 2d numpy array  
               Lower Threshold value of the pixel count for each pixel above which we have to fit the correction for non-linearity

    UThresh   : 2d numpy array 
               Upper Threshold value of the pixel count for each pixel upto which we have to fit the correction for non-linearity.
               It is better to keep this slightly below hard saturation. And where the coefficients are valid.

    tile      : Tuple   (optional)
               The tile of the image we should estimated the coefficients of each pixel.
               To be used for large number of NDRs and if you don't have enough RAM memory to load everything together.
               Format of tuple is (Xbegin,Xend,Ybegin,Yend)  Eg: (0,512,0,512) for a quadrant. 

    Returns
    -----------
    imgcube   : numpy ma masked 3d cube array
              Outputs the non linearity corrected 3d cube image and numpy masked ma array 
              and with all values which was above UThresh masked.

    """
    if type(image)==str : imgcube=LoadDataCube(image,RampNDRsuffix,tile=tile)
    else : imgcube=image

    # If NLcoeffs is a string it must be the .npy file of the non linear polynomial correction coefficients.
    if type(NLcoeffs)==str : NLcoeffs=np.load(NLcoeffs)   #Load the coeffs 3d array.

    order=NLcoeffs.shape[0]-1    # Order of polynomial
    
    UpperMask=np.ma.getmaskarray(np.ma.masked_greater(imgcube,UThresh))
    LowerMask=np.ma.getmaskarray(np.ma.masked_less(imgcube,LThresh))
    # Masking both linear and hard saturation region to apply Non linear correction
    imgcube=np.ma.array(imgcube,mask=np.ma.mask_or(UpperMask, LowerMask))
    
    summ=imgcube.copy()
    summ[~np.ma.getmaskarray(summ)]=0   #Setting all unmasked elements to zero
    for i in range(order+1):
        summ+=NLcoeffs[i][np.newaxis,:,:]*(np.power(imgcube,(order-i)))

    imgcube=summ
    #Removing all the maskes..
    imgcube.mask = np.ma.nomask
    #Putting back the upper hard Saturation region's Mask
    imgcube=np.ma.array(imgcube,mask=UpperMask)
    return imgcube

    
def main():
    import os
    import os.path
    import sys
    import re
    import pyfits.convenience

    DarkNDRS=112    # Specify the number of NDRS to be loaded for average Dark cube HERE....
    FullframeMax=46   # Maximum NDRS we can load as a single frame within memory constrains
    UthreshFactor=0.9   #Fraction of Uthresh to be actually used as upper cuttoff
    Uthreshnpyfile='/media/PlanetX/TIRSPEC1strun/ReductionWorkarea/20130620/UpperThreshold.npy'

    if len(sys.argv) < 2 or not os.path.isdir(sys.argv[1]):
        print("-"*30)
        print("Usage  : "+sys.argv[0]+" Night_Directory")
        print("Eg: "+sys.argv[0]+" /data/20130621/")
        print("-"*30)
        exit(1)
    os.chdir(sys.argv[1])
    listOFimgsT= [f for f in os.listdir('.') if re.match(r'^(?!(.*debug.*|.*[Tt]est.*|.*[Ff]ocus.*)).*\.fits', f)]
    listOFimgs= [f for f in listOFimgsT if CheckNDRexist(f,RampNDRsuffix)]
    for lostimgs in list(set(listOFimgsT)-set(listOFimgs)) :
        print('Warning: Discarded %s because of missing NDR files'%lostimgs)

    Darklist=[f for f in listOFimgs if re.match(r'^[Dd]ark.*', f)]
    Darks=[]
    #For Specific NDR long darks
    print('Darks with %d NDRS used for averaging'%DarkNDRS)
    for i in Darklist :
        if int(pyfits.convenience.getval(i,'NDRS'))>=DarkNDRS : 
            Darks.append(i)
            print(i,pyfits.convenience.getval(i,'NDRS'),pyfits.convenience.getval(i,'ITIME'))
    DarkAvgcube=AverageWithCRrej(Darks,NoofNDRs=DarkNDRS) #Average Dark with CR hits rejected
#    np.save('AVGdarkforthenight',DarkAvgcube)    #Saving to save time later
#    DarkAvgcube=np.load('AVGdarkforthenight.npy')
    #Load the latest saturation levels. (Upper threshold file)
    Uthresh=np.load(Uthreshnpyfile)  
    listwithoutdark=[f for f in listOFimgs if not re.match(r'^[Dd]ark.*', f) ]
    sortedlist=sorted(listwithoutdark, key=lambda k: int(k[:-5].split('-')[-1]))
    logfile=open('SlopeimagesLog.txt','w')
    for img in sortedlist:
        print(img)
        hdulist=pyfits.open(img)
        prihdr= hdulist[0].header
        NDRS2read=min(prihdr['NDRS'],DarkNDRS)
        for hkeys in  ['TARGET','TCOMMENT']:   #To capture bug of these keywords missing in header
            if hkeys not in prihdr : prihdr.update(hkeys,'-NA-')

        logfile.write('%s %s "%s" %d %.2f %s %s %s %s %s %s \n'%('Slope-'+img,prihdr['TIME'],prihdr['TARGET'],prihdr['NDRS'],prihdr['ITIME'],prihdr['UPPER'],prihdr['LOWER'],prihdr['SLIT'],prihdr['CALMIR'],prihdr['DATE'],prihdr['TCOMMENT']))
        if NDRS2read <= FullframeMax :
            tilelist=[(0,1024,0,1024)]  #Loading out as a single frame ; else we load each quadrent seperately
        else: tilelist=[(0,512,0,512),(0,512,512,1024),(512,1024,0,512),(512,1024,512,1024)]
        for sec in tilelist:
            gc.collect()
            print("Working on tile :"+str(sec))
            imgcube=LoadDataCube(img,RampNDRsuffix,NoofNDRs=NDRS2read,tile=sec)
            imgcube=np.ma.masked_greater(imgcube,Uthresh[sec[0]:sec[1],sec[2]:sec[3]]*UthreshFactor)
            imgcube-=DarkAvgcube[:NDRS2read,sec[0]:sec[1],sec[2]:sec[3]]  #Dark subtraction
            if imgcube.shape[0] >= 60 :  #Throw away first 5 pixels before fitting dlope
                imgcube=imgcube[5:,:,:]
                prihdr.add_history('Discarded first 5 NDRS')
            elif imgcube.shape[0] >= 20 and prihdr['LOWER']=='G' : #Spectras more than 20 NDRS
                imgcube=imgcube[5:,:,:]
                prihdr.add_history('Discarded first 5 NDRS') 
            elif imgcube.shape[0] >= 15 : #Atleast 15 NDRS are there...
                imgcube[:5,np.ma.count(imgcube,axis=0) >= 15]=np.ma.masked
                prihdr.add_history('Discarded first 5 NDRS of non saturating pixels')
            time=np.arange(imgcube.shape[0],dtype=np.int)
            beta,alpha=FitSlope(imgcube,time)#, CRcorr=True)
            hdulist[0].data[sec[0]:sec[1],sec[2]:sec[3]]=beta.data
            del imgcube
        prihdr.add_history('Dark subtracted using:DarkAvgcube')
        prihdr.add_history('Calculated slope removing saturation')# and CR hits')
        hdulist.writeto('Slope-'+img)
    logfile.close()
    print('Succesfully created slope images for this night \n Copy out all\033[91m Slope*\033[0m files')


if __name__ == "__main__":
    main()
