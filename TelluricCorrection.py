#!/usr/bin/env python
#This script will help in telluric correction of spectra with standard star spectra interactively.
import readline
import numpy as np
import numpy.ma
import scipy.interpolate as interp
import scipy.constants
import matplotlib.pyplot as plt

def BlackBodyFlux(Wavelengths,Temp):
    """ Returns an array of median scaled flux from balckbody correponding to the input list of wavelengths, and given Temperature 
    Input Wavelengths : a list or  numpy array of wavelengths in Angstrom
          Temp        : Temperature of black body in Kelvin  """
    h=scipy.constants.h
    k=scipy.constants.k
    c=scipy.constants.c
    WLfactor=10**-10   #The factor to convert input wavelength in Angstrom to SI unit meter
    Wavelengths=np.array(Wavelengths)*WLfactor
    numerator=2*h*c**2
    expfactor=h*c/(k*Temp*Wavelengths)
    denom=np.power(Wavelengths,5)*(np.exp(expfactor)-1)
    BBflux=numerator/denom
    return BBflux/np.median(BBflux)   #Return the median scaled BB spectrum

def AskAndMaskSpectra(spec):
    """ Displays the spectra in index coordinates and masks user requested points """
    if not np.ma.isMaskedArray(spec): spec=np.ma.array(spec)

    while True:
        plt.close()
        plt.plot(spec[:,1])    
        plt.plot(spec[:,1],'.')
        plt.grid(True)
        plt.show(block=False)
        print('-'*30)
        print("Following example shows various input option available to you.")
        print(" 215 218    #This will remove points from 215 to 217. (including 215 and 217)")
        print(" c 215 218  #This will unmask, if already masked the points from 215 to 217")
        print(" q          #This will finish the masking and exit the masking subroutine")
        choice=raw_input("Enter your option : ").strip(' ')
        print('-'*30)
        try:
            if choice[0] == 'q': #Exiting
                break
            elif choice[0] == 'c': #clear mask
                start=int(choice.split()[1])
                end=int(choice.split()[2])
                spec.mask[start:end,:]=np.ma.nomask
            elif len(choice.split()) == 2:  #To mask
                start=int(choice.split()[0])
                end=int(choice.split()[1])
                spec[start:end,:]=np.ma.masked
            else :
                print("Error: Unknown choice, please retype correctly.")
        except (IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")
    #Now return the masked array.
    return spec
    
def AskAndCleanStellarLines(spec):
    """ Interactively helps user to fit a convolved function of standard star spectral line and remove it """
    print("Not yet implemented")
    return spec

def TelluricCorrect(Sci,Std):
    """ Interactively asks user to align, scale, and remove spectral lines to correct telluric lines and fringes.
    Input : Sci is the 2d numpy array with 1st column wavelength and 2nd column flux of Science target
    Input : Std is the 2d numpy array with 1st column wavelength and 2nd column flux of Standard star """

    ScaleSci=1.0/np.median(Sci[:,1])  #Default scaling by median
    ScaleStd=1.0/np.median(Std[:,1])
    Sci=np.ma.array(Sci)  #Converting to masked array for future
    Std=np.ma.array(Std)
    StdWLshift=0  # Wavelength shift in Standard star spectra
    while True:
        #First we generate a spline interpolation of standard star spectra
        Std_tck=interp.splrep(np.ma.compressed(Std[:,0])+StdWLshift,np.ma.compressed(Std[:,1])*ScaleStd,s=0)
        #Now we generate the binned values of std spectra for science star.
        Std4Sci=interp.splev(Sci[:,0],Std_tck,der=0)
        
        plt.close() #Close any previously open plot window
        plt.figure(1)
        plt.subplot(211) #Top plot in the two layer subplot window
        plt.plot(Sci[:,0],Sci[:,1]*ScaleSci,label="Science")
        plt.plot(Sci[:,0],Std4Sci,label="Telluric Std")
        plt.legend()
        plt.grid(True)
        # Now generate the divided spectra
        CorrectedSci=(Sci[:,1]*ScaleSci)/Std4Sci
        plt.subplot(212) #And plot in the second subplot
        plt.plot(Sci[:,0],CorrectedSci,label=r'$\frac{Science}{Telluric}$')
        plt.legend()
        plt.grid(True)
        #Display
        plt.show(block=False)
        
        print('-'*40)
        print("Following examples shows the options to chose :")
        print(" sh 23     #To shift the Telluric spectra by +23 Angstrom")
        print(" sc 1.1    #To scale the Telluric spectra by factor of 1.1")
        print(" m         #To start a window on masking regions in spectra of Telluric Std spectra")
        print(" f         #To start a window on fitting and removing stellar lines in Telluric Std spectra")
        print(" q         #To exit, returning the latest telluric spectra and corrected science target")
        choice=raw_input("Enter your choice : ").strip(' ')
        print('-'*40)
        try:
            if choice[0] == 'q': #Exiting
                break
            elif choice[0] == 'f': #Stellar line fitting..
                Std=AskAndCleanStellarLines(Std)
            elif choice[0] == 'm': #Masking of unwanted region in std spectra
                Std=AskAndMaskSpectra(Std)
            elif choice[0:2] == 'sc': #Scaling of Telluric spectra
                print("Previous scaling factor was %f"%(ScaleStd))
                if float(choice.split()[1]) != 0 : ScaleStd*=float(choice.split()[1])
                else: print("You cannot scale by 0 factor")
                print("New scaling factor is %f"%(ScaleStd))
            elif choice[0:2] == 'sh': #Shifting of Telluric spectra
                print("Previous spectral shift was %f Angstrom"%(StdWLshift))
                StdWLshift+=float(choice.split()[1])
                print("New spectral shift is %f Angstrom"%(ScaleWLshift))
            else :
                print("Error: Unknown choice, please retype correctly.")
        except (IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")
    #Return the cleaned spectra
    return CorrectedSci,Sci[:,0],Std
            
            

def AskAndLoadAsciiData(toask="Enter Filename : ",Skipascii=69):
    """ Asks user ascii filename and return the table loaded as numpy array. 
    toask is string which contains the prompt to ask. 
    Skipasciii is an integer which tells the number of initial lines to skip while reading the data from file """

    while True:  #Keep asking till a file is properly loaded
        try :
            filename=raw_input(toask).strip(' ')
            Data=np.loadtxt(filename,skiprows=Skipascii)
            break
        except IOError :
            print("Error: Cannot find the file %s, Please give the correct file name."%(filename))
        except ValueError as e:
            print("Non float values in file? ")
            print(e)
            print("Present number of lines to skip in ascii file is %d "%Skipascii) 
            Skipascii=int(raw_input("Enter number of initial lines to skip in file:").strip(' '))
    return Data,filename

def is_number(s):   # A funtion to check wheter string s is a number or not.
    try:
        float(s)
        return True
    except ValueError:
        return False


def main():
    """ Asks user for the standard star spectra and science target spectra files and start telluric correction """
    StdStarRaw,stdfname=AskAndLoadAsciiData(toask="Enter the name of standard star spectrum ascii file : ",Skipascii=69)
    ScienceStarRaw,scifname=AskAndLoadAsciiData(toask="Enter the name of science star spectrum ascii file : ",Skipascii=69)
    CorrSci,CorrSciWL,Telluric=TelluricCorrect(ScienceStarRaw,StdStarRaw)
    BBtemp='SkipCC'
    BBtemp=raw_input("Enter Temperature of Std star (Kelvin) to correct continuum (to skip press enter): ")
    if is_number(BBtemp): #User entered the temperature of Std star
        T=float(BBtemp)
        CorrSci=CorrSci*BlackBodyFlux(CorrSciWL,T)
    else:
        print("Not doing continuum correction using BB curve.")
    plt.close()
    plt.plot(CorrSciWL/10000.0,CorrSci,drawstyle='steps',color='black')
    plt.xlabel(r'Wavelength ($\mu m$)')
    plt.ylabel('Normalised counts')
    plt.show()
    #Save the generated spectras values as numpy arrays
    scifname=scifname.split('.')[0]
    stdfname=stdfname.split('.')[0]
    np.save(scifname+'_Y_',np.array(np.ma.compressed(CorrSci)))  #Adding extra np.array to fix the bug in numpy till a fix 
    np.save(scifname+'_X_',np.array(np.ma.compressed(CorrSciWL)))
    np.save(stdfname+'_XY_',np.array(np.ma.compress_rows(Telluric)))

if __name__ == "__main__":
    main()
