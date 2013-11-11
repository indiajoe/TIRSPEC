#!/usr/bin/env python
#This script will help in telluric correction of spectra with standard star spectra interactively.
import readline
import numpy as np
import numpy.ma
import scipy.interpolate as interp
import scipy.constants
import matplotlib.pyplot as plt

def RichardsonLucy(Signal,PSF,IniKernel,niter):
    """ By Richardson-Lucy Algorithm, retrun the deconvolution of Signal using input PSF.
    Input: Signal : 1d numpy array, which is to be deconvolved
           PSF    : 1d numpy array, which is the PSF to be used for deconvolution
           IniKernel: 1d numpy array, which is the initial guess of original signal
           niter    : integer, Number of iterations of RL deconvolution steps to be run """
    signal=Signal.copy()   #Copying to protect input
    psf=PSF.copy()
    inikernel=IniKernel.copy()
    # Algorithm and syntax based on http://en.wikipedia.org/wiki/Richardson-Lucy_deconvolution
    psf_hat=psf[::-1] #psf hat, following wikipedia's syntax
    for i in range(niter):
        est_conv = np.convolve(inikernel,psf,mode='same')
        relative_blur = signal/est_conv
        error_est = np.convolve(relative_blur,psf_hat,mode='same') 
        inikernel = inikernel* error_est
        
    return inikernel
    

def BlackBodyFlux(InWavelengths,Temp):
    """ Returns an array of median scaled flux from balckbody correponding to the input list of wavelengths, and given Temperature 
    Input InWavelengths : a list or  numpy array of wavelengths in Angstrom
          Temp        : Temperature of black body in Kelvin  """
    Wavelengths=InWavelengths.copy() #Copying to protect input
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

def AskAndMaskSpectra(Inspec):
    """ Displays the spectra in index coordinates and masks user requested points """
    spec=Inspec.copy() #Copying to protect input
    if not np.ma.isMaskedArray(spec): spec=np.ma.array(spec)

    while True:
        plt.close()
        plt.plot(spec[:,1])    
        plt.plot(spec[:,1],'.')
        plt.grid(True)
        plt.show(block=False)
        print('M'+'-'*30+'M')
        print("Following example shows various input option available to you.")
        print(" 215 218    #This will remove points from 215 to 217. (including 215 and 217)")
        print(" c 215 218  #This will unmask, if already masked the points from 215 to 217")
        print(" q          #This will finish the masking and exit the masking subroutine")
        choice=raw_input("Enter your option : ").strip(' ')
        print('M'+'-'*30+'M')
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


def ContinuumNormSpec(Inspec):
    """ Returns the continuum normalised spectra by interactively masking lines
    Input: spec : a 2d numpy array, with first column wavelength and second column flux 
    Output: Continuum normalise spec, flux of line in arbitrary unit and Fitted continuum """
    spec=Inspec.copy() #Copying to protect input
    continuum=np.ma.array(spec)
    useradd=0
    Domask=True
    while True:
        if Domask : 
            print("Mask away the regions which are not part of continuum")
            continuum=AskAndMaskSpectra(continuum)
            Domask=False
            plt.close()
        l=np.ma.count(continuum,axis=0)[1]  #Number of points to fit continuum
        sm=(l+np.sqrt(2*l))*6+useradd   #Smoothing factor for fitting.
        cont_tck=interp.splrep(np.ma.compressed(continuum[:,0]),np.ma.compressed(continuum[:,1]),s=sm)
        CleanCont=interp.splev(spec[:,0],cont_tck,der=0)
        scale=np.median(spec[:,1])
        plt.plot(spec[:,0],spec[:,1]/scale)
        plt.plot(spec[:,0],CleanCont/scale)
        plt.plot(spec[:,0],spec[:,1]/CleanCont)
        plt.grid(True)
        plt.show(block=False)
        print('-'*10)
        print("Present smoothing factor was %f "%(sm))
        print("Following example shows various input options available for you.")
        print(" 200   #Add 200 to default smoothing parameter")
        print(" m     #Mask more points for fitting continuum")
        print(" q     #Return this procedure with normalised spectra and eqw inside masked region")
        choice=raw_input("Enter your option : ").strip(' ')
        print('-'*10)
        try:
            if choice[0] == 'q': #Exiting 
                break
            elif choice[0] == 'm': #Continue Masking
                Domask=True
            elif is_number(choice) : # Add to smoothing factor
                useradd=float(choice)
            else :
                print("Error: Unknown choice, please retype correctly.")
        except(IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")
    #Calculate the flux inside the masked areas of continuum
    spec[:,1]=spec[:,1]/CleanCont   #Normalised spectrum
    flux=np.sum(spec[np.ma.getmaskarray(continuum[:,1]),1]-1)
    return spec,flux,CleanCont
         
def InterpolateRemoveFringes(InNspec,norm=True):
    """ Interactively lets user shift and mearge fringes from continuum on the spectral line to clean fringes on spectral line 
    Input: InNspec is 2d numpy array with first column wavelength and 2nd column flux
            norm is Boolean variable, set it to False, if continuum normalisation has to be done from this function
    Output: Fringe removed spectra , Continuum """ 
    Nspec=InNspec.copy() #Copying to protect input
    if not norm : Nspec,flux,continuum=ContinuumNormSpec(Nspec)  #Normalise first if norm is set to False
    else : continuum=np.ones(len(Nspec[:,1]))  # Continuum to return is just 1.
    FringeDic=dict()  #Dictionary to store fringe template arrays
    ShiftDic=dict()  #Dictionary to store shifts in each fringe templates
    while True:
        plt.close()
        L=len(Nspec[:,1])
        plt.plot(range(L),Nspec[:,1],marker='o',alpha=0.5)
        avgFringe=np.zeros(L)
        count=np.zeros(L)
        for i in FringeDic:   #Ploting all the fringe templates
            plt.plot(range(ShiftDic[i],ShiftDic[i]+len(FringeDic[i])),FringeDic[i],'-.',label=i)
            #Add to the combined fringe template to divide the original spectra
            avgFringe[ShiftDic[i]:ShiftDic[i]+len(FringeDic[i])]+=FringeDic[i]
            count[ShiftDic[i]:ShiftDic[i]+len(FringeDic[i])]+=1
        #Replace Zeros with 1 before division
        count[count==0]=1
        avgFringe/=count
        avgFringe[avgFringe==0]=1  #Replace all zeros with 1...
        #Plot the Fringe normalised spectrum
        plt.plot(range(L),Nspec[:,1]/avgFringe,'--',color='black',label='Result')
        plt.grid(True)
        plt.legend()
        plt.show(block=False)
        print('-'*25)
        print("Following example shows various input option available to you.")
        print(" n a 10 41  #Selects the region 10 to 40 as new fringe template with label 'a' ")
        print(" m a 40     #Moves the fringe template 'a' to location 40 ")
        print(" wq         #Exit by returning the present fringe removed profile ")        
        print(" q          #Discards everything and exit ")
        choice=raw_input("Enter your option : ").strip(' ')
        print('-'*25)
        try:
            if choice[0] == 'q': #Exiting discarding everything
                confirm=raw_input("Do you want to exit, Discarding all the fringe fitting? (Enter y to exit) : ").strip(' ')
                if confirm.lower() == 'y' :
                    print("Exiting, discarding everything. No changes made to spectral line")
                    break
            elif choice[0] == 'n' and len(choice.split())==4:  # Selecting new fringe template
                label=choice.split()[1]
                start=int(choice.split()[2])
                end=int(choice.split()[3])
                FringeDic[label]=Nspec[start:end,1]
                ShiftDic[label]=start
                print("Selected fringe template: %s of length %d, at location %d"%(label,len(FringeDic[label]),ShiftDic[label]))
            elif choice[0] == 'm' and len(choice.split())==3:  # moving fringe template
                label=choice.split()[1]
                print("Previous location of fringe template: %s was %d"%(label,ShiftDic[label]))
                start=int(choice.split()[2])
                ShiftDic[label]=start
                print("Current location of fringe template: %s is %d"%(label,ShiftDic[label]))
            elif choice[0:2] == 'wq' :  # Exit returning the fringe remove line profile
                Nspec[:,1]=Nspec[:,1]/avgFringe
                break
            else :
                print("Error: Unknown choice, please retype correctly.")
        except (IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")
    #Return the fringe removed (or not!) Nspec

    return Nspec,continuum
        

def RatioToScaleEQW(InSpec1,InSpec2):
    """ Interactively chooses the ratio of equivalent width of line in first spectra to second spectra.
    Input: Two 2d array of flux vectors.   
    IMP: For ratio to make sense, both the spectra should have been sampled at same wavelengths """
    Spec1=InSpec1.copy()
    Spec2=InSpec2.copy()
    print("Measuring the first spectra")
    Nspec1,flux1,continuum1=ContinuumNormSpec(Spec1)  #Normalise the spectras and estimate the EQW in arbitary units
    print("Measuring the first spectra")
    Nspec2,flux2,continuum2=ContinuumNormSpec(Spec2)  #Normalise the spectras and estimate the EQW in arbitary units
    print("EQW in arbitrary units of first spectra is %f and second spectra is %f"%(flux1,flux2))
    ORatio=flux1/flux2
    Ratio=ORatio
    print("Recommended Ratio of them : %f"%(Ratio)) 
    while True:
        plt.close()
        plt.plot(Nspec1[:,0],Nspec1[:,1])
        plt.plot(Nspec2[:,0],Nspec2[:,1])
        #Plot the ratio of the Second spectra by first spectra after scaling with the Ratio
        plt.plot(Nspec1[:,0],(((Nspec2[:,1]-1)*Ratio)+1)/Nspec1[:,1],'--',drawstyle='steps',color='black')
        plt.grid(True)
        plt.show(block=False)
        print('-'*10)
        print("Following example shows various input option available to you.")
        print(" 1.5 Change the ratio of equivalent with to 1.5 instead of original %f"%(ORatio))
        print(" q   Exit returning the finally decided ratio of eqw.")
        choice=raw_input("Enter your option : ").strip(' ')
        print('-'*10)
        try:
            if choice[0] == 'q': #Exiting returning the ratio
                break
            elif is_number(choice) : # Add to smoothing factor
                Ratio=float(choice)
            else :
                print("Error: Unknown choice, please retype correctly.")
        except(IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")
    #Return the finally accepted Ratio
    return Ratio

def FindConvolutionKernel(InSpec1,InSpec2,wlshift=0,Ratio=1):
    """ Interactively find the convolution kernel to convolve First spectra to match Second spectra 
    Input: Two 2d array of flux vectors. """
    Spec1=InSpec1.copy()  # Normalised high resolution star spectra
    Spec2=InSpec2.copy()
    #First we have to sample the second spectra exactly to same axis of first spectra
    tck2=interp.splrep(Spec2[:,0],Spec2[:,1],s=0)
    Spec2=Spec1.copy()  #Making the array same size of first spectra
    Spec2[:,1]=interp.splev(Spec1[:,0]-wlshift,tck2,der=0)
    #We have to normalise the continuum of Second spectra. (First is believed to be already normalised).
    Spec2,flux2,continuum2=ContinuumNormSpec(Spec2)
    #We also have to scale the flux of the line spectra to Ratio
    Spec1[:,1]=((Spec1[:,1]-1)*Ratio)+1
    #Ask user, whether they want to remove fringes.
    userinput=raw_input("Enter y if you want to remove fringes first, before deconvolveing: ").strip(' ')
    if userinput[0].lower()== 'y': #Run fringe cleaning algorithm
        print("Remove all fringes from the black dashed plot")
        Spec2,junkContinuum=InterpolateRemoveFringes(Spec2,norm=True)
    else : print("Okay, No fringe removal done")
    print("Starting Deconvolution...")
    #We can start with a kernel which is of the same shape of True spectral line.Not important
    FirstKernel=1/Spec1[:,1] -1  #Flipping and subtracting 1
    FirstKernel/=np.sum(FirstKernel) # Normalising the total flux inside the kernel
    Kernel=FirstKernel.copy()
    CSpec1=Spec1.copy()
    RLcount=0
    while True:
        CSpec1[:,1]=np.convolve(Spec1[:,1],Kernel,mode='same')
        plt.close()
        plt.figure(1)
        plt.subplot(211) #Top plot in the two layer subplot window
        plt.plot(Spec2[:,0],Spec2[:,1],label='Star')
        plt.plot(CSpec1[:,0],CSpec1[:,1],label='Convolved profile')
        plt.plot(Spec1[:,0],Spec2[:,1]/CSpec1[:,1],'--',drawstyle='steps',color='black', label='Ratio')
        plt.grid(True)
        plt.legend()
        plt.subplot(212) #And plot latest Kernel in the second subplot.
        plt.plot(Kernel,label='Kernel iter='+str(RLcount))
        plt.legend()
        plt.grid(True)
        plt.show(block=False)
        print('**'+'-'*20)
        print("No of Lucy iterations done so far = %d"%(RLcount))
        print("Following example shows various input option available to you.")
        print(" 5   #Run 5 more iterations of Richardson-Lucy deconvolution to generate a new kernel from present one.")
        print(" n   #Normalise the kernel.")
        print(" r   #Reset the kernel to the Initial one.")
        print(" q   #Exit returning the latest generated kernel for convolution.")
        choice=raw_input("Enter your option : ").strip(' ')
        print('**'+'-'*20)
        try:
            if choice[0] == 'q': #Exiting returning the present kernel
                break
            elif choice[0] == 'r': #Reset kernel to initial one. 
                print("Reseting the kernel to first one.")
                Kernel=FirstKernel.copy()
                RLcount=0
            elif choice[0] == 'n': #Normalsing the kernel. 
                print("Normalising the kernel to unit flux.")
                Kernel/=np.sum(Kernel)
            elif is_number(choice) : # Run RL more times
                RLtorun=int(choice)
                Kernel=RichardsonLucy(Spec2[:,1],Spec1[:,1],Kernel,RLtorun)
                RLcount+=RLtorun
            else :
                print("Error: Unknown choice, please retype correctly.")
        except (IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")
        
    return Kernel
        
              
def AskAndCleanStellarLines(Inspec):  
    """ Interactively helps user to fit a convolved function of standard star spectral line and remove it """
    spec=Inspec.copy() #Copying to protect input
    spec[:,1]/=np.median(spec[:,1])  #Scaling for neatness of plot.
    Star='A0'
    StdNormSpecdic=dict([('A0','Vega_R2000_WR9000to30000_RS10_Norm.txt')])
    StdContBreakdic=dict([('A0','Cont_break_Vega_R2000_WR9000to30000_RS10_Norm.txt')])
    StdNorm=np.loadtxt(StdNormSpecdic[Star])
#    StdContBreak=np.loadtxt(StdContBreakdic[Star])
    StdConvolved=False
    WLshift=0.0
    ratio=1
    while True:
        #First we generate a spline interpolation of True standard star spectra
        Std_tck=interp.splrep(StdNorm[:,0]+WLshift,((StdNorm[:,1]-1)*ratio)+1,s=0)
        #Now we generate the sampled values of std spectra for input spectra.
        Std4Spec=interp.splev(spec[:,0],Std_tck,der=0)
        plt.close()
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(spec[:,1],marker='.')    
        ax1.plot(Std4Spec)
        if StdConvolved : ax1.plot(spec[:,1]/Std4Spec,'--',drawstyle='steps',color='black')
        ax1.grid(True)
        ax2 = ax1.twiny()   # Adding the wavelength axis at top of the plot.
        ax2.plot(spec[:,0],spec[:,1],alpha=0)
        plt.show(block=False)
        print('-'*20)
        print("Following example shows various input option available to you.")
        print(" me 210 251  #Select region from 210 to 250 for matching the equivalent width of line")
        print(" fk 210 251  #Select region from 210 to 250 for finding convolution kernel and apply")
        print(" sh 20       #Shift stellar line profiles by +20 Angstrom")
        print(" wq          #Exit this section returning the cleaned Telluric line spectra")
        print(" q           #Exit this section discarding everything")
        choice=raw_input("Enter your option : ").strip(' ')
        print('-'*20)
        try:
            if choice[0] == 'q': #Exiting discarding everything
                confirm=raw_input("Do you want to exit, Discarding all the line fitting? (Enter y to exit) : ").strip(' ')
                if confirm.lower() == 'y' :
                    print("Exiting, discarding everything. No changes made to Telluric std spectra")
                    break
            elif choice[0] == 'm' and len(choice.split())==3: #Match equivalent width
                start=int(choice.split()[1])
                end=int(choice.split()[2])
                ratio=RatioToScaleEQW(Std4Spec[start:end],spec[start:end,:])
            elif choice[0] == 'f' and len(choice.split())==3: #Find the instrumental convolution kernel
                start=int(choice.split()[1])
                end=int(choice.split()[2])
                Wstart=spec[start,0]+WLshift
                Wend=spec[end,0]+WLshift
                SNstart=np.argmax(StdNorm[:,0] > Wstart)
                SNend=np.argmax(StdNorm[:,0] > Wend)
                kernel=FindConvolutionKernel(StdNorm[SNstart:SNend,:],spec[start:end,:],wlshift=WLshift,Ratio=ratio)
                #Applying the kernel
                confirm=raw_input("Do you want to apply the newly made kernel? (Enter y to apply) : ").strip(' ')
                if confirm.lower() == 'y' :
                    print("Convolving the line profiles with newly made kernal")
                    StdNorm[:,1]=np.convolve(StdNorm[:,1],kernel,mode='same')
                    StdConvolved=True
                else : print("Discarding the kernel generated for convolution.")
            elif choice[0] == 's': #Shifting of Telluric spectra
                print("Previous spectral shift was %f points"%(WLshift))
                WLshift+=float(choice.split()[1])
                print("New spectral shift is %f points"%(WLshift))
            elif choice[0:2] == 'wq' :  # Exit returning the Spectral line removed profile
                spec[:,1]=spec[:,1]/Std4Spec
                break
            else :
                print("Error: Unknown choice, please retype correctly.")
        except (IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")

    return spec

def TelluricCorrect(InSci,InStd):
    """ Interactively asks user to align, scale, and remove spectral lines to correct telluric lines and fringes.
    Input : Sci is the 2d numpy array with 1st column wavelength and 2nd column flux of Science target
    Input : Std is the 2d numpy array with 1st column wavelength and 2nd column flux of Standard star """
    Sci=InSci.copy() #Copying to protect input
    Std=InStd.copy()
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
                ScaleStd=1.0/np.median(Std[:,1])
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
                print("New spectral shift is %f Angstrom"%(StdWLshift))
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
