#!/usr/bin/env python
#This script will help in telluric correction of spectra with standard star spectra interactively.
import readline
import sys
import numpy as np
import numpy.ma
import scipy.interpolate as interp
import scipy.constants
import scipy.optimize
import matplotlib.pyplot as plt
from astropy.io import fits  #Needed only to load fits spectra

def RichardsonLucy(Signal,PSF,IniKernel,niter):
    """ By Richardson-Lucy Algorithm, retrun the deconvolution of Signal using input PSF.
    Input: Signal : 1d numpy array, which is to be deconvolved
           PSF    : 1d numpy array, which is the PSF to be used for deconvolution
           IniKernel: 1d numpy array, which is the initial guess of original signal
           niter    : integer, Number of iterations of RL deconvolution steps to be run """
    signal=Signal.copy()   #Copying to protect input
    psf=PSF/np.sum(PSF)   #Normalised PSF
    inikernel=IniKernel.copy()
    # Algorithm and syntax based on http://en.wikipedia.org/wiki/Richardson-Lucy_deconvolution
    psf_hat=psf[::-1] #psf hat, following wikipedia's syntax
    for i in range(niter):
        est_conv = np.convolve(inikernel,psf,mode='same')
        relative_blur = signal/est_conv
        error_est = np.convolve(relative_blur,psf_hat,mode='same') 
        inikernel = inikernel* error_est
        
    return inikernel



#------------------------------------1D gaussian fitting functions...
def moments1D(inpData):
    """ Returns the (amplitude, center,sigma, bkgC, bkgSlope) estimated from moments in the 1d input array inpData """

    bkgC=inpData[0]  #Taking first point as intercept and fitting a straght line for slope
    bkgSlope=(inpData[-1]-inpData[0])/(len(inpData)-1.0)
    
    Data=np.ma.masked_less(inpData-(bkgC+bkgSlope*np.arange(len(inpData))),0)   #Removing the background for calculating moments of pure gaussian
    #We also masked any negative values before measuring moments

    amplitude=Data.max()

    total=float(Data.sum())
    Xcoords=np.arange(Data.shape[0])

    center=(Xcoords*Data).sum()/total
 
    sigma=np.sqrt(np.ma.sum((Data*(Xcoords-center)**2))/total)

    return amplitude,center,sigma,bkgC,bkgSlope

def Gaussian1D(amplitude, center,sigma,bkgC,bkgSlope):
    """ Returns a 1D Gaussian function with input parameters. """
    Xc=center  #Center
    #Now lets define the 1D gaussian function
    def Gauss1D(x) :
        """ Returns the values of the defined 1d gaussian at x """
        return amplitude*np.exp(-(((x-Xc)/sigma)**2)/2) +bkgC+bkgSlope*x

    return Gauss1D

def FitGauss1D(Data,ip=None):
    """ Fits 1D gaussian to Data with optional Initial conditions ip=(amplitude, center, sigma, bkgC, bkgSlope)
    Example: 
    >>> X=np.arange(40,dtype=np.float)
    >>> Data=np.exp(-(((X-25)/5)**2)/2) +1+X*0.5
    >>> FitGauss1D(Data)
    (array([  1. ,  25. ,   5. ,   1. ,   0.5]), 2)
     """
    if ip is None:   #Estimate the initial parameters from moments 
        ip=moments1D(Data)
    
    def errfun(ip):
        return np.ravel(Gaussian1D(*ip)(*np.indices(Data.shape)) - Data)

    p, success = scipy.optimize.leastsq(errfun, ip)

    return p,success

#------------------------------------------------------------#

def FitGaussianLineProfile(XYDataToFit,absorption=True,displayfit=True):
    """ Fits a 1 D gaussian to XYDataToFit numpy array. And returns the best fitted gaussian
    Input :
         XYDataToFit : 2 coulumn numpy array with X and Y data.
         absorption : True to fit absorption line gaussian, and False to fit emission line gaussian 
         displayfit : True to finally display the best fit plot, False to not show any fit at the end.
    Output:
         BestFitPureGaussian : Best Fit pure gaussian of on the input data without any background
         p                   : Full set of parameters of best fit Gaussian1D
    """
    Ydata=XYDataToFit[:,1].copy()
    if absorption : Ydata*=-1   # To make it emission
    
    p,succ=FitGauss1D(Ydata,ip=None)

    if absorption : 
        p[0]*=-1  #Amplitude
        p[3]*=-1  #Background
        p[4]*=-1  #Background's Slope

    centerW=XYDataToFit[p[1],0]
    sigma=p[2]*(XYDataToFit[1,0]-XYDataToFit[0,0])
    print('Amplitude:{0}, Center:{1}, Sigma:{2}, Bkg:{3}, Bkg Slope:{4}'.format(p[0],centerW,sigma,p[3],p[4]))
    pureGauss=[p[0],p[1],[2],0,0]  # ie. with Bkg and Bkg Slope set to Zero
    BestFitGaussian=Gaussian1D(*p)(np.arange(len(Ydata)))
    BestFitPureGaussian=Gaussian1D(*pureGauss)(np.arange(len(Ydata)))
    if displayfit:
        plt.close() #Closing any previously open window
        fig = plt.figure()
        ax=fig.add_subplot(1,1,1)
        ax.plot(XYDataToFit[:,0],XYDataToFit[:,1],marker='.',alpha=0.5) 
        ax.plot(XYDataToFit[:,0],BestFitGaussian,color='k') 
        ax.plot(XYDataToFit[:,0],XYDataToFit[:,1]-BestFitPureGaussian,color='red') 
        plt.show()

    return BestFitPureGaussian,p

def BlackBodyFlux(InWavelengths,Temp):
    """ Returns an array of median scaled flux from blackbody corresponding to the input list of wavelengths, and given Temperature 
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
    plt.close() #Closing the window, just checking whether this removes the error while calling m
    fig = plt.figure()
    ax=fig.add_subplot(1,1,1)
    fullline,=ax.plot(spec[:,1],marker='.',alpha=0.5)  #Saving the input full plot
    plt.close()  #Close any previous windows and clear fig
    while True:
        if not plt.get_fignums():  #If window is not open
            fig = plt.figure()
            ax=fig.add_subplot(1,1,1)
            fig.canvas.set_window_title('Masking')
            ax.lines.append(fullline)  #Adding the original line.
            ax.grid(True)
            line,=ax.plot(spec[:,1],marker='.')
            plt.show(block=False)
        
        ax.lines.remove(line)
        line,=ax.plot(spec[:,1],marker='.')
        plt.draw()
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
    plt.close()  #close window 
    return spec

def AdjustAlphaPlots(ax,linelist,PlotHistLength,OnlyOneLabel=False):
    """ In plot axis ax, this subroutine adjusts the transparency (alpha) of lines in linelist in increasing order.
    It will also reduce the linelist to size of PlotHistLength by removing from beginning of list.
    Those removed lines will also be removed from ax plot.
    If OnlyOneLabel = True, then the labels of lines except the last one will be removed."""
    if len(linelist) > PlotHistLength : 
        for line in linelist[:-1*PlotHistLength] : 
            ax.lines.remove(line)  #Removing older plots
        linelist=linelist[-1*PlotHistLength:]   #Removing old lines from list also
    #Adjusting the alpha of lines in decreasing order.
    Nooflines=len(linelist)
    for i,line in enumerate(linelist) :  
        line.set_alpha((1.0+i)/Nooflines)
    
    if OnlyOneLabel :  #If all the labels except last one is asked to be removed
        for line in linelist[:-1] : line.set_label('')
    #Return the modified plots and line list
    return ax,linelist
    

def ContinuumNormSpec(Inspec):
    """ Returns the continuum normalised spectra by interactively masking lines
    Input: spec : a 2d numpy array, with first column wavelength and second column flux 
    Output: Continuum normalise spec, flux of line in arbitrary unit and Fitted continuum """
    spec=Inspec.copy() #Copying to protect input
    continuum=np.ma.array(spec)
    useradd=0
    scale=np.median(spec[:,1])
    PlotHistLength=3   #Number of previous history plots to show.
    Domask=True
    plt.close()
    while True:
        if Domask : 
            print("Mask away the regions which are not part of continuum")
            continuum=AskAndMaskSpectra(continuum)
            Domask=False
        l=np.ma.count(continuum,axis=0)[1]  #Number of points to fit continuum
        sm=(l+np.sqrt(2*l))*6+useradd   #Smoothing factor for fitting.
        cont_tck=interp.splrep(np.ma.compressed(continuum[:,0]),np.ma.compressed(continuum[:,1]),s=sm)
        CleanCont=interp.splev(spec[:,0],cont_tck,der=0)

        if not plt.get_fignums():  #If window is not already open
            fig = plt.figure()
            ax=fig.add_subplot(1,1,1)
            fig.canvas.set_window_title('Normalising Continuum')
            ax.plot(spec[:,0],spec[:,1]/scale)
            ax.grid(True)
            plt.show(block=False)
            ContPlotlist=[]
            NspecPlotlist=[]

        NewContPlotline,=ax.plot(spec[:,0],CleanCont/scale)
        NewNspecPlotline,=ax.plot(spec[:,0],spec[:,1]/CleanCont)
        ContPlotlist.append(NewContPlotline)
        NspecPlotlist.append(NewNspecPlotline)
        # Adjust the transparencies and remove old plots
        ax,ContPlotlist=AdjustAlphaPlots(ax,ContPlotlist,PlotHistLength)
        ax,NspecPlotlist=AdjustAlphaPlots(ax,NspecPlotlist,PlotHistLength)

        plt.draw()
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

    plt.close()
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
    plt.close()
    PlotHistLength=3   #Number of previous history plots to show.
    L=len(Nspec[:,1])
    while True:
        if not plt.get_fignums():  #If window is not already open
            fig = plt.figure()
            ax=fig.add_subplot(1,1,1)
            fig.canvas.set_window_title('Fringe Removal')
            ax.plot(range(L),Nspec[:,1],marker='o',alpha=0.5)
            ax.grid(True)
            plt.show(block=False)
            ResultPlotlist=[]

        avgFringe=np.zeros(L)
        count=np.zeros(L)
        for i in FringeDic:   #Ploting all the fringe templates
            ax.lines = [ line for line in ax.lines if line.get_label() != i ]  #Removing any previous plots of this label
            ax.plot(range(ShiftDic[i],ShiftDic[i]+len(FringeDic[i])),FringeDic[i],'-.',label=i)
            #Add to the combined fringe template to divide the original spectra
            avgFringe[ShiftDic[i]:ShiftDic[i]+len(FringeDic[i])]+=FringeDic[i]
            count[ShiftDic[i]:ShiftDic[i]+len(FringeDic[i])]+=1
        #Replace Zeros with 1 before division
        count[count==0]=1
        avgFringe/=count
        avgFringe[avgFringe==0]=1  #Replace all zeros with 1...
        #Plot the Fringe normalised spectrum
        NewResultline,=ax.plot(range(L),Nspec[:,1]/avgFringe,'--',color='black',label='Result')
        ResultPlotlist.append(NewResultline)
        ax,ResultPlotlist=AdjustAlphaPlots(ax,ResultPlotlist,PlotHistLength,OnlyOneLabel=True)
        ax.legend()
        plt.draw()
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
    plt.close()
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
    plt.close()
    PlotHistLength=3   #Number of previous history plots to show.
    while True:
        if not plt.get_fignums():  #If window is not already open
            fig = plt.figure()
            ax=fig.add_subplot(1,1,1)
            fig.canvas.set_window_title('EQW Ratio')
            ax.plot(Nspec1[:,0],Nspec1[:,1])
            ax.plot(Nspec2[:,0],Nspec2[:,1])
            ax.grid(True)
            plt.show(block=False)
            RatioPlotlist=[]

        #Plot the ratio of the Second spectra by first spectra after scaling with the Ratio
        NewRatioline,=ax.plot(Nspec1[:,0],(((Nspec2[:,1]-1)*Ratio)+1)/Nspec1[:,1],'--',drawstyle='steps',color='black')
        RatioPlotlist.append(NewRatioline)
        ax,RatioPlotlist=AdjustAlphaPlots(ax,RatioPlotlist,PlotHistLength)
        plt.draw()
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
    plt.close()
    return Ratio

def FindConvolutionKernel(InSpec1,InSpec2,wlshift=0,Ratio=1):
    """ Interactively find the convolution kernel to convolve First spectra to match Second spectra 
    Input: Two 2d array of flux vectors. """
    Spec1=InSpec1.copy()  # Normalised high resolution star spectra
    Spec2=InSpec2.copy()
    #First we have to sample the second spectra exactly to same axis of first spectra
    tck2=interp.splrep(Spec2[:,0],Spec2[:,1],s=0)
    Spec2=Spec1.copy()  #Making the array same size of first spectra
    Spec2[:,1]=interp.splev(Spec1[:,0]+wlshift,tck2,der=0)
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
    plt.close()
    PlotHistLength=3   #Number of previous history plots to show.
    while True:
        CSpec1[:,1]=np.convolve(Spec1[:,1],Kernel,mode='same')
        if not plt.get_fignums():  #If window is not already open
            fig = plt.figure()
            fig.canvas.set_window_title('Convolution Kernel')
            ax1 = fig.add_subplot(211) #Top plot in the two layer subplot window
            ax1.plot(Spec2[:,0],Spec2[:,1],label='Star')
            Convplotline,=ax1.plot(CSpec1[:,0],CSpec1[:,1],label='Convolved profile')
            ax1.grid(True)
            ax1.legend()
            ax2= fig.add_subplot(212) #For ploting latest Kernel in the second subplot.            
            ax2.grid(True)
            plt.show(block=False)
            RatioPlotlist=[]
            KernelPlotlist=[]

        ax1.lines.remove(Convplotline)  #Removing previous convoluted line plot
        Convplotline,=ax1.plot(CSpec1[:,0],CSpec1[:,1],label='Convolved profile')  #New one
        NewRatioLine,=ax1.plot(Spec1[:,0],Spec2[:,1]/CSpec1[:,1],'--',drawstyle='steps',color='black', label='Ratio')
        RatioPlotlist.append(NewRatioLine)
        ax1,RatioPlotlist=AdjustAlphaPlots(ax1,RatioPlotlist,PlotHistLength,OnlyOneLabel=True)
        ax1.legend()
        #And plot latest Kernel in the second subplot.
        NewKernelline,=ax2.plot(Kernel,label='Kernel iter='+str(RLcount))
        KernelPlotlist.append(NewKernelline)
        ax2,KernelPlotlist=AdjustAlphaPlots(ax2,KernelPlotlist,PlotHistLength)
        ax2.legend()
        plt.draw()
        print('**'+'-'*20)
        print("No of Lucy iterations done so far = %d"%(RLcount))
        print("Following example shows various input option available to you.")
        print(" 5   #Run 5 more iterations of Richardson-Lucy deconvolution to generate a new kernel from present one.")
#        print(" n   #Normalise the kernel.")
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
            # elif choice[0] == 'n': #Normalsing the kernel. 
            #     print("Normalising the kernel to unit flux.")
            #     Kernel/=np.sum(Kernel)
            elif is_number(choice) : # Run RL more times
                RLtorun=int(choice)
                Kernel=RichardsonLucy(Spec2[:,1],Spec1[:,1],Kernel,RLtorun)
                Kernel/=np.sum(Kernel) # "Normalising the kernel to unit flux."                 
                RLcount+=RLtorun
            else :
                print("Error: Unknown choice, please retype correctly.")
        except (IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")
    
    plt.close()
    return Kernel
        
              
def AskAndCleanStellarLines(Inspec):  
    """ Interactively helps user to fit a convolved function of standard star spectral line and remove it """
    spec=Inspec.copy() #Copying to protect input
    spec[:,1]/=np.median(spec[:,1])  #Scaling for neatness of plot.
    Star='A0'
    StdNormSpecdic=dict([('A0','Vega_R2000_WR9000to30000_RS10_Norm.txt')])
    StdContBreakdic=dict([('A0','Cont_break_Vega_R2000_WR9000to30000_RS10_Norm.txt')])
    try :
        StdNorm=np.loadtxt(StdNormSpecdic[Star])
#    StdContBreak=np.loadtxt(StdContBreakdic[Star])
    except IOError as e:
        print(e)
        print('File : '+StdNormSpecdic[Star]+' not found in present directory.')
        print('Please copy that file to present directory and proceed')
        sys.exit(1)

    StdConvolved=False
    GaussSubtracted=False
    FittedGaussians=dict()
    WLshift=0.0
    ratio=1
    plt.close()
    PlotHistLength=3   #Number of previous history plots to show.
    while True:
        #First we generate a spline interpolation of True standard star spectra
        Std_tck=interp.splrep(StdNorm[:,0]+WLshift,((StdNorm[:,1]-1)*ratio)+1,s=0)
        #Now we generate the sampled values of std spectra for input spectra.
        Std4Spec=interp.splev(spec[:,0],Std_tck,der=0)
        if not plt.get_fignums():  #If window is not already open
            fig = plt.figure()
            ax1=fig.add_subplot(1,1,1)
            fig.canvas.set_window_title('Spectral Line Removal')
            ax1.plot(spec[:,1],marker='.')
            stdspecplotline,=ax1.plot(Std4Spec)
            ax1.grid(True)
            ax2 = ax1.twiny()   # Adding the wavelength axis at top of the plot.
            ax2.plot(spec[:,0],spec[:,1],alpha=0)
            plt.show(block=False)
            RatioPlotlist=[]

        ax1.lines.remove(stdspecplotline)  #Removing previous spectral line plot
        stdspecplotline,=ax1.plot(Std4Spec)
        if StdConvolved : 
            NewRatioLine,=ax1.plot(spec[:,1]/Std4Spec,'--',drawstyle='steps',color='black')
            RatioPlotlist.append(NewRatioLine)
            ax1,RatioPlotlist=AdjustAlphaPlots(ax1,RatioPlotlist,PlotHistLength)
        elif GaussSubtracted :
            TempSpec=spec.copy()
            for x,y in FittedGaussians.keys(): TempSpec[x:y,1]-=FittedGaussians[(x,y)]
            NewRatioLine,=ax1.plot(TempSpec[:,1],'--',drawstyle='steps',color='black')
            RatioPlotlist.append(NewRatioLine)
            ax1,RatioPlotlist=AdjustAlphaPlots(ax1,RatioPlotlist,PlotHistLength)
            
        plt.draw()
        print('-'*20)
        print("Following example shows various input option available to you.")
        print("Note: fk and gs are mutually exclusive options, you can chose only one out of them")
        print(" me 210 251  #Select region from 210 to 250 for matching the equivalent width of line")
        print(" fk 210 251  #Select region from 210 to 250 for finding convolution kernel and apply")
        print(" sh 20       #Shift stellar line profiles by +20 Angstrom")
        print(" --OR--")
        print(" gs 210 251  #Select region from 210 to 250 for fitting Gaussian profile and subtract")
        print(" gs H        #Automatically fit all Hydrogen lines with Gaussian profiles and subtract")
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
                if not GaussSubtracted: 
                    start=int(choice.split()[1])
                    end=int(choice.split()[2])
                    Wstart=spec[start,0]-WLshift
                    Wend=spec[end,0]-WLshift
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
                else:
                    print('You have already choosen to do gaussian subtraction, so this step is not allowed. \n However you can quit this section and return back to choose this process')
            elif choice[0] == 'g' and len(choice.split())==3: #Gaussian fitting of absorption profiles
                if not StdConvolved:
                    start=int(choice.split()[1])
                    end=int(choice.split()[2])
                    Gaussianline,params=FitGaussianLineProfile(spec[start:end,:],absorption=True)
                    confirm=raw_input("Do you want to subtract the fitted gaussian? (Enter y to subtract) : ").strip(' ')
                    if confirm.lower() == 'y' :
                        print("Subtracting the best fit gaussian line profile.")
                        FittedGaussians[(start,end)]=Gaussianline   #Storing the line for subtracting in future
                        GaussSubtracted=True
                    else : print("Discarding the fitted gaussian line.")
                else:
                    print('You have already choosen to do star template division, so this step is not allowed. \n However you can quit this section and return back to choose this process')

            elif choice[0] == 'g' and choice.split()[1].upper()=='H': #Gaussian fitting of all hydrogen lines
                print('Not yet implemented..: Please use "gs start end" for each line.')
                pass

            elif choice[0] == 's': #Shifting of Telluric spectra
                print("Previous spectral shift was %f points"%(WLshift))
                WLshift+=float(choice.split()[1])
                print("New spectral shift is %f points"%(WLshift))
            elif choice[0:2] == 'wq' :  # Exit returning the Spectral line removed profile
                if StdConvolved:
                    spec[:,1]=spec[:,1]/Std4Spec
                elif GaussSubtracted:
                    for x,y in FittedGaussians.keys(): spec[x:y,1]-=FittedGaussians[(x,y)]
                break
            else :
                print("Error: Unknown choice, please retype correctly.")
        except (IndexError,ValueError):
            print("Error: Unknown choice, please retype correctly..")

    plt.close()
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
    plt.close() #Close any previously open plot window
    PlotHistLength=3   #Number of previous history plots to show.
    while True:
        #First we generate a spline interpolation of standard star spectra
        Std_tck=interp.splrep(np.ma.compressed(Std[:,0])+StdWLshift,np.ma.compressed(Std[:,1])*ScaleStd,s=0)
        #Now we generate the binned values of std spectra for science star.
        Std4Sci=interp.splev(Sci[:,0],Std_tck,der=0)
        if not plt.get_fignums():  #If window is not already open
            fig = plt.figure()
            fig.canvas.set_window_title('Telluric Line Removal')
            ax1 = fig.add_subplot(211) #Top plot in the two layer subplot window
            ax1.plot(Sci[:,0],Sci[:,1]*ScaleSci,label="Science")
            Stdplotline,=ax1.plot(Sci[:,0],Std4Sci,label="Telluric Std")
            ax1.grid(True)
            ax1.legend()
            ax2= fig.add_subplot(212) #For ploting latest ratio in second subplot.            
            ax2.grid(True)
            plt.show(block=False)
            RatioPlotlist=[]


        ax1.lines.remove(Stdplotline)  #Removing previous standard star line plot        
        Stdplotline,=ax1.plot(Sci[:,0],Std4Sci,label="Telluric Std")
        ax1.legend()
        # Now generate the divided spectra
        CorrectedSci=(Sci[:,1]*ScaleSci)/Std4Sci

        NewRatioLine,=ax2.plot(Sci[:,0],CorrectedSci,label=r'$\frac{Science}{Telluric}$')
        RatioPlotlist.append(NewRatioLine)
        ax2,RatioPlotlist=AdjustAlphaPlots(ax2,RatioPlotlist,PlotHistLength,OnlyOneLabel=True)
        ax2.legend()
        #Display
        plt.draw()
        
        print('-'*40)
        print("Following examples shows the options to chose :")
        print(" sh 23     #To shift the Telluric spectra by +23 Angstrom")
        print(" sc 1.1    #To scale the Telluric spectra by factor of 1.1")
        print(" m         #To start a window on masking regions in spectra of Telluric Std spectra")
        print(" f         #To start a window on fitting and removing stellar lines in Telluric Std spectra")
        print(" wq        #To exit, returning the latest telluric spectra and corrected science target to save")
        print(" q         #Exit and close, discarding everything")
        choice=raw_input("Enter your choice : ").strip(' ')
        print('-'*40)
        try:
            if choice[0] == 'q': #Exiting discarding everything
                confirm=raw_input("Do you want to exit, Discarding everything done so far? (Enter y to exit) : ").strip(' ')
                if confirm.lower() == 'y' :
                    print("Exiting, discarding everything.")
                    sys.exit(0)
            if choice[0] == 'w': #Exiting to save the new spectra
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
    plt.close()
    #Return the cleaned spectra
    return CorrectedSci,Sci[:,0],Std
            

def LoadFitsSpectra(filename):
    """Load and return the wavelength calibrated input fits spectra as 2 column numpy array. """
    fitsfile=fits.open(filename)
    flux=fitsfile[0].data

    try :
        ref_pixel = fitsfile[0].header['CRPIX1']
        coord_ref_pixel = fitsfile[0].header['CRVAL1']
        wave_per_pixel = fitsfile[0].header['CDELT1']
    except KeyError as e :
        print('Error: Missing keywords in fits header to do wavelength calibration')
        print(e)
        print('You might have entered wrong file name. Hence I am raising IOError')
        print('Enter the fits file name which is wavelength calibrated.')
        raise IOError
    else:
        w_start=coord_ref_pixel - ((ref_pixel-1) * wave_per_pixel)  #Starting wavelength
        Wavelengths=w_start+np.arange(len(flux))*wave_per_pixel

        return np.vstack((Wavelengths,flux)).T
    

def LoadAsciiSpectra(filename,Skipascii=69):
    """Load and return the ascii spectra as 2 column numpy array. 
    Skipasciii is an integer which tells the number of initial lines to skip while reading the data from file """
    while True:  #Keep asking number of lines to skip in the ascii if error comes
        try:
            Data=np.loadtxt(filename,skiprows=Skipascii)
            break
        except ValueError as e:
            print("Non float values in file? ")
            print(e)
            print("Present number of lines to skip in ascii file is %d "%Skipascii) 
            Skipascii=int(raw_input("Enter number of initial lines to skip in file:").strip(' '))

    return Data
    

def AskAndLoadData(toask="Enter Filename : ",Skipascii=69):
    """ Asks user to enter filename and return the table loaded as numpy array. 
    User can input, a fits file, ascii spectra or even a numpy .npy array.
    toask is string which contains the prompt to ask. 
    Skipasciii is an integer which tells the number of initial lines of header to skip, incase user inputs an ascii spectra file instead of a fits data file """

    while True:  #Keep asking till a file is properly loaded
        try :
            filename=raw_input(toask).strip(' ')
            if filename[-5:] == '.fits' : #User has inputed a fits file
                Data=LoadFitsSpectra(filename)
            elif filename[-4:] == '.npy' : #User has inputed a numpy .npy array file
                Data=np.load(filename)
            else:   #We assume it is an ascii spectra
                Data=LoadAsciiSpectra(filename,skiprows=Skipascii)
            break
        except IOError :
            print("Error: Cannot find the file %s, Please give the correct file name."%(filename))

    return Data,filename

def is_number(s):   # A funtion to check wheter string s is a number or not.
    try:
        float(s)
        return True
    except ValueError:
        return False


def main():
    """ Asks user for the standard star spectra and science target spectra files and start telluric correction """
    StdStarRaw,stdfname=AskAndLoadData(toask="Enter the name of standard star spectrum file : ",Skipascii=69)
    ScienceStarRaw,scifname=AskAndLoadData(toask="Enter the name of science star spectrum file : ",Skipascii=69)
    CorrSci,CorrSciWL,Telluric=TelluricCorrect(ScienceStarRaw,StdStarRaw)
    BBtemp='SkipCC'
    BBtemp=raw_input("Enter Temperature of Std star (Kelvin) to correct continuum (to skip press enter): ")
    if is_number(BBtemp): #User entered the temperature of Std star
        T=float(BBtemp)
        CorrSci=CorrSci*BlackBodyFlux(CorrSciWL,T)
    else:
        print("Not doing continuum correction using BB curve.")
    plt.close()
    fig=plt.figure()
    fig.canvas.set_window_title('Final Spectrum')
    ax = fig.add_subplot(1,1,1)
    ax.plot(CorrSciWL/10000.0,CorrSci,drawstyle='steps',color='black')
    ax.set_xlabel(r'Wavelength ($\mu m$)')
    ax.set_ylabel('Normalised counts')
    print("Final spectrum will be saved after window is closed.")
    plt.show()
    #Save the generated spectras values as numpy arrays
    fnameSuffix=raw_input("Enter any suffix, you want to add to the output spectra file names: ")
    scifname=scifname.split('.')[0]+'_'+fnameSuffix
    stdfname=stdfname.split('.')[0]+'_'+fnameSuffix
    np.save(scifname+'_Y_',np.array(np.ma.compressed(CorrSci)))  #Adding extra np.array to fix the bug in numpy till a fix 
    np.save(scifname+'_X_',np.array(np.ma.compressed(CorrSciWL)))
    np.save(stdfname+'_XY_',np.array(np.ma.compress_rows(Telluric)))
    print('The spectra saved by the following names')
    print(scifname+'_Y_.npy',scifname+'_X_.npy',stdfname+'_XY_.npy')
    print('Enjoy!!_________________________indiajoe')

if __name__ == "__main__":
    main()
