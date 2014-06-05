#!/usr/bin/env python
""" This script is to calculate the exposure times for up-the-ramp readout of TIRSPEC observation """
import numpy as np

def SlopeSigma(Count,N,readnoise=4.0,gain=5.0):
    """ Returns the error in slope of up-the-ramp line fitting
    Input :
        Count : The ADU count of flux per unit readout time. (slope)
        N  : Number of readout points used in slope fitting.
        readnoise : rms Readnoise sigma of a single readout in ADU.
        gain  : ADU to electron convertion factor
    Output :
        sigma : The expected sigma error in slope (Count)
    """
    EffReadNoiseVar=readnoise**2 *12.0/(N*(N-1)*(N+1))    #in ADU units
    PhotonNoiseVar= (Count*gain *6.0*(N**2 +1)/(5*N*(N**2 -1)))/gain**2  #dividing by gain square to get in ADU units
    sigma=np.sqrt(EffReadNoiseVar+PhotonNoiseVar)
    return sigma


def StarProfile(flux,fwhm,bkg=0,size=None):
    """ Returns a size x size square matrix filled with star profile intensity and also an array of same size with radius from center. It assumes a gaussian profile.
    Input :
       flux: Total flux from star
       fwhm: FWHM of star profile
       bkg (optional): background flux. default=0
       size (optional): size of the side square matrix to return. default=5*fwhm
    Output:
      Profile: Matrix populated with flux of star profile centered in the matrix.
      Radius: Radius of each pixel from the center."""

    if size is None : size=int(5*fwhm)
    sigma=fwhm/2.35482
    amplitude=flux/(2*np.pi*sigma*sigma)
    X,Y=np.indices((size,size),dtype=np.float)
    center=size/2.0
    Profile=amplitude*np.exp(-(((X-center)/sigma)**2 +((Y-center)/sigma)**2)/2) +bkg
    Radius=np.sqrt((X-center)**2 +(Y-center)**2)
    return Profile,Radius

def Signal2NoiseRatio(N,Profile,Radius,apperture,bkg,bkgsigma):
    """ Returns S/N ratio for various apperture size photometry 
    Input:
       N  : Number of readout points used in slope fitting.
       Profile: 2d matrix which has star profile in it.
       Radius : 2d matrix which has radius from center of profile in it.
       apperture: a tuple which contains list of appertures to use for summing flux and noise.
       bkg: the background estimate per pixel.
       bkgsigma: the sigma in background estimate per pixel.
    Output:
       SbyN : S/N ratios corresponding to input appertures
    """
    slopeSigma=SlopeSigma(Profile,N*np.ones(Profile.shape))
    Sigma=np.sqrt(slopeSigma**2 + bkgsigma**2)
    SbyN=[]
    for app in apperture:
        if 2*app > Profile.shape[0] : 
            print "Apperture {0} too large than profile.".format(app)
            continue
            
        signal=np.sum(Profile[Radius<app]-bkg)
        noise=np.sqrt(np.sum(Sigma[Radius<app]**2))
        SbyN.append(signal/noise)
    return tuple(SbyN)


def PhotometrySbyN(flux,fwhm,N,bkg,apperture,bkgsigma=None):
    """ Returns S/N ratio for various apperture size photometry for gaussian stars of given flux. 
    Input:
       flux: Total flux from star
       fwhm: FWHM of star profile
       N  : Number of readout points used in slope fitting.
       bkg: the background estimate per pixel.
       apperture: a tuple which contains list of appertures to use for summing flux and noise.
       bkgsigma (optional): the sigma in background estimate per pixel. default S/N correponding to bkg value.

    Output:
       SbyN : S/N ratios corresponding to input appertures
    """
    size=int(max(apperture)*2+1)
    Profile,Rad=StarProfile(flux,fwhm,bkg=bkg,size=size)
    if bkgsigma is None : bkgsigma=SlopeSigma(bkg,N)
    SbyN=Signal2NoiseRatio(N,Profile,Rad,apperture,bkg,bkgsigma)
    return SbyN


def Mag2Counts(Mag,Filter):
    """ Returns the ADU/sec counts corresponding to Magnitude and filter, based on zero point of the instrument 
    Input:
       Mag: Magnitude of source in Vega system
       Filter: Options are = 'J', 'H' or 'KS'
    
    Output:
       CountsPerSec: ADU/sec 
    """
    ZeroPoint=dict()
    # All Filter names Keys should be in uppercase
    ZeroPoint['J']=20.8
    ZeroPoint['H']=20.9
    ZeroPoint['KS']=20.15

    Filter=Filter.upper()
    try:
        CountsPerSec= 10**((ZeroPoint[Filter]-Mag)/2.5)
    except KeyError:
        print('ERROR: Unknown filter :'+Filter)
        print('Available filters for calculation are '+','.join(ZeroPoint.keys()))
        raise
    else:
        return CountsPerSec

def SbyN_of_StarPhotometry(StarMag,Filter,itime,seeing=1.6,SkyMag=None,apperture=None,bkgsigma=0,NDRtime=0.9):
    """ Returns the S/N ratio for a star in itime exposure through a filter
    Input:
       StarMag: Magnitude of star in Vega system
       Filter: Photometric Filter used for observation
       itime : Exposure time given
       seeing  : Seeing of the night in units of arcsec (default: 1.6 arcsec)
       SkyMag : Sky Brightness Magnitude in mag/arcsec^2 (default: typical estimate at HCT)
       apperture: Appertures used for photometry (default= 1,2 and 3 times seeing)
       bkgsigma : Error in estimate of background (default=0)
       NDRtime: Time required for one Non-distructive up-the-ramp readout.(deafult=0.9 sec)
    Output:
       SbyN  : Signal to Noise Ratio of photometry
    """

    Filter=Filter.upper()
    #HCT Sky brightness estimates on a good night in mag/arcsec^2
    SkyMagHCT=dict()
    SkyMagHCT['J']=16.14
    SkyMagHCT['H']=13.7
    SkyMagHCT['KS']=13.35


    NoOfFrames=itime/(NDRtime*1.0)
    StarCountsPerNDR=Mag2Counts(StarMag,Filter)*NDRtime

    if SkyMag is None: SkyMag=SkyMagHCT[Filter]

    SkyCountsPerNDRPerPixel=Mag2Counts(SkyMag,Filter)*NDRtime*(0.3*0.3)
    fwhm=seeing/0.3
    if apperture is None : apperture=(1*fwhm,2*fwhm,3*fwhm)

    SbyN=PhotometrySbyN(StarCountsPerNDR,fwhm,NoOfFrames,SkyCountsPerNDRPerPixel,apperture,bkgsigma=bkgsigma)

    return SbyN

    
if __name__ == "__main__":
    import matplotlib
    matplotlib.use('TkAgg')
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    # implement the default mpl key bindings
    from matplotlib.backend_bases import key_press_handler
    from matplotlib.figure import Figure

    import Tkinter as Tk

    ExpRange=range(4,31)

    def RedrawPlot():
        ax.clear()
        seeing=float(SeeingV.get())
        StarMag=float(StarMagV.get())
        Filter=FilterV.get()
        N=int(NoFV.get())
        appertures=(.5*seeing/0.3,1*seeing/0.3,2*seeing/0.3,3*seeing/0.3)
        plotedlines=ax.plot(ExpRange,CalculateSbyForItimes(ExpRange,N,StarMag,Filter,seeing,appertures))
        ax.legend(plotedlines,('0.5 FWHM','1 FWHM','2 FWHM','3 FWHM'))
        ax.set_title('Net S/N versus exposure per frame in {0} filter'.format(Filter))
        ax.set_xlabel('Itime per frame (sec)')
        ax.set_ylabel('Net S/N')
        ax.grid()
        fig.canvas.draw_idle()

    def CalculateSbyForItimes(ITime,N,StarMag,Filter,seeing,appertures):
        """ Returns S/N for input Itime array """
        Nsqrt=np.sqrt(N)
        SbyNList=[]
        for itime in ITime :
            SbyNList.append(tuple([SbyN*Nsqrt for SbyN in SbyN_of_StarPhotometry(StarMag,Filter,itime,seeing=seeing,apperture=appertures)]))

        return SbyNList
        

    root = Tk.Tk()
    root.wm_title("TIRSPEC Exposure Time Calculator")
    
    FilterLabel=Tk.Label(root, text='Choose Filter -->')
    FilterLabel.grid(row=1,column=0,columnspan=2)

    FilterV = Tk.StringVar()
    FiltButton1 = Tk.Radiobutton(root, text="J", variable=FilterV, value='J')
    FiltButton1.grid(row=1,column=2)
    FiltButton2 = Tk.Radiobutton(root, text="H", variable=FilterV, value='H')
    FiltButton2.grid(row=1,column=3)
    FiltButton3 = Tk.Radiobutton(root, text="Ks", variable=FilterV, value='KS')
    FiltButton3.grid(row=1,column=4)
    FilterV.set('J')

    StarMagLabel=Tk.Label(root, text='   Mag:')
    StarMagLabel.grid(row=2,column=0)
    StarMagV = Tk.StringVar()
    MagInpBox= Tk.Entry(root,textvariable=StarMagV)
    MagInpBox.grid(row=2,column=1)
    StarMagV.set('0')

    SeeingLabel=Tk.Label(root, text=' Seeing("):')
    SeeingLabel.grid(row=2,column=2)
    SeeingV = Tk.StringVar()
    SeeingInpBox= Tk.Entry(root,textvariable=SeeingV)
    SeeingInpBox.grid(row=2,column=3)
    SeeingV.set('1.6')

    NoFLabel=Tk.Label(root, text=' No# Frames:')
    NoFLabel.grid(row=2,column=4)
    NoFV = Tk.StringVar()
    NoFInpBox= Tk.Entry(root,textvariable=NoFV)
    NoFInpBox.grid(row=2,column=5)
    NoFV.set('15')

    fig = Figure(figsize=(8,5), dpi=80)
    ax = fig.add_subplot(111)
    
    ReplotButton = Tk.Button(root, text='Refresh', command=RedrawPlot)
    ReplotButton.grid(row=3,column=1,columnspan=6)
    # a tk.DrawingArea
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.show()
    canvas.get_tk_widget().grid(row=4,column=0, columnspan=6)

    Tk.mainloop()
