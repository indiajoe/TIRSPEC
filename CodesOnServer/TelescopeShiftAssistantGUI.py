#!/usr/bin/env python
""" This is a telescope shift assistant to be used while observing with TIRSPEC 
If you need to any help, feel free to contact me at indiajoe@gmail.com
-----------------------Enjoy Observing!! -  J.P.Ninan"""
import Tkinter
import ttk
import math

class AGtab(ttk.Frame):
    """ This Tab is the interface for calculating shifts to be given in Autoguider while tracking"""
    def __init__(self,parent):
        ttk.Frame.__init__(self,parent)
        self.parent=parent
        self.title="Auto Guider Control"
        self.grid()
        
        # The Top Clear Button
        self.ClearButton = Tkinter.Button(self,text=u"Clear All",command=self.OnClearButtonClick)
        self.ClearButton.grid(column=0,row=0,columnspan=2)

        #Current TIRSPEC's X Y position entry boxes
        self.CurrentXYEntryFrame = ttk.Labelframe(self, text='Current X Y')
        self.CurrentXYEntryFrame.grid(column=0,row=1,columnspan=2)
        self.InpCurrentX = Tkinter.StringVar()
        self.InpCurrentY = Tkinter.StringVar()
        self.CurrXentrybox = Tkinter.Entry(self.CurrentXYEntryFrame,textvariable=self.InpCurrentX)
        self.CurrYentrybox = Tkinter.Entry(self.CurrentXYEntryFrame,textvariable=self.InpCurrentY)
        self.CurrXentrybox.grid(column=0,row=0,sticky='EW')
        self.CurrYentrybox.grid(column=1,row=0,sticky='EW')
        self.CurrXentrybox.selection_range(0, Tkinter.END)        
        self.CurrYentrybox.selection_range(0, Tkinter.END)        


        #Target TIRSPEC's X Y position entry boxes
        self.TargetXYEntryFrame = ttk.Labelframe(self, text='Target X Y')
        self.TargetXYEntryFrame.grid(column=0,row=2,columnspan=2)
        self.InpTargetX = Tkinter.StringVar()
        self.InpTargetY = Tkinter.StringVar()
        self.TarXentrybox = Tkinter.Entry(self.TargetXYEntryFrame,textvariable=self.InpTargetX)
        self.TarYentrybox = Tkinter.Entry(self.TargetXYEntryFrame,textvariable=self.InpTargetY)
        self.TarXentrybox.grid(column=0,row=0,sticky='EW')
        self.TarYentrybox.grid(column=1,row=0,sticky='EW')
        self.TarXentrybox.selection_range(0, Tkinter.END)        
        self.TarYentrybox.selection_range(0, Tkinter.END)        


        #Current AutoGuider's X Y position entry boxes
        self.CurAGXYEntryFrame = ttk.Labelframe(self, text='AutoGuider X Y')
        self.CurAGXYEntryFrame.grid(column=0,row=3,columnspan=2)
        self.InpCurAGX = Tkinter.StringVar()
        self.InpCurAGY = Tkinter.StringVar()
        self.AGXentrybox = Tkinter.Entry(self.CurAGXYEntryFrame,textvariable=self.InpCurAGX)
        self.AGYentrybox = Tkinter.Entry(self.CurAGXYEntryFrame,textvariable=self.InpCurAGY)
        self.AGXentrybox.grid(column=0,row=0,sticky='EW')
        self.AGYentrybox.grid(column=1,row=0,sticky='EW')
        self.AGXentrybox.selection_range(0, Tkinter.END)        
        self.AGYentrybox.selection_range(0, Tkinter.END)        


        # The Calculate Autoguider Shift Button
        self.CalculateButton = Tkinter.Button(self,text=u"Calculate",command=self.OnCalculateButtonClick)
        self.CalculateButton.grid(column=0,row=4,columnspan=2)
 
        # Output Label showing The New Autoguider positions to shift to.
        self.OutputText = Tkinter.StringVar()
        self.OutputBox = Tkinter.Label(self,textvariable=self.OutputText,anchor="w",fg="white",bg="black")
        self.OutputBox.grid(column=0,row=5,columnspan=2,sticky='EW')
        self.OutputText.set(u"Fill the boxes and press Calculate")
        

        self.pack(fill='both', expand='yes')


        
    def OnClearButtonClick(self):
        """ Clears all the Entries in each entry box"""
        self.InpCurrentX.set('')
        self.InpCurrentY.set('')
        self.InpTargetX.set('')
        self.InpTargetY.set('')
        self.InpCurAGX.set('')
        self.InpCurAGY.set('')

    def OnCalculateButtonClick(self):
        """ Calculates the shift to be given in Autoguider """
        try:
            cur_X=float(self.InpCurrentX.get())
            cur_Y=float(self.InpCurrentY.get())
            target_X=float(self.InpTargetX.get())
            target_Y=float(self.InpTargetY.get())
            curAG_X=float(self.InpCurAGX.get())
            curAG_Y=float(self.InpCurAGY.get())
        except(ValueError):
            self.OutputText.set(u"Incorrect coordinate entered in one of the boxes!!")
        else :
            V=(target_X-cur_X,target_Y-cur_Y)   # Required displacement vector in TIRSPEC frame
            V_AG=self.LinearTransform(V)  # Calculated displacement vector in Auto guider frame
            
            targetAG_X=curAG_X + V_AG[0]  #Targeted Auto Guider star's X coordinate 
            targetAG_Y=curAG_Y + V_AG[1]  #Targeted Auto Guider star's Y coordinate 
            self.OutputText.set("Move AG star to X={0:.2f}, Y={1:.2f}".format(targetAG_X,targetAG_Y))


    def LinearTransform(self,V,T=(-0.64472891,-0.69205482,0.66728573,-0.65989201)):
        """ Linear Transform the input 1d vector V. Input V is a python 2 element tuple, V=(X,Y). 
    Input T is the Linear Transformation Matrix given in the form of a tuple (a11,a12,a21,a22).  Function will return the Linear transformed  1d vector as another tuple"""
        x,y=V  
        xout=x *T[0] + y *T[1]  #Linear Transforming
        yout=x *T[2] + y *T[3]
        return xout,yout


class OStab(ttk.Frame):
    """ This Tab is the interface for calculating shifts to be given in Observatory Server."""
    def __init__(self,parent):
        ttk.Frame.__init__(self,parent)
        self.parent=parent
        self.title="Observatory Server Control"
        self.grid()
        
        # The Top Clear Button
        self.ClearButton = Tkinter.Button(self,text=u"Clear All",command=self.OnClearButtonClick)
        self.ClearButton.grid(column=0,row=0,columnspan=2)

        # Present Telescope Offset entry boxes
        self.OffsetXYEntryFrame = ttk.Labelframe(self, text='Current Offset Ra Dec')
        self.OffsetXYEntryFrame.grid(column=0,row=1,columnspan=2)
        self.InpOffsetX = Tkinter.StringVar()
        self.InpOffsetY = Tkinter.StringVar()
        self.OffXentrybox = Tkinter.Entry(self.OffsetXYEntryFrame,textvariable=self.InpOffsetX)
        self.OffYentrybox = Tkinter.Entry(self.OffsetXYEntryFrame,textvariable=self.InpOffsetY)
        self.OffXentrybox.grid(column=0,row=0,sticky='EW')
        self.OffYentrybox.grid(column=1,row=0,sticky='EW')
        self.OffXentrybox.selection_range(0, Tkinter.END)        
        self.OffYentrybox.selection_range(0, Tkinter.END)        
        self.InpOffsetX.set('0')
        self.InpOffsetY.set('0')

        # Present Telescope Rotation value
        self.RotEntryFrame = ttk.Labelframe(self, text='Telescope Rotation')
        self.RotEntryFrame.grid(column=0,row=2,columnspan=1)
        self.InpRot = Tkinter.StringVar()
        self.Rotentrybox = Tkinter.Entry(self.RotEntryFrame,textvariable=self.InpRot)
        self.Rotentrybox.grid(column=0,row=0)
        self.Rotentrybox.selection_range(0, Tkinter.END)        
        self.InpRot.set('0')
        

        #Current TIRSPEC's X Y position entry boxes
        self.CurrentXYEntryFrame = ttk.Labelframe(self, text='Current X Y')
        self.CurrentXYEntryFrame.grid(column=0,row=3,columnspan=2)
        self.InpCurrentX = Tkinter.StringVar()
        self.InpCurrentY = Tkinter.StringVar()
        self.CurrXentrybox = Tkinter.Entry(self.CurrentXYEntryFrame,textvariable=self.InpCurrentX)
        self.CurrYentrybox = Tkinter.Entry(self.CurrentXYEntryFrame,textvariable=self.InpCurrentY)
        self.CurrXentrybox.grid(column=0,row=0,sticky='EW')
        self.CurrYentrybox.grid(column=1,row=0,sticky='EW')
        self.CurrXentrybox.selection_range(0, Tkinter.END)        
        self.CurrYentrybox.selection_range(0, Tkinter.END)        


        #Target TIRSPEC's X Y position entry boxes
        self.TargetXYEntryFrame = ttk.Labelframe(self, text='Target X Y')
        self.TargetXYEntryFrame.grid(column=0,row=4,columnspan=2)
        self.InpTargetX = Tkinter.StringVar()
        self.InpTargetY = Tkinter.StringVar()
        self.TarXentrybox = Tkinter.Entry(self.TargetXYEntryFrame,textvariable=self.InpTargetX)
        self.TarYentrybox = Tkinter.Entry(self.TargetXYEntryFrame,textvariable=self.InpTargetY)
        self.SlitPositionButton = Tkinter.Button(self.TargetXYEntryFrame,text=u"Get Slit Pos",command=self.OnGetSlitPositionButtonClick)
        self.TarXentrybox.grid(column=0,row=0,sticky='EW')
        self.TarYentrybox.grid(column=1,row=0,sticky='EW')
        self.SlitPositionButton.grid(column=2,row=0)
        self.TarXentrybox.selection_range(0, Tkinter.END)        
        self.TarYentrybox.selection_range(0, Tkinter.END)        


        # The Calculate Telescope Shift command Button
        self.CalculateButton = Tkinter.Button(self,text=u"New Shift Values",command=self.OnCalculateButtonClick)
        self.CalculateButton.grid(column=0,row=5,columnspan=2)
 
        # Output Label showing The New Autoguider positions to shift to.
        self.OutputText = Tkinter.StringVar()
        self.OutputBox = Tkinter.Label(self,textvariable=self.OutputText,anchor="w",fg="white",bg="black")
        self.OutputBox.grid(column=0,row=6,columnspan=2,sticky='EW')
        self.OutputText.set(u"Fill the boxes and press New Shift Values")
        

        self.pack(fill='both', expand='yes')


        
    def OnClearButtonClick(self):
        """ Clears all the Entries in each entry box"""
        self.InpOffsetX.set('0')
        self.InpOffsetY.set('0')
        self.InpRot.set('0')
        self.InpCurrentX.set('')
        self.InpCurrentY.set('')
        self.InpTargetX.set('')
        self.InpTargetY.set('')

    def OnGetSlitPositionButtonClick(self):
        """ Updates the Target Box entry with slit position coordinate """
        try:
            target_X,target_Y=self.GetSlitPosition(self.InpTargetY.get(),X=self.InpTargetX.get())
        except(KeyError):
            self.OutputText.set(u"Unknown slit name entered in Target Y box.")
        else :
            self.InpTargetX.set('{0:.2f}'.format(target_X))
            self.InpTargetY.set('{0:.2f}'.format(target_Y))
        
    def OnCalculateButtonClick(self):
        """ Calculates the shift to be given in telescope offset command of Observatory server """
        try:
            cur_X=float(self.InpCurrentX.get())
            cur_Y=float(self.InpCurrentY.get())
            target_X=float(self.InpTargetX.get())
            target_Y=float(self.InpTargetY.get())
            OffsetRA=float(self.InpOffsetX.get())
            OffsetDec=float(self.InpOffsetY.get())
            Rot=float(self.InpRot.get())
        except(ValueError):
            self.OutputText.set(u"Incorrect coordinate entered in one of the boxes!!")
        else :
            currpix=(cur_X,cur_Y)
            destpix=(target_X,target_Y)
            OffsetRA,OffsetDec=self.CalculateNewOffset(currpix,destpix,OffsetRA,OffsetDec,rot=Rot)
            self.OutputText.set(u"RA Offset: {0:.2f}  Dec Offset: {1:.2f}".format(OffsetRA,OffsetDec))
            self.InpOffsetX.set('{0:.2f}'.format(OffsetRA))
            self.InpOffsetY.set('{0:.2f}'.format(OffsetDec))


    def is_number(self,s): 
        """A funtion to check wheter string s is a number or not."""
        try:
            float(s)
            return True
        except(ValueError,TypeError):
            return False
        

    def GetSlitPosition(self,slit,X=None):
        """ Returns the X,Y coordinate to position star inside the slit . If the slit position changes update them inside this function."""
        slit=slit.upper()
                
        # Using the X coordinate if user has already specified it else default value 510
        PosX= float(X) if self.is_number(X) else 510.0
        
        #Slit position coefficents of a 3 degree polynomial x^2 *C[0] + x*C[1] + C[2]
        S1pC=(0, 1.45341055e-02,   5.50352476e+02)
        S2pC=(0, 1.71960827e-02,   5.48740933e+02) 
        S3pC=(0, 4.02612613e-02,   5.33352515e+02)
        S4pC=(0, 2.40138996e-02,   5.44706099e+02) 
        S5pC=(0, -1.77628424e-01,   6.45344322e+02) 

        L1pC=(4.81281533e-06,   1.20302096e-02,   5.50744714e+02)
        L2pC=(-4.06345776e-06,   1.52502824e-02,   5.53112549e+02)
        L3pC=(4.59518777e-06,   1.08348073e-02,   5.51389802e+02) 
        L4pC=(-5.05144414e-06,   2.01066378e-02,   5.51651474e+02)
        L5pC=(4.27444671e-06,   9.47632824e-03,   5.49040853e+02) 

        SlitDict={'S1':S1pC, 'S2':S2pC, 'S3':S3pC, 'S4':S4pC, 'S5': S5pC,'L1':L1pC, 'L2':L2pC, 'L3':L3pC, 'L4':L4pC, 'L5': L5pC}
        
        PosY=SlitDict[slit][0]*PosX**2 +SlitDict[slit][1]*PosX +SlitDict[slit][2]

        return PosX,PosY
        
    def RotateScaleVector(self,V,S=(1,1),rot=0):
        """ Rotated the input 1d vector V. Input V is a python 2 element tuple, V=(X,Y). 
        Input S is the scaling vector to be applied to input 1d vector before rotation.
        Input Rot is the rotation angle in radians  
        Function will return the rotated  1d vector as another tuple"""
        x=S[0]*V[0]   #Scaling
        y=S[1]*V[1]
        xout=x *math.cos(rot) - y *math.sin(rot)  #Rotating
        yout=x *math.sin(rot) + y *math.cos(rot)
        return xout,yout


    def CalculateNewOffset(self,Start,Destination,OffsetRA,OffsetDec,rot=0):
        """ Start is a tuple (X,Y) and Destination is the final position tuple (X,Y)
        OffsetRA and OffsetDec is previous offset in RA and Dec.
        Rot is the rotation given to telescope by  telescopesetrotator command in degrees.
        This function returns the new Offsets in RA and Dec"""
        PixelScale=0.3081
        RApos,DECpos=self.RotateScaleVector((Start[0]-Destination[0],Start[1]-Destination[1]),S=(PixelScale,PixelScale),rot=math.radians(rot))
        RApos+=OffsetRA
        DECpos+=OffsetDec
        return RApos,DECpos
        


class MainWindow(Tkinter.Tk):
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()
    
    def initialize(self):
        FullNotebook = ttk.Notebook(self)
        Tab4AG = AGtab(FullNotebook)
        Tab4OS = OStab(FullNotebook)
        FullNotebook.add(Tab4AG, text=Tab4AG.title)
        FullNotebook.add(Tab4OS, text=Tab4OS.title)
        FullNotebook.pack(fill='both', expand='yes')
        self.grid()

        self.update()
        self.geometry(self.geometry())   

if __name__ == "__main__":
    GUI = MainWindow(None)
    GUI.title('TIRSPEC Telescope Shift Assistant')
    GUI.mainloop()
