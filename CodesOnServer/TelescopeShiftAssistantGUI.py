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
        ClearButtonToolTip = ToolTip(self.ClearButton, follow_mouse=1, text="Resets all entries")
        self.ClearButton.grid(column=0,row=0,columnspan=2)

        #Current TIRSPEC's X Y position entry boxes
        self.CurrentXYEntryFrame = ttk.Labelframe(self, text='Current X Y')
        self.CurrentXYEntryFrame.grid(column=0,row=1,columnspan=2)
        self.InpCurrentX = Tkinter.StringVar()
        self.InpCurrentY = Tkinter.StringVar()
        self.CurrXentrybox = Tkinter.Entry(self.CurrentXYEntryFrame,textvariable=self.InpCurrentX)
        CurrXentryboxToolTip = ToolTip(self.CurrXentrybox, follow_mouse=1, text="Enter Current X pixel coordinate of star.")
        self.CurrYentrybox = Tkinter.Entry(self.CurrentXYEntryFrame,textvariable=self.InpCurrentY)
        CurrYentryboxToolTip = ToolTip(self.CurrYentrybox, follow_mouse=1, text="Enter Current Y pixel coordinate of star.")
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
        TarXentryboxToolTip = ToolTip(self.TarXentrybox, follow_mouse=1, text="Enter Targeted X pixel coordinate.")
        self.TarYentrybox = Tkinter.Entry(self.TargetXYEntryFrame,textvariable=self.InpTargetY)
        TarYentryboxToolTip = ToolTip(self.TarYentrybox, follow_mouse=1, text="Enter Targeted Y pixel coordinate.")
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
        AGXentryboxToolTip = ToolTip(self.AGXentrybox, follow_mouse=1, text="Enter Current AutoGuider's X coordinate.")
        self.AGYentrybox = Tkinter.Entry(self.CurAGXYEntryFrame,textvariable=self.InpCurAGY)
        AGYentryboxToolTip = ToolTip(self.AGYentrybox, follow_mouse=1, text="Enter Current AutoGuider's Y coordinate.")
        self.AGXentrybox.grid(column=0,row=0,sticky='EW')
        self.AGYentrybox.grid(column=1,row=0,sticky='EW')
        self.AGXentrybox.selection_range(0, Tkinter.END)        
        self.AGYentrybox.selection_range(0, Tkinter.END)        


        # The Calculate Autoguider Shift Button
        self.CalculateButton = Tkinter.Button(self,text=u"Calculate",command=self.OnCalculateButtonClick)
        CalculateButtonToolTip = ToolTip(self.CalculateButton, follow_mouse=1, text="Calculate new AutoGuider X Y coordinates you should shift to.")
        self.CalculateButton.grid(column=0,row=4,columnspan=2)
 
        # Output Label showing The New Autoguider positions to shift to.
        self.OutputText = Tkinter.StringVar()
        self.OutputBox = Tkinter.Label(self,textvariable=self.OutputText,anchor="w",fg="white",bg="black")
        self.OutputBox.grid(column=0,row=5,columnspan=2,sticky='EW')
        self.OutputText.set(u"Fill the boxes and press Calculate")
        



        
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
        ClearButtonToolTip = ToolTip(self.ClearButton, follow_mouse=1, text="Resets all entries")
        self.ClearButton.grid(column=0,row=0,columnspan=2)

        # Present Telescope Offset entry boxes
        self.OffsetXYEntryFrame = ttk.Labelframe(self, text='Current Offset Ra Dec')
        self.OffsetXYEntryFrame.grid(column=0,row=1,columnspan=2)
        self.InpOffsetX = Tkinter.StringVar()
        self.InpOffsetY = Tkinter.StringVar()
        self.OffXentrybox = Tkinter.Entry(self.OffsetXYEntryFrame,textvariable=self.InpOffsetX)
        OffXentryboxToolTip = ToolTip(self.OffXentrybox, follow_mouse=1, text="Enter previous RA offset given to telescope.")
        self.OffYentrybox = Tkinter.Entry(self.OffsetXYEntryFrame,textvariable=self.InpOffsetY)
        OffYentryboxToolTip = ToolTip(self.OffYentrybox, follow_mouse=1, text="Enter previous Dec offset given to telescope.")
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
        RotentryboxToolTip = ToolTip(self.Rotentrybox, follow_mouse=1, text="Enter Rotation given to telescope.")
        self.Rotentrybox.grid(column=0,row=0)
        self.Rotentrybox.selection_range(0, Tkinter.END)        
        self.InpRot.set('0')
        

        #Current TIRSPEC's X Y position entry boxes
        self.CurrentXYEntryFrame = ttk.Labelframe(self, text='Current X Y')
        self.CurrentXYEntryFrame.grid(column=0,row=3,columnspan=2)
        self.InpCurrentX = Tkinter.StringVar()
        self.InpCurrentY = Tkinter.StringVar()
        self.CurrXentrybox = Tkinter.Entry(self.CurrentXYEntryFrame,textvariable=self.InpCurrentX)
        CurrXentryboxToolTip = ToolTip(self.CurrXentrybox, follow_mouse=1, text="Enter Current X pixel coordinate of star.")
        self.CurrYentrybox = Tkinter.Entry(self.CurrentXYEntryFrame,textvariable=self.InpCurrentY)
        CurrYentryboxToolTip = ToolTip(self.CurrYentrybox, follow_mouse=1, text="Enter Current Y pixel coordinate of star.")
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
        TarXentryboxToolTip = ToolTip(self.TarXentrybox, follow_mouse=1, text="Enter Targeted X pixel coordinate.")
        self.TarYentrybox = Tkinter.Entry(self.TargetXYEntryFrame,textvariable=self.InpTargetY)
        TarYentryboxToolTip = ToolTip(self.TarYentrybox, follow_mouse=1, text="Enter Targeted Y pixel coordinate. \nOr name of slit [+/- optional offset] to get slit position.\n Eg:    L2 +3.1")
        self.SlitPositionButton = Tkinter.Button(self.TargetXYEntryFrame,text=u"Get Slit Pos",command=self.OnGetSlitPositionButtonClick)
        SlitPositionButtonToolTip = ToolTip(self.SlitPositionButton, follow_mouse=1, text="Update the boxes with slit position.")
        self.TarXentrybox.grid(column=0,row=0,sticky='EW')
        self.TarYentrybox.grid(column=1,row=0,sticky='EW')
        self.SlitPositionButton.grid(column=2,row=0)
        self.TarXentrybox.selection_range(0, Tkinter.END)        
        self.TarYentrybox.selection_range(0, Tkinter.END)        


        # The Calculate Telescope Shift command Button
        self.CalculateButton = Tkinter.Button(self,text=u"New Shift Values",command=self.OnCalculateButtonClick)
        CalculateButtonToolTip = ToolTip(self.CalculateButton, follow_mouse=1, text="Calculate new telescope offsets, and update the boxes.")
        self.CalculateButton.grid(column=0,row=5,columnspan=2)
 
        # Output Label showing The New Autoguider positions to shift to.
        self.OutputText = Tkinter.StringVar()
        self.OutputBox = Tkinter.Label(self,textvariable=self.OutputText,anchor="w",fg="white",bg="black")
        self.OutputBox.grid(column=0,row=6,columnspan=2,sticky='EW')
        self.OutputText.set(u"Fill the boxes and press New Shift Values")
        



        
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
            self.InpCurrentX.set(self.InpTargetX.get())  # Automatically updating Current position to Target position.
            self.InpCurrentY.set(self.InpTargetY.get())  # This will prevent values from changing in accidental double click.


    def is_number(self,s): 
        """A funtion to check wheter string s is a number or not."""
        try:
            float(s)
            return True
        except(ValueError,TypeError):
            return False
        

    def GetSlitPosition(self,slit,X=None):
        """ Returns the X,Y coordinate to position star inside the slit . If the slit position changes update them inside this function. You can give Y axis offset also to the slit name. Eg: L1+3.1 """
        slit=slit.upper()
        offset=0
        #If user has given any offest along with slit name
        if ('+' in slit) and self.is_number(slit.split('+')[-1]) : 
            offset=float(slit.split('+')[-1])
            slit=slit.split('+')[0].strip()
        elif ('-' in slit) and self.is_number(slit.split('-')[-1]) :  
            offset=-1*float(slit.split('-')[-1])
            slit=slit.split('-')[0].strip()            
                
        # Using the X coordinate if user has already specified it else default value 510
        PosX= float(X) if self.is_number(X) else 510.0
        
        #Slit position coefficents of a 3 degree polynomial x^2 *C[0] + x*C[1] + C[2]
        S1pC=(7.78861088e-06,  -2.54277935e-02,   5.71895049e+02)#(0, 1.45341055e-02,   5.50352476e+02)
        S2pC=(-2.20255755e-06,   3.74638507e-02,   5.35639543e+02)#(0, 1.71960827e-02,   5.48740933e+02) 
        S3pC=(-6.72455275e-06,   2.72933649e-02,   5.47496337e+02)#(0, 4.02612613e-02,   5.33352515e+02)
        S4pC=(-4.35836377e-06,   2.87768257e-02,   5.46352457e+02)#(0, 2.40138996e-02,   5.44706099e+02) 
        S5pC=(0,   0,   5.59e+02)#(0, -1.77628424e-01,   6.45344322e+02) 

        L1pC=(4.53580932e-06,   1.06911199e-02,   5.53509356e+02)#(4.81281533e-06,   1.20302096e-02,   5.50744714e+02)
        L2pC=(-4.05115016e-06,   1.50115644e-02,   5.55404684e+02)#(-4.06345776e-06,   1.52502824e-02,   5.53112549e+02)
        L3pC=(4.82560669e-06,   1.04908938e-02,   5.54319255e+02)#(4.59518777e-06,   1.08348073e-02,   5.51389802e+02) 
        L4pC=(-4.92560367e-06,   1.84712028e-02,   5.55983179e+02)#(-5.05144414e-06,   2.01066378e-02,   5.51651474e+02)
        L5pC=(4.27444671e-06,   9.47632824e-03,   5.49040853e+02) 

        SlitDict={'S1':S1pC, 'S2':S2pC, 'S3':S3pC, 'S4':S4pC, 'S5': S5pC,'L1':L1pC, 'L2':L2pC, 'L3':L3pC, 'L4':L4pC, 'L5': L5pC}
        
        PosY=SlitDict[slit][0]*PosX**2 +SlitDict[slit][1]*PosX +SlitDict[slit][2] +offset

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
        RApos,DECpos=self.RotateScaleVector((Start[0]-Destination[0],Start[1]-Destination[1]),S=(PixelScale,-1*PixelScale),rot=math.radians(rot))
        DECpos*=-1     # Flipping the axis back..
        RApos+=OffsetRA
        DECpos+=OffsetDec
        return RApos,DECpos


# Below is a Contributed widget in TKinter  to show Tips in each box.
# http://tkinter.unpythonic.net/wiki/ContributedWidgets?highlight=%28%28ToolTip%29%29
#--------------------------------------------------------------------------ToolTip Addon-----Starts HERE.....
'''Michael Lange <klappnase (at) freakmail (dot) de>
The ToolTip class provides a flexible tooltip widget for Tkinter; it is based on IDLE's ToolTip
module which unfortunately seems to be broken (at least the version I saw).
INITIALIZATION OPTIONS:
anchor :        where the text should be positioned inside the widget, must be on of "n", "s", "e", "w", "nw" and so on;
                default is "center"
bd :            borderwidth of the widget; default is 1 (NOTE: don't use "borderwidth" here)
bg :            background color to use for the widget; default is "lightyellow" (NOTE: don't use "background")
delay :         time in ms that it takes for the widget to appear on the screen when the mouse pointer has
                entered the parent widget; default is 1500
fg :            foreground (i.e. text) color to use; default is "black" (NOTE: don't use "foreground")
follow_mouse :  if set to 1 the tooltip will follow the mouse pointer instead of being displayed
                outside of the parent widget; this may be useful if you want to use tooltips for
                large widgets like listboxes or canvases; default is 0
font :          font to use for the widget; default is system specific
justify :       how multiple lines of text will be aligned, must be "left", "right" or "center"; default is "left"
padx :          extra space added to the left and right within the widget; default is 4
pady :          extra space above and below the text; default is 2
relief :        one of "flat", "ridge", "groove", "raised", "sunken" or "solid"; default is "solid"
state :         must be "normal" or "disabled"; if set to "disabled" the tooltip will not appear; default is "normal"
text :          the text that is displayed inside the widget
textvariable :  if set to an instance of Tkinter.StringVar() the variable's value will be used as text for the widget
width :         width of the widget; the default is 0, which means that "wraplength" will be used to limit the widgets width
wraplength :    limits the number of characters in each line; default is 200

WIDGET METHODS:
configure(**opts) : change one or more of the widget's options as described above; the changes will take effect the
                    next time the tooltip shows up; NOTE: follow_mouse cannot be changed after widget initialization

Other widget methods that might be useful if you want to subclass ToolTip:
enter() :           callback when the mouse pointer enters the parent widget
leave() :           called when the mouse pointer leaves the parent widget
motion() :          is called when the mouse pointer moves inside the parent widget if follow_mouse is set to 1 and the
                    tooltip has shown up to continually update the coordinates of the tooltip window
coords() :          calculates the screen coordinates of the tooltip window
create_contents() : creates the contents of the tooltip window (by default a Tkinter.Label)
'''
# Ideas gleaned from PySol

import Tkinter

class ToolTip:
    def __init__(self, master, text='Your text here', delay=1500, **opts):
        self.master = master
        self._opts = {'anchor':'center', 'bd':1, 'bg':'lightyellow', 'delay':delay, 'fg':'black',\
                      'follow_mouse':0, 'font':None, 'justify':'left', 'padx':4, 'pady':2,\
                      'relief':'solid', 'state':'normal', 'text':text, 'textvariable':None,\
                      'width':0, 'wraplength':200}
        self.configure(**opts)
        self._tipwindow = None
        self._id = None
        self._id1 = self.master.bind("<Enter>", self.enter, '+')
        self._id2 = self.master.bind("<Leave>", self.leave, '+')
        self._id3 = self.master.bind("<ButtonPress>", self.leave, '+')
        self._follow_mouse = 0
        if self._opts['follow_mouse']:
            self._id4 = self.master.bind("<Motion>", self.motion, '+')
            self._follow_mouse = 1
    
    def configure(self, **opts):
        for key in opts:
            if self._opts.has_key(key):
                self._opts[key] = opts[key]
            else:
                KeyError = 'KeyError: Unknown option: "%s"' %key
                raise KeyError
    
    ##----these methods handle the callbacks on "<Enter>", "<Leave>" and "<Motion>"---------------##
    ##----events on the parent widget; override them if you want to change the widget's behavior--##
    
    def enter(self, event=None):
        self._schedule()
        
    def leave(self, event=None):
        self._unschedule()
        self._hide()
    
    def motion(self, event=None):
        if self._tipwindow and self._follow_mouse:
            x, y = self.coords()
            self._tipwindow.wm_geometry("+%d+%d" % (x, y))
    
    ##------the methods that do the work:---------------------------------------------------------##
    
    def _schedule(self):
        self._unschedule()
        if self._opts['state'] == 'disabled':
            return
        self._id = self.master.after(self._opts['delay'], self._show)

    def _unschedule(self):
        id = self._id
        self._id = None
        if id:
            self.master.after_cancel(id)

    def _show(self):
        if self._opts['state'] == 'disabled':
            self._unschedule()
            return
        if not self._tipwindow:
            self._tipwindow = tw = Tkinter.Toplevel(self.master)
            # hide the window until we know the geometry
            tw.withdraw()
            tw.wm_overrideredirect(1)

            if tw.tk.call("tk", "windowingsystem") == 'aqua':
                tw.tk.call("::tk::unsupported::MacWindowStyle", "style", tw._w, "help", "none")

            self.create_contents()
            tw.update_idletasks()
            x, y = self.coords()
            tw.wm_geometry("+%d+%d" % (x, y))
            tw.deiconify()
    
    def _hide(self):
        tw = self._tipwindow
        self._tipwindow = None
        if tw:
            tw.destroy()
                
    ##----these methods might be overridden in derived classes:----------------------------------##
    
    def coords(self):
        # The tip window must be completely outside the master widget;
        # otherwise when the mouse enters the tip window we get
        # a leave event and it disappears, and then we get an enter
        # event and it reappears, and so on forever :-(
        # or we take care that the mouse pointer is always outside the tipwindow :-)
        tw = self._tipwindow
        twx, twy = tw.winfo_reqwidth(), tw.winfo_reqheight()
        w, h = tw.winfo_screenwidth(), tw.winfo_screenheight()
        # calculate the y coordinate:
        if self._follow_mouse:
            y = tw.winfo_pointery() + 20
            # make sure the tipwindow is never outside the screen:
            if y + twy > h:
                y = y - twy - 30
        else:
            y = self.master.winfo_rooty() + self.master.winfo_height() + 3
            if y + twy > h:
                y = self.master.winfo_rooty() - twy - 3
        # we can use the same x coord in both cases:
        x = tw.winfo_pointerx() - twx / 2
        if x < 0:
            x = 0
        elif x + twx > w:
            x = w - twx
        return x, y

    def create_contents(self):
        opts = self._opts.copy()
        for opt in ('delay', 'follow_mouse', 'state'):
            del opts[opt]
        label = Tkinter.Label(self._tipwindow, **opts)
        label.grid()
############
#----------------------------------------------------------------ToolTip Addon END here----------------

class MainWindow(Tkinter.Tk):
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()
    
    def initialize(self):
        FullNotebook = ttk.Notebook(self)
        Tab4OS = OStab(FullNotebook)
        Tab4AG = AGtab(FullNotebook)
        FullNotebook.add(Tab4OS, text=Tab4OS.title)
        FullNotebook.add(Tab4AG, text=Tab4AG.title)
        FullNotebook.grid()
        self.grid()

        self.update()
        self.geometry(self.geometry())   

if __name__ == "__main__":
    GUI = MainWindow(None)
    GUI.title('TIRSPEC Telescope Shift Assistant')
    GUI.mainloop()
