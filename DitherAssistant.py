#!/usr/bin/env python
#This script gives the telescope offset command to be given to Observatory server for a particular pattern of dithering.
#----------------------------------------------indiajoe
import readline
import os
#Initializing
XshiftOLD=0.0
YshiftOLD=0.0

Xshift=0.0
Yshift=0.0

#Some viewer which will update the images automatically once overwritten. Eg: gthumb
Imageviewer='display'
#Imageviewer='gthumb'

#Size of image set to 600x600 pixel with dot at center. I arc sec taken to be 5 pixels
Clearconvertargs=" -size 600x600 xc:black  -fill white -draw 'circle 300,300 300,303' "
convertargs=Clearconvertargs
Outputimage=" ~/Desktop/DitherPattern.png"
os.system('convert '+convertargs+' '+Outputimage)
os.system(Imageviewer+' '+Outputimage+' &')

defaultT='c,r15,u15,l15,d15'


def Moveinstruction(Xfinal,Yfinal,convertargs):
    """ Update the image with new diagram aswell as print the telescopeOffset command"""
    os.system('convert '+convertargs+' '+Outputimage)    
    if (Xfinal==0 and Yfinal == 0): print("======Clearing offsets to Zero=========")
    print('\033[94m>>\033[92m telescopeOffset 0 0 \033[91m %f %f \033[0m 0 0 0 0'%(-1*Xfinal,Yfinal))

def Appendconvertargs(XshiftOLD,YshiftOLD,Xshift,Yshift,convertargs):
    """ Appends the required arguments for convert command """
    
#    Pixelscale=5
    Center=(300,300)
    Xbeg=str(Center[0]+XshiftOLD)
    Ybeg=str(Center[1]-YshiftOLD)
    Xend=str(Center[0]+Xshift)
    Yend=str(Center[1]-Yshift)
    PointCircle=Xend+','+Yend+' '+Xend+','+str(Center[1]-Yshift+2)
    toappend=" -stroke white -draw 'line "+Xbeg+","+Ybeg+" "+Xend+","+Yend+"' -fill red -draw 'circle "+PointCircle+"'"
    return convertargs+toappend

while 1 :
    try:
        InpPattern=raw_input('Enter dither pattern (default: %s) :'%defaultT).strip(' ')  #strips away leading and trailing whitespaces.
    except (KeyboardInterrupt, EOFError):
        print("\n Exiting.. \n ")
        exit(0)
    if not InpPattern : InpPattern=defaultT
    else: defaultT=InpPattern

    Positions2Move=InpPattern.split(',')
    for pos in Positions2Move :
        targetco=pos.split(' ')
        for coord in targetco:
            try:
                if coord[0] == 'c' : #Clear everything to center
                    Xshift=0
                    Yshift=0
                    XshiftOLD=0
                    YshiftOLD=0
                    convertargs=Clearconvertargs
                elif coord[0] == 'r' : #Moving towards right
                    xinc=eval(coord[1:])
                    Xshift+=xinc
                elif coord[0] == 'l' : #Moving towards left
                    xinc=eval(coord[1:])
                    Xshift-=xinc
                elif coord[0] == 'u' : #Moving towards up
                    yinc=eval(coord[1:])
                    Yshift+=yinc
                elif coord[0] == 'd' : #Moving towards down
                    yinc=eval(coord[1:])
                    Yshift-=yinc
            except(SyntaxError,NameError):
                print("\n You typed input in wrong syntax. \n Please correct it and reenter. \n ")
                continue
                
        #Appending the arguments for convert command for drawing diagram.
        convertargs=Appendconvertargs(XshiftOLD,YshiftOLD,Xshift,Yshift,convertargs)
        XshiftOLD=Xshift
        YshiftOLD=Yshift
        Moveinstruction(Xshift,Yshift,convertargs)
        
        
                


                

        
        
