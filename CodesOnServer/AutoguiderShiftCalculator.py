#!/usr/bin/env python
#This script is to calculate the shift user has to give
# in order to shift the guiding source in the TIRSPEC frame.
import math     # NOT USED while running
import sys
def RotateScaleVector(V,S=(1,1),rot=0):   # NOT USED while running
    """ Rotated the input 1d vector V. Input V is a python 2 element tuple, V=(X,Y). 
    Input S is the scaling vector to be applied to input 1d vector before rotation.
    Input Rot is the rotation angle in radians  
    Function will return the rotated  1d vector as another tuple"""
    x=S[0]*V[0]   #Scaling
    y=S[1]*V[1]
    xout=x *math.cos(rot) - y *math.sin(rot)  #Rotating
    yout=x *math.sin(rot) + y *math.cos(rot)
    return xout,yout

def LinearTransform(V,T=(-0.64472891,-0.69205482,0.66728573,-0.65989201)):
    """ Linear Transform the input 1d vector V. Input V is a python 2 element tuple, V=(X,Y). 
    Input T is the Linear Transformation Matrix given in the form of a tuple (a11,a12,a21,a22).
    Function will return the Linear transformed  1d vector as another tuple"""
    x,y=V  
    xout=x *T[0] + y *T[1]  #Linear Transforming
    yout=x *T[2] + y *T[3]
    return xout,yout


while True :
    try:
        CurrentTxy=raw_input("Present X Y coordinate in TIRSPEC (Eg: 450 556) : ").strip(' ') #strips away leading and trailing whitespaces.
        TargetTxy=raw_input("Targeted X Y coordinate in TIRSPEC (Eg: 490 556) : ").strip(' ')
        CurrentAGxy=raw_input("Present X Y coordinate in AUTOGUIDER (Eg: 1799.25 234.35) : ").strip(' ')
    except (KeyboardInterrupt, EOFError):
        print("\n Exiting.. \n ")
        sys.exit(0)

    try :
        cur_X, cur_Y = tuple(float(i) for i in CurrentTxy.split())   #Current TIRSPEC star's X,Y coordinate
        target_X, target_Y = tuple(float(i) for i in TargetTxy.split())    #Targeted TIRSPEC star X,Y coordinate

        curAG_X, curAG_Y = tuple(float(i) for i in CurrentAGxy.split())   #Current Auto Guider star's X,Y coordinate 
    except (ValueError):
        print("Error: Wrong input given. \n---------------------------------------------")
        continue

    V=(target_X-cur_X,target_Y-cur_Y)   # Required displacement vector in TIRSPEC frame
    V_AG=LinearTransform(V)  # Calculated displacement vector in Auto guider frame

    targetAG_X=curAG_X + V_AG[0]  #Targeted Auto Guider star's X coordinate 
    targetAG_Y=curAG_Y + V_AG[1]  #Targeted Auto Guider star's Y coordinate 

    print("***** Move AutoGuider star to X={0}, Y={1}".format(targetAG_X,targetAG_Y))
