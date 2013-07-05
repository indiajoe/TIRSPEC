#!/usr/bin/env python
#This script takes current position of source and gives the input to be given to keystone telscope control for shifting star to slit.
import readline

S1p=(515,558.4)
S2p=(515,558.5) #To Verify
S3p=(517,554.4)
S4p=(511,557.5) #To Verify
S5p=(513,556.0) #To Verify
L1p=(520,558.6)
L2p=(520,560.0)
L3p=(520,558.5) #To Verify
L4p=(520,559.1)
L5p=(520,556.0) #To Verify

defaultT='L1'

SlitDict={'S1':S1p, 'S2':S2p, 'S3':S3p, 'S4':S4p, 'S5': S5p,'L1':L1p, 'L2':L2p, 'L3':L3p, 'L4':L4p, 'L5': L5p}

def Moveinstruction(Start,Destination):
    """ Start is a tuple (X,Y) and Destination is the final position tuple (X,Y)"""
    PixelScale=0.3081
    print('\033[94m>>\033[92m Along Axis:\033[91m %f\033[0m arcsec \033[92m Across Axis:\033[91m %f\033[0m arcsec'%((Destination[0]-Start[0])*PixelScale,(Destination[1]-Start[1])*PixelScale))

while 1 :
    try:
        Target=raw_input('Targeted Location (default: %s) :'%defaultT).strip(' ')  #strips away leading and trailing whitespaces.
        Current=raw_input('Current Location (Eg: X Y) :').strip(' ')
    except (KeyboardInterrupt, EOFError):
        print("\n Exiting.. \n ")
        exit(0)
    if not Target : Target=defaultT
    else: defaultT=Target

    if len(Current.split()) != 2 : #Sanity check
        print('INPUT Error: Current pixel location should be in format "X Y"')
        continue
    currpix=tuple([float(x) for x in Current.split()])
    if Target in SlitDict.keys():
        print('To move to the slit %s '%Target +':'+str(SlitDict[Target]))
        Moveinstruction(currpix,SlitDict[Target])
    elif len(Target.split()) == 2 :
        destpix=tuple([float(x) for x in Target.split()])
        print('To move to the pixel '+Target)
        Moveinstruction(currpix,destpix)
    else :
        print('INPUT Error: Target pixel location should be in format "X Y" or it should be name of a slit like S1,S2.. or L1,L2 etc..')
        continue
        
        