#!/usr/bin/env python
# IMP:  Keep ds9 open 
# This script will guide you through the  images in directory.
# Tell the script which combination of images to align , combine etc.
# Usage: ./AlignCombineImages.py Directory1 [directory2 directory3 ..] 'n*.fits'
# Enjoy!!!---------------------- indiajoe@gmail.com

import os
import sys, traceback 
import glob
import readline

if len(sys.argv) < 3 or not os.path.isdir(sys.argv[1]):
    print("-"*30)
    print("Usage : "+sys.argv[0]+" Night_Directories 'wildcards_to_choose_files'")
    print("Eg: "+sys.argv[0]+" 20130621/ 20130622/ 'n*.fits' ")
    print("-"*30)
    exit(1)

#import pyfits
import astropy.io.fits as pyfits
from pyraf import iraf
iraf.imexamine.unlearn()

# If zenity is not working on you machine, use the cat command instead by uncommnting it.
DISPLAY_MENU_CMD = 'zenity --text-info --title "Images in {imgdir}" --filename={imglist} --width 600 --height 250 &' 
#DISPLAY_MENU_CMD = 'echo "Images in {imgdir}" ; cat {imglist}' 

IMGListFile = 'imgs.text'  # This file will be overwritten
pwd=os.getcwd()
wildcard=sys.argv[-1]
dirs=sys.argv[1:-1]
for imgdir in dirs :
    iraf.cd(imgdir)
    while 1 :
        print imgdir
        images=glob.glob(wildcard)
        images.sort()
        imgs_txt=open(IMGListFile,'w')
        i=0
        for img in images:
            hdulist=pyfits.open(img)
            Exptime=hdulist[0].header.get('ITIME')
            Commentx=hdulist[0].header.get('TCOMMENT')
            Obj=hdulist[0].header.get('TARGET')
            Comment=hdulist[0].header.get('DATE')
            Filter=hdulist[0].header.get('UPPER')
            hdulist.close()
            imgs_txt.write(str(i)+' '+img+' '+str(Obj)+' '+str(Exptime)+' '+str(Commentx)+' '+str(Comment)+' '+str(Filter)+'\n')
            i=i+1
    #   Now list is ready, continuing with what to do
        imgs_txt.close()

        junk = os.system(DISPLAY_MENU_CMD.format(**{"imgdir":imgdir,"imglist":IMGListFile}))

        print ('What do you want to do? Enter a = Align ; c = Combine ; s = Weighted average of Exptime ; i = Inspect ; e = Exit and go to Next filter directory')
        choice=raw_input('Enter your Choice (a,c,s,i or e) :')
        if ( choice == "e" ) : break   # breaking to move to next filter image directory
        elif (choice == "i") :
            while 1 :
                imgEX=input('Enter the Sl.No# of Image you want to Examine [Enter a large number not in image list to exit] : ')
                if imgEX >= len(images) : break
                iraf.display(images[imgEX],1)
                iraf.imexam()
        elif (choice == "c" ) :  # Simply combining 
            imgINP=map(int, raw_input('Enter the Sl.No#s of all Images to combine [ Eg: 1,2,3 ] : ').split(','))
            i=len(imgINP)
            if i < 2 : break    # minimum 2 images to combine are not there..
            for j in imgINP : #Sanity check
                if j >= len(images) : 
                    print('Wrong serial number '+str(j)+' given.. exiting the combine task..')
                    break
            inpVAR=','.join([images[j] for j in imgINP])
            name=raw_input('Enter the name for combined image without .fits extension : ')
            iraf.imcombine(input=inpVAR, output=name+'.fits',combine="average")

        elif (choice == "s" ) :  # Summing and finally dividing by total EXPTIME 
            imgINP=map(int, raw_input('Enter the Sl.No#s of all Images to average by weighted sum of EXPTIME [ Eg: 1,2,3 ] : ').split(','))
            i=len(imgINP)
            if i < 2 : break    # minimum 2 images to combine are not there..
            for j in imgINP :
                if j >= len(images) :  #Sanity check
                    print('Wrong serial number '+str(j)+' given.. exiting the combine task..')
                    break
                hdulist=pyfits.open(images[j])
                Exptime=hdulist[0].header.get('EXPTIME')
                Observatory=hdulist[0].header.get('OBSERVAT')
                hdulist.close()
                if Observatory == 'IGO' : iraf.hedit(images[j],'IntSec',Exptime/1000.0,add=1, ver=0)
                else : iraf.hedit(images[j],'IntSec',Exptime,add=1, ver=0)

            inpVAR=','.join([images[j] for j in imgINP])
            name=raw_input('Enter the name for combined image without .fits extension : ')
            iraf.imcombine(input=inpVAR, output=name+'.fits',combine="sum",scale="exp",zero="mode",weight="exp",expname="IntSec")

        elif ( choice == "a" ) :    # Moving forward to do alignment
            inp=input('Enter the Sl.No# of Reference Image :')
            Refimage=images[inp]
            iraf.display(Refimage,1) 
            print ('Press _a_ over some good stars to align, u can press s also, but DONT press r \n')
            imx=iraf.imexam(Stdout=1)
            Xref=eval(imx[2].split()[0])    #Using the first star as ref star for cose shift calculation
            Yref=eval(imx[2].split()[1])
            foo=open(Refimage+'.coo','w')
            i=2
            while i < len(imx) :               
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                i=i+2
            foo.close()
        #   Now asking for images to align with the Reference image
            imgINP=[]
            foo2=open('shifts.in','w')
            imgINP=map(int, raw_input('Enter the Sl.No#s of all Images to shift and align with Ref image[ Eg: 1,2,3 ] : ').split(','))
            i=len(imgINP)
            for j in imgINP :
                if j >= len(images) : 
                    print('Wrong serial number given.. exiting the align task..')
                    break
                iraf.display(images[j],1)
                print ('Press _a_ over the first star you selected in Reference image \n')
                imx=iraf.imexam(Stdout=1)
                Xin=eval(imx[2].split()[0])    #Using the first star as ref star for corse shift calculation
                Yin=eval(imx[2].split()[1])
                foo2.write(str(Xref-Xin)+'   '+str(Yref-Yin)+'\n')
            foo2.close()

            if i < 1 : break    # No images to align
            inVAR=""#Refimage
            outVAR=""# "s"+Refimage
            imgs2align=open('imgs2align.list','w')
            imgs2alignOUT=open('imgs2alignOUT.list','w')
            for j in imgINP:
                imgs2align.write(images[j]+'\n')
                imgs2alignOUT.write('s'+images[j]+'\n')
        #       inVAR=inVAR+','+images[j]
        #       outVAR=outVAR+',s'+images[j]
            imgs2align.close()
            imgs2alignOUT.close()
     #      inVAR=inVAR[1:]    # To remove leading comma
     #      outVAR=outVAR[1:]
            inVAR="@imgs2align.list"
            outVAR="@imgs2alignOUT.list"

            try :
                iraf.imalign(input=inVAR, reference=Refimage, coords=Refimage+".coo", output=outVAR, shifts="shifts.in", interp_type="nearest",boundary_type="constant",trimimages="no")
                ask=raw_input('Do you want to combine these images? [Enter _y_ for yes]: ')
                if ask == "y" : 
                    name=raw_input('Enter name for the combined image without .fits extension: ')
                    iraf.imcombine(input=Refimage+','+outVAR, output=name+'.fits',combine="average")
                else : print ("No combining done..\n")
            except iraf.IrafError, e :
                print ('IRAF ERROR : Some image might be having problem. Remove it and try')
                print e
                print('-'*60)
                traceback.print_exc(file=sys.stdout)
                print('-'*60)



    print ('Moving to next directory\n')
    iraf.cd(pwd)

print ("All aligning and combining of Images Done \n Enjoy !!!----------- indiajoe@gmail.com \n")
