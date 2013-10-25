#!/usr/bin/env python
#This script is to semi-automate basic data reduction of TIRSPEC data.
# IMP:  Keep ds9 open 
#---------------------------------------------indiajoe
import os
import os.path
import glob
#import pyfits  #Depricated and mearged to astropy
import astropy.io.fits as pyfits
import pyfits.convenience
import sys, traceback 
import numpy as np
import warnings
import re
import shlex
import readline
import shutil

# Other modules inserted inside the code are
# pyraf.iraf
# from astropy.io import ascii
# import astropy.table as table  #IMP: Requires astropy version >= 0.3 for vstack function in Task #9 Photmetry


def Photometry():
    iraf.noao(_doprint=0)     #Loading packages noao digiohot apphot daophot
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.ptools(_doprint=0)  #for txdump

    iraf.images(_doprint=0) 
    iraf.immatch(_doprint=0) #Loading for xyxymatch, geomap, geotran of coords
    iraf.imfilter(_doprint=0)  #Loading packages for convolutiong of gauss

    iraf.phot.unlearn()   #Setting everything to default
    iraf.psf.unlearn()
    iraf.allstar.unlearn()
    iraf.daopars.unlearn()
    iraf.datapar.unlearn()
    iraf.fitskypar.unlearn()
    iraf.photpar.unlearn()
    iraf.findpars.unlearn()

    try:
        imgfile=open(MotherDIR+'/Images4Photo.in','r')
    except IOError,e:
        print('Cannot open Images4Photo.in file. Run Task #6 ')
        print(e)
        print('-'*60)
        traceback.print_exc(file=sys.stdout)
        print('-'*60)
        exit

    # Setting flag by checking wheter the size of qphotinput.txt is Zero or not.
    if os.stat(MotherDIR+"/qphotinput.txt")[6]!=0 : QPHOT_todo='Y'
    else : QPHOT_todo='N'


    imgNo=0
    for imgline in imgfile.readlines() :
        imgline=imgline.rstrip()
        imgline=shlex.split(imgline)
        wdir=imgline[0].split('/')[-2]
        if wdir != os.getcwd().split('/')[-1] : #If we are not already in img dir
            iraf.cd(MotherDIR)  #Going back to parent directory
            DIRtogo="/".join(imgline[0].split('/')[:-1]) #Now going to dir of img
            iraf.cd(DIRtogo)
            foo=open(MotherDIR+'/'+OUTPUTfile,'a')    #Appending into tbe log file The begining of a new directory
            foo.write(DIRtogo+' ---------------------------------------- \n')  # '-'*40  To mark begining of a DIrectory
            foo.close()

        
        img=imgline[0].split('/')[-1]
        filterr=imgline[1]
        intime=imgline[2]
        threshold=imgline[3]
        StartUT=pyfits.convenience.getval(img,UTHDR)  #Reading from Headers
        h,m,s=StartUT.split(':')
        StartUT=int(h)*60*60+int(m)*60+int(s)   #Converting to seconds
        print(wdir, img)
        TrueSigma=0.01      #Will be updated later for every img. Here it should be very small quantity
        yxdim=pyfits.getdata(img).shape  #Storing the (Ymax,Xmax) of image
        leftover=glob.glob(img+'.*')
        if len(leftover) > 0 :
            os.system('mkdir -p Leftover')
            for lft in leftover :
                os.system('mv '+lft+' Leftover/ ')


        #Calling Sextracter And find coordinates
        if not os.path.isfile(MotherDIR+"/sextractor.sex") :  #Incase the parameter file is not already created
            Sextractor_subrout(img=imgline[0])

        os.system("sextractor "+img+" -c "+MotherDIR+"/sextractor.sex -PARAMETERS_NAME "+MotherDIR+"/default.param -FILTER_NAME "+MotherDIR+"/default.conv")
#        os.system("awk 'NR>7{print $3,$5,$6}' test.cat | sort -nr | head -30 | awk '{print $2,$3}' > Bright30.coo")
        os.system("awk 'NR>7{if ($7 == 0){print $3,$5,$6}}' test.cat | sort -nr | cut -d' ' -f 2,3 | head -30 > Bright30.coo")
        #Runing xyxymatch and geo map and geotran to create new SourceT.coo , GoodStarsT.coo, BlankSky.coo
        Nmatch=32
        num_lines=0
        while num_lines < 6 :       # Set the number of stars it should atlest mach here....
            os.system("rm "+img+"xymatch.out 2> /dev/null")
            Nmatch=Nmatch-2
            iraf.xyxymatch.unlearn()
            iraf.xyxymatch(input="Bright30.coo",reference=MotherDIR+"/FirstImageTop30.coo",output=img+"xymatch.out", toler=3, matching="triangles",nmatch=Nmatch)
            # Making sure atleast a few stars were matched. otherwise geoxytran will exit with error.
            os.system("gawk '{if ($1 >0){print $0}}' "+img+"xymatch.out > matchedstars.txt") #Removing headers
            num_lines = sum(1 for line in open("matchedstars.txt"))  #Counting number of lines in the xymatch output. it should be more than 15 (size of header)
            print("Number of stars Matched= "+str(num_lines))
            if Nmatch < 5 : 
                print("Failed to find the coordinates for "+img)
                exit
        
        iraf.geomap(input=img+"xymatch.out", database=img+"rtran.db", xmin=1, xmax=yxdim[1], ymin=1, ymax=yxdim[0], interactive=0)
        iraf.geoxytran(input=MotherDIR+"/GoodStars.coo", output=img+"GoodStarsT.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        iraf.geoxytran(input=MotherDIR+"/Source.coo", output=img+"SourceT.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        iraf.geoxytran(input=MotherDIR+"/BlankSky.coo", output=img+"BlankSky.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        if QPHOT_todo=='Y' :
            iraf.geoxytran(input=MotherDIR+"/qphotinput.txt", output=img+"qphotinput.txt",database=img+"rtran.db",transforms=img+"xymatch.out")


        # Sanity check: To remove any new coordinates calculated lying outside image in *.coo 
        coofileLIST=['GoodStarsT.coo','SourceT.coo','BlankSky.coo']
        if QPHOT_todo=='Y' : coofileLIST.append('qphotinput.txt') 
        for coofile in coofileLIST :
            fooIN=open(img+coofile,'r')
            fooOUT=open(img+coofile+'TEMP','w')
            for star in fooIN.readlines():
                if float(star.split()[0]) > 1 and float(star.split()[0]) < yxdim[1] and float(star.split()[1]) > 1 and float(star.split()[1]) < yxdim[0] : fooOUT.write(star)
                else: print(star +": Outside the image field \n")
            fooIN.close()
            fooOUT.close()
            os.system('mv  '+img+coofile+'TEMP  '+img+coofile)

        #---------------------------------
        # Due to small error in calculation of star position, we need to create a more accurate GoodStars.coo and Source.coo
        # Plus we have to remove any saturated stars also.
        imx=iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=img+'GoodStarsT.coo',Stdout=1)
        foo=open(img+'GoodStars.coo','w')    #Creating good stars coords files
        starlist=" "
        i=2
        while i < len(imx) :
            if (imx[i+1].split()[4] != 'INDEF') and (float(imx[i+1].split()[4]) > 10*np.sqrt(float(imx[i+1].split()[3]))) and (float(imx[i+1].split()[4])+float(imx[i+1].split()[3])) < float(DATAMAX) : # (Peak-sky) is not INDEF and (Peak-sky) > 10*sqrt(sky) and  Peak is not saturated
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                starlist=starlist+str(i/2)+"_"   #Saving the string of good stars survived.
            else : print('Discarded: '+str(i/2)+' th number star not good of '+DIRtogo+' '+img)
            i=i+2
        foo.close()

        try :
            imx=iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=img+'SourceT.coo',Stdout=1)
            Xprim=eval(imx[2].split()[0])  
            Yprim=eval(imx[2].split()[1])
            foo=open(img+'Source.coo','w')    #Creating text file contiaing coords of v1647
            i=2
            while i < len(imx) :
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                i=i+2
            foo.close()
        except iraf.IrafError, e :
            print('Iraf Error: going forward with first estimated Source.coo')
            shutil.copy(img+'SourceT.coo',img+'Source.coo')
    
        #---------------------------------END of recalculation of coordinates------


        #Creating all the gauss convolved images if the CONVOLVEIMG variable is _not_ set to NO
        OriginalIMG=img
        convIMGS=[img]
        if CONVOLVEIMG != 'NO' :  # If the CONVOLVEIMG variable is not set to NO
            for si in eval(CONVOLVEIMG) :
                iraf.gauss(input=img,output=img+'_'+str(si)+'.fits',sigma=si)
                convIMGS.append(img+'_'+str(si)+'.fits')       #List of all convolved images on which we have to do photometry
        # Now the loop of  doing the photometry for all the convolved) images
        for img in convIMGS :
            imx=iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=OriginalIMG+'GoodStars.coo',Stdout=1)
            print imx           #DEBUGGING---------------------------------------**
            #Calculating median FWHM
            fwhmlist=[]
            i=3
            while i < len(imx) :               
                fwhmlist.append(eval(imx[i].split()[10]))
                i=i+2
            #Median FWHM is
            fwhm=np.median(fwhmlist)  
            print('Setted value of FWHM =' + str(fwhm))
            #Calculating sky mean and stdev
            imx=iraf.imexam(input=img,frame=1,use_display=0,defkey='m',imagecur=OriginalIMG+'BlankSky.coo',Stdout=1)
            print imx            #DEBUGGING--------------------------------------**
            #Calculating average Sigma and mean of sky
            skylist=[]
            sigmalist=[]
            i=1
            while i < len(imx) :               
                skylist.append(eval(imx[i].split()[2]))
                sigmalist.append(eval(imx[i].split()[4]))
                i=i+1
            #Average mean,sigma and datamin are
            mean=np.median(skylist)
            sigma=np.median(sigmalist)
            datamin= mean - 10*max(sigma,TrueSigma)  #Setting lowerlimit to sigma from going less than TrueSigma of original img
            print('Mean sky = '+str(mean))
            print('Mean sigma = '+str(sigma))
            print('Datamin = '+str(datamin))

            foo=open(MotherDIR+'/'+OUTPUTfile,'a')    #Appending into the log file to write output of photometry
            foo.write(img +'  '+str(round(fwhm,2)) + '  "'+filterr+'"  '+str(intime) +'  '+str(StartUT)+'  '+ str(round(mean,3)) +'  ' + str(round(sigma,3)) +'  '+str(round(datamin,3)) + ' | ')
            foo.close()
            #Now starts photometry processes....
            # First is daopar, Press :w then :q to continue if everything is fine
            psfradi= 4*fwhm +1

            iraf.daopars.setParam('matchrad',fwhm)
            iraf.daopars.setParam('psfrad',psfradi)
            iraf.daopars.setParam('fitrad',fwhm)

            iraf.datapar.setParam('fwhmpsf',fwhm)
            iraf.datapar.setParam('sigma',sigma)
            iraf.datapar.setParam('datamin',datamin)
            iraf.datapar.setParam('datamax',DATAMAX)
            iraf.datapar.setParam('readnoi',READNOISE)
            iraf.datapar.setParam('epadu',EPADU)
            iraf.datapar.setParam('itime',intime)
#            iraf.datapar.setParam('ifilter',filterr)

            iraf.fitskypar.setParam('annulus',eval(ANNULUS))
            iraf.fitskypar.setParam('dannulu',eval(DANNULUS))
            
            aperture = eval(APPERTURE)  # 4*fwhm
            iraf.photpar.setParam('apertur',aperture)
            iraf.findpars.setParam('threshold',threshold)
            if OriginalIMG == img :
                TrueSigma=sigma   #Setting the correct sigma of sky
                iraf.daofind(image=img,output="default",verify=VER)
            else :
                shutil.copy(OriginalIMG+'.coo.1',img+'.coo.1')
            #Going forward to do phot
            iraf.phot(image=img,coords="default",output="default",verify=VER)

            magtable=ascii.read(img+'.mag.1')
            goodstarsFILE=open(OriginalIMG+'GoodStars.coo','r')
            tablelist=[]
            for goodstarXY in goodstarsFILE.readlines():
                gsX,gsY=goodstarXY.rstrip().split()
                gsX=eval(gsX)
                gsY=eval(gsY)
                tablelist.append(magtable[((magtable['XCENTER']-gsX)**2<4) & ((magtable['YCENTER']-gsY)**2<4)]['XCENTER','YCENTER','MAG','ID'])
            goodstarsFILE.close()
            goodstarsTable=table.vstack(tablelist)
            if DOPSF == 'YES' :  # IF PSF photometry has to be done...
                #Creating the imcommands file by finding Star IDs
#                os.system(MotherDIR+'/Finding_StarID_Curser_File.sh ' + img +' '+OriginalIMG+'GoodStars.coo' )  #OLD Way...
                icomFILE=open('icommands.in','w')
                icomFILE.write(':a '+'\n:a '.join([str(sid) for sid in goodstarsTable['ID']])+'\n')
                icomFILE.write('f \nw \nq \n') # Adding lines f w q at the end.
                icomFILE.close()
                print ("Doing psf, Non-interactively.. Using Coords of good star")
                iraf.psf(image=img, pstfile="", photfile="default", psfimage="default", opstfile="default", groupfil="default", icommands='icommands.in', verify=VER)
            #    print ("Doing psf, Interactively.. Use the a a ... f w q  sequence..")
            #    iraf.psf(image=img, pstfile="", photfile="default", psfimage="default", opstfile="default", groupfil="default")

                iraf.allstar(image=img, photfile="default", psfimage="default", allstarf="default", rejfile="default", subimage="default" ,verify=VER )
                print ("Psf photometry over")
                print ("--------------------------------------")

            #Doing the phot again on Source, Just in case Daofind didn't detect it and Good stars...
            iraf.datapar.setParam('datamin',mean-5*max(sigma,TrueSigma))
            iraf.phot(image=img,coords=OriginalIMG+'Source.coo',output="default",verify=VER)
            Sourcemagtable=ascii.read(img+'.mag.2')
#            iraf.phot(image=img,coords=OriginalIMG+'GoodStars.coo',output="default",verify=VER)
#            SecondPhotresults=iraf.txdump(textfiles=img+".mag.2",fields="XCENTER,YCENTER,MAG",expr="yes",Stdout=1)
#            SecondPhotresults.extend(iraf.txdump(textfiles=img+".mag.3",fields="XCENTER,YCENTER,MAG",expr="yes",Stdout=1))

            iraf.hedit(img, "intime", intime, add=1, ver=0)
            #Doing qphot at all the points in qphotinput.txt with the corresponding parameters.
            if QPHOT_todo=='Y' :  #If there exist some qphot sources
                foo=open(OriginalIMG+"qphotinput.txt",'r')
                for qphotobj in foo.readlines():
                    qphotobj=qphotobj.rstrip()
                    obj=qphotobj.split()
                    foo2=open('qphotSource.Tcoo','w')
                    foo2.write(obj[0]+'  '+obj[1])
                    foo2.close()
                    iraf.qphot(image=img , coords='qphotSource.Tcoo', cbox=5, annulus=obj[3], dannulus=obj[4], aperture=obj[2], exposur="intime", epadu=EPADU ,interactive=0 )
                
                foo.close()
            #Now, Writing the Mag to output file
            foo=open(MotherDIR+'/'+OUTPUTfile,'a')
#            os.system(MotherDIR+'/Creating_Log_File.sh '+img+' '+OriginalIMG+'GoodStars.coo'+' '+OriginalIMG+'Source.coo'+' '+MotherDIR+'/'+OUTPUTfile ) 

            #First append the qphot magnitudes to the Photometry output file
            magfiles=glob.glob(img+'.mag.*')
            magfiles.sort()
            for filemag in magfiles:
                if eval(filemag.split('.')[-1]) > 2 : # The qphot output files
                    qphottable=ascii.read(filemag)
                    foo.write(' '+qphottable['MAG'][-1]) #Writing the last Mag in file
            foo.write(' | ')  #Adding a seperator after qphot mags
            # If PSF photometry as done, adding those mags to the file.
            if DOPSF == 'YES' :  # IF PSF photometry was done...
                #First the mags of Source stars
                alstable=ascii.read(img+'.als.1')
                SourcestarsFILE=open(OriginalIMG+'Source.coo','r')
                tablelist=[]
                for sourstarXY in SourcestarsFILE.readlines():
                    ssX,ssY=sourstarXY.rstrip().split()
                    ssX=eval(ssX)
                    ssY=eval(ssY)
                    tablelist.append(alstable[((alstable['XCENTER']-ssX)**2<4) & ((alstable['YCENTER']-ssY)**2<4)]['XCENTER','YCENTER','MAG'])
                SourcestarsFILE.close()
                SourcestarsALSTable=table.vstack(tablelist)
                for rows in SourcestarsALSTable: 
                    foo.write(' %f %f %f'%(rows['XCENTER'],rows['YCENTER'],rows['MAG']))
                foo.write(' | ')  #Adding a seperator after Source Mags
                #Now the psf magnitudes of the good stars
                goodstarsFILE=open(OriginalIMG+'GoodStars.coo','r')
                tablelist=[]
                for goodstarXY in goodstarsFILE.readlines():
                    gsX,gsY=goodstarXY.rstrip().split()
                    gsX=eval(gsX)
                    gsY=eval(gsY)
                    tablelist.append(alstable[((alstable['XCENTER']-gsX)**2<4) & ((alstable['YCENTER']-gsY)**2<4)]['XCENTER','YCENTER','MAG'])
                goodstarsFILE.close()
                goodstarsALSTable=table.vstack(tablelist)
                for rows in goodstarsALSTable:
                    foo.write(' %f %f %f'%(rows['XCENTER'],rows['YCENTER'],rows['MAG']))
                foo.write(' | ')  #Adding a seperator after Good stars X Y Mags
                
            else:
                foo.write(' | | ')
            
            # Writing the pure phot results we calculated into the list before closing the line
            for rows in Sourcemagtable:  #Source Mags
                foo.write(' %f %f %f'%(rows['XCENTER'],rows['YCENTER'],rows['MAG']))
            foo.write(' | ')  #Adding a seperator after Source Mags
            for rows in goodstarsTable:  #Good Stars Mags
                foo.write(' %f %f %f'%(rows['XCENTER'],rows['YCENTER'],rows['MAG']))
            foo.write(' | ')  #Adding a seperator after Good Star Mags

            foo.write(' '+starlist+' \n') # Ending this image line with Good star's ID.
            foo.close()
           
            print ("Photometry of "+img+" over. \n Now proceeding to next image")
            #END of the photometry of convolved images set..
        imgNo=imgNo+1
        foo=open(MotherDIR+'/'+OUTPUTfile,'a')    #Appending into the log file to write output of photometry
        foo.write('-------------------------------------------- \n')  # '-'*44  To mark end of an image
        foo.close()

    #All photometry over
    imgfile.close()
    print("Great...Photometry of all "+str(imgNo)+" images are over...")
    print("Enjoy!!! ---------------------indiajoe@gmail.com ")

def is_number(s):   # A funtion to check wheter string s is a number or not.
    try:
        float(s)
        return True
    except ValueError:
        return False

def Sextractor_subrout(img=None,N=30):
    """ Calls the Sextractor and create the sextractor parameter files if it doesn't already exists. And also create coord file of the brightest N=30 number of stars."""
    N=str(N)
    backupPWD=os.getcwd()
    iraf.cd(MotherDIR)  #Going back to parent directory
    if not os.path.isfile("sextractor.sex") : #If a config file doesn't exist already
        os.system("sextractor -d > sextractor.sex")
        if os.path.isfile("/usr/share/sextractor/default.param") :
            os.system("cp /usr/share/sextractor/default.param /usr/share/sextractor/default.conv .")
        else : print ("Error: cannot find default.param (.conv) in /usr/share/sextractor/ . \n You might have installed sextracter somewhere else")
        print("Sextractor Config file sextractor.sex and default.parm and default.conv created. \n If required u can edit it before calling Photometry")

    if img is None : # If No img is given, then using the first image in Images4Photo.in file
        try:
            imgfile=open(MotherDIR+'/Images4Photo.in','r')
        except IOError,e:
            print('Cannot open Images4Photo.in file. Run Task #6 ')
            print(e)
            print('-'*60)
            traceback.print_exc(file=sys.stdout)
            print('-'*60)
            exit

        imgline=imgfile.readlines()[0]
        imgline=imgline.rstrip()
        img=imgline.split()[0]
        imgfile.close()

    os.system("sextractor "+img+" -c sextractor.sex")
    os.system("awk 'NR>7{if ($7 == 0){print $3,$5,$6}}' test.cat | sort -nr | cut -d' ' -f 2,3 | head -"+N+" > FirstImageTop"+N+".coo")
#    os.system("awk 'NR>7{print $3,$5,$6}' test.cat | sort -nr | head -"+N+" | awk '{print $2,$3}' > FirstImageTop"+N+".coo")
    print("Brightest "+N+" stars coordinates of first image created in FirstImageTop"+N+".coo")
    iraf.cd(backupPWD)


def Star_sky_subrout(img=None) :
    """ Opens the image and create Source.coo, GoodStars.coo, BlankSky.coo, Polygon.coo files"""
    backupPWD=os.getcwd()
    iraf.cd(MotherDIR)  #Going back to parent directory

    if img is None : # If No img is given, then using the first image in Images4Photo.in file
        try:
            imgfile=open(MotherDIR+'/Images4Photo.in','r')
        except IOError,e:
            print('Cannot open Images4Photo.in file. Run Task #6 ')
            print(e)
            print('-'*60)
            traceback.print_exc(file=sys.stdout)
            print('-'*60)
            exit

        imgline=imgfile.readlines()[0]
        imgline=imgline.rstrip()
        img=imgline.split()[0]
        imgfile.close()

    if not ( os.path.isfile("Source.coo") and os.path.isfile("GoodStars.coo") and os.path.isfile("BlankSky.coo") )  : #If the .coo files doesn't exist already
        iraf.display(img,1)
        print ('\n For taking coordinates of Source. Press _a_ over Primary Sources.')
        imx=iraf.imexam(Stdout=1)
        foo=open('Source.coo','w')    #Creating text file contiaing coords of v1647
        i=2
        while i < len(imx) :               
            foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
            i=i+2
        foo.close()
        print ('\n For taking coordinates of good stars. Press _a_ over some good stars. \n Nonsaturated among them will be used for psf fitting.')
        print ('IMP: Press coordinate of Stars in standard required order')
        imx=iraf.imexam(Stdout=1)
        foo=open('GoodStars.coo','w')    #Creating good stars coords files
        i=2
        while i < len(imx) :               
            foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
            i=i+2
        foo.close()
        shutil.copy('GoodStars.coo','GoodStars.cooORIGINAL')   #Keeping BACKUP....
        shutil.copy('Source.coo','Source.cooORIGINAL')   #Keeping BACKUP....
        print ('\n For taking coordinates of good sky. Press _x_ over blank sky areas.')
        imx=iraf.imexam(Stdout=1)
        foo=open('BlankSky.coo','w')    #Creating blank sky coords files
        i=0
        while i < len(imx) :               
            foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
            i=i+1
        foo.close()
        print('\n Use the ds9 to Mark circles centered at locations to do qphot.')
        print('Enter the center X Y and radius of apperture for qphot annulus and dannulus for sky')
        print('Enter values space seperated in the format below. Enter "q" to exit')
        print('X   Y   Apperture  Annulus  Dannulus ')
        foo=open('qphotinput.txt','w')    #Creating the qphot partameter file
        qphot_inp="junk"
        while (qphot_inp != "q") :
            qphot_inp=raw_input("|> ")
            boolvar=True
            for i in qphot_inp.split() : boolvar = boolvar and is_number(i)
            if boolvar and (len(qphot_inp.split()) == 5) : foo.write(qphot_inp+' \n')
            elif (qphot_inp != "q") : print("Wrong Entry. Please enter properly the 5 values. q is to stop.")
        foo.close()
    #Finished all first images data collection. Now going forward..

    print("\n All required human input of coordinates taken..")
    iraf.cd(backupPWD)


def Createlist_subrout():
    """ Creates the Images4Photo.in containing the image name , filter, exposure time, threshold """
    fooOUT=open(MotherDIR+'/Images4Photo.in','w')
    directories=LoadDirectories(CONF=False)
    Exptime=1  #By default, all TIRSPEC images are scaled to per second.
    for night in directories:
        print('Working on night: '+night)
        try:
            #First we load a dictionary of raw images to their filters
            FiltrFILE=open(night+'/AllObjects.List','r')
            Filtrfiledic=dict([(filtset.split()[0],shlex.split(filtset.rstrip())[1]) for filtset in FiltrFILE.readlines()])  #Dictionary of filterset for each image.
            FiltrFILE.close()
            #Secondly, we load a dictionary of Dither Frame to Raw images
            DitherFILE=open(night+'/FirstoneANDcombinedImages.List','r')
            Ditherfiledic=dict([(ditherset.rstrip().split()[1],ditherset.split()[0]) for ditherset in DitherFILE.readlines() if len(ditherset.split()) == 2])  #Dictionary of First image of each Dither set.
            DitherFILE.close()
            #Now Read and write the images to do photometry one by one.
            ImgsFILE=open(night+'/FirstoneANDalignNcombinedImages.List','r')
        except IOError,e:
            print('Cannot open the image file list.')
            print(e)
            print('So skipping this directory.')
            print('-'*60)
            continue

        for imgline in ImgsFILE.readlines():
            img=imgline.rstrip().split()[1]
            imgfilter=Filtrfiledic[Ditherfiledic[imgline.split()[0]]]
            fooOUT.write(night+'/'+img+'  "'+imgfilter+'"  '+str(Exptime)+'  '+str(threshold)+' \n')
        ImgsFILE.close()
    fooOUT.close()
    print('All nights over...')


def AlignNcombine_subrout(method="average"):
    """ This will align and combine the images in each paragraph in /FirstoneANDcombinedImages.List file for photometry images """
    if TODO == 'S' :
        print("You are doing Spectroscopy. So no align and combining to do on this raw images. \n Skipping...")
        return()
    
    iraf.imcombine.unlearn()
    directories=LoadDirectories(CONF=False)
    for night in directories:
        print('Working on night: '+night)

        #Load all the X,Y coords of star indexed for every file already
        XYFILE=open(night+'/AllObjects2Combine.List','r')
        XYfiledic=dict([(XYset.split()[0],XYset.rstrip().split()[1:]) for XYset in XYFILE.readlines() if len(XYset.split()) == 3 ])  #Dictionary of XY coords for each image.
        XYFILE.close()
        if len(XYfiledic) == 0 : #No images this night..
            print('No images to work on this night. skipping...')
            continue

        Obj2CombFILE=open(night+'/FirstoneANDcombinedImages.List','r')
        #Firstly generate the list of lists of images to combine. Also a dict which maps the combined images to first image.
        ListofLists=[[]]
        Comb2Firstdic=dict()
        for imgline in Obj2CombFILE.readlines():
            if len(imgline.split()) == 0 and len(ListofLists[-1]) != 0 :  ListofLists.append([])  #Start a new list at end
            elif len(imgline.split()) > 0 : 
                ListofLists[-1].append(imgline.rstrip().split()[1]) #Append to the last list
                Comb2Firstdic[imgline.rstrip().split()[1]]=imgline.split()[0] #Adding the mapping to the dictionary
        Obj2CombFILE.close()
        #Now iterate through every list of images to combine
        outlogFILE=open(night+'/FirstoneANDalignNcombinedImages.List','w')
        for imglist in ListofLists:
            if len(imglist) == 1 : #Single image. nothing to align and combine
                OutCombimg=imglist[0]
            elif len(imglist) > 1 :
                OutCombimg=imglist[0][:-5]+'_align'+method+'_'+imglist[-1][:-5]+'.fits'
                OutCoofile=OutCombimg[:-5]+'.GScoo'
                Refimage=imglist[0]
                Xref=eval(XYfiledic[Comb2Firstdic[Refimage]][0])
                Yref=eval(XYfiledic[Comb2Firstdic[Refimage]][1])
                iraf.display(night+'/'+Refimage,1) 
                print ('Press _a_ over some good stars to align, u can press s also, but DONT press r \n')
                imx=iraf.imexam(Stdout=1)
                foo=open(night+'/'+OutCoofile,'w')
                i=2
                while i < len(imx) :               
                    foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                    i=i+2
                foo.close()
                #Now enter the crude shifts for other images from our dic in the file. And also create text files containing images to align and aligned output
                foo2=open(night+'/shifts.in','w')
                alignInpfname=night+'/'+OutCombimg[:-5]+'.ditherList'
                alignOutfname=night+'/'+OutCombimg[:-5]+'.AlignedditherList'
                imgs2align=open(alignInpfname,'w')
                imgs2alignOUT=open(alignOutfname,'w')
                for img in imglist[1:]:
                    Xin=eval(XYfiledic[Comb2Firstdic[img]][0])  
                    Yin=eval(XYfiledic[Comb2Firstdic[img]][1])
                    foo2.write(str(Xref-Xin)+'   '+str(Yref-Yin)+'\n')
                    imgs2align.write(night+'/'+img+'\n')
                    imgs2alignOUT.write(night+'/'+'s'+img+'\n')
                foo2.close()
                imgs2align.close()
                imgs2alignOUT.close()
                try :  #Now align and if succeded combine those images....
                    iraf.imalign(input='@'+alignInpfname, reference=night+'/'+Refimage, coords=night+'/'+OutCoofile, output='@'+alignOutfname, shifts=night+'/shifts.in', interp_type="nearest",boundary_type="constant",trimimages="no")
                    iraf.imcombine(input=night+'/'+Refimage+','+'@'+alignOutfname, output=night+'/'+OutCombimg,combine=method,reject="sigclip")
                except iraf.IrafError, e :
                    print ('IRAF ERROR : Some image might be having problem. Remove it and try later')
                    print e
                    print('-'*60)
                    traceback.print_exc(file=sys.stdout)
                    print('-'*60)
            
            outlogFILE.write(imglist[0]+' '+OutCombimg+'\n')
        outlogFILE.close()
    print('All nights over...')             
                                
                
def FixBadPixels(images,nightdir):
    """ This will run iraf task proto.fixpix to interpolate badpixels """
    if TODO=='P' : PixelMask=nightdir+'/'+PhotBadPixelMaskName
    elif TODO=='S' : PixelMask=nightdir+'/'+SpecBadPixelMaskName
    else : 
        print("What are you doing? (S or P)")
        return
    if not os.path.isfile(PixelMask):
        print("No Bad Pixel Mask file found by the name "+ PixelMask)
        print("Hence skipping Bad pixel interpolation")
        return

    iraf.proto(_doprint=0)
    iraf.fixpix.unlearn()
    iraf.fixpix(images=images,masks=PixelMask)

def CombDith_FlatCorr_subrout(method="median",FlatStatSection='[200:800,200:800]'):
    """ This will combine (default=median) with avsigclip the images in single dither and also create corresponding normalized flats and divide for flat correction """
    iraf.imcombine.unlearn()
    directories=LoadDirectories(CONF=False)
    for night in directories:
        print('Working on night: '+night)
        #Load all the Flat indexing file data
        FlatFILE=open(night+'/AllObjects-FinalFlat.List','r')
        Flatfiledic=dict([(flatset.split()[0],flatset.rstrip().split()[1:]) for flatset in FlatFILE.readlines()])  #Dictionary of flats list for each image.
        FlatFILE.close()
        if len(Flatfiledic) == 0 : #No images this night..
            print('No images to work on this night. skipping...')
            continue
        if TODO=='P':  #Load all the FilterSet indexing file data
            FiltrFILE=open(night+'/AllObjects.List','r')
            Filtrfiledic=dict([(filtset.split()[0],shlex.split(filtset.rstrip())[1]) for filtset in FiltrFILE.readlines()])  #Dictionary of filterset for each image.
            FiltrFILE.close()
            NewFiltSet='(Blah,Blah,Blah)'

        Obj2CombFILE=open(night+'/AllObjects2Combine.List','r')
        #Secondly generate the list of lists of images to combine.
        ListofLists=[[]]
        for imgline in Obj2CombFILE.readlines():
            if len(imgline.split()) == 0 and len(ListofLists[-1]) != 0 :  ListofLists.append([])  #Start a new list at end
            elif len(imgline.split()) > 0 : ListofLists[-1].append(imgline.split()[0]) #Append to the last list
        Obj2CombFILE.close()
        #Now iterate through every list of images to combine
        outlogFILE=open(night+'/FirstoneANDcombinedImages.List','w')
        for imglist in ListofLists:
            if len(imglist) == 1 : #Single image. no need to combine
                OutCombimg=imglist[0]
            elif len(imglist) > 1 :
                OutCombimg=imglist[0][:-5]+'_'+method+'_'+imglist[-1][:-5]+'.fits'
                inpVAR=','.join([night+'/'+img for img in imglist])
                iraf.imcombine(input=inpVAR, output=night+'/'+OutCombimg,combine=method,reject="sigclip")
            #Now make list of flats to be combined for this image set
            Flats2Comb=[]
            for img in imglist:
                Flats2Comb+=Flatfiledic[img]  #Adding all the flat lists
            Flats2Comb=set(Flats2Comb)  #Making a set to remove duplicates
            #Write all these flat names to a file.
            imgflatlistfname=night+'/'+OutCombimg[:-5]+'.flatlist'
            imgflatlistFILE=open(imgflatlistfname,'w')
            imgflatlistFILE.write('\n'.join([night+'/'+fla for fla in Flats2Comb]))
            imgflatlistFILE.close()
            
#To Do Later:: #If spectroscopy of cross disperse mode : change the FlatStatSection
            
            outflatname=night+'/'+OutCombimg[:-5]+'_flat.fits'
            iraf.imcombine (input='@'+imgflatlistfname, output=outflatname, combine="median", scale="median",reject="sigclip", statsec=FlatStatSection)
            statout=iraf.imstatistics(outflatname+FlatStatSection,fields='mode',Stdout=1)
            mode=float(statout[1])
            #We will normalise this flat with the mode of pixels in FlatStatSection
            Noutflatname=night+'/'+OutCombimg[:-5]+'_Nflat.fits'
            iraf.imarith(operand1=outflatname,op="/",operand2=mode,result=Noutflatname)
            #Now divide by flat...
            OutFCimg=OutCombimg[:-5]+'_FC.fits'
            iraf.imarith(operand1=night+'/'+OutCombimg,op="/",operand2=Noutflatname,result=night+'/'+OutFCimg)
            #Now interpolate the bad pixels in the final image.
            FixBadPixels(night+'/'+OutFCimg,night)
            if TODO=='P':
                Oldfiltset=NewFiltSet
                NewFiltSet=Filtrfiledic[imglist[0]]
                if Oldfiltset != NewFiltSet : outlogFILE.write('\n')  #Entering a blank line to show change of filters in Photometry
            # elif TODO=='S':
            #     outlogFILE.write('\n')  #Entering a blank line no matter what. We will ask user to change is they want to move and add.
            outlogFILE.write(imglist[0]+' '+OutFCimg+'\n')
        outlogFILE.close()
        if TODO=='P': print('Edit the spaces between image sets in file '+night+'/FirstoneANDcombinedImages.List'+' to align and combine them in next step.')
    print('All nights over...')             
                
                

def Manual_InspectFlat_subrout():
    """ This will display Flats and Argons one after other, and based on user input select/reject """
    directories=LoadDirectories(CONF=True)
    filelist=['AllObjects-Flat.List']
    outfilelist=['AllObjects-FinalFlat.List']
    if TODO=='S' :
        filelist.append('AllObjects-Argon.List')
        outfilelist.append('AllObjects-FinalArgon.List')
    for night in directories:
        for inpfile,outfile in zip(filelist,outfilelist):
            print('Files in:'+night+'/'+inpfile)
            inFILE=open(night+'/'+inpfile,'r')
            ouFILE=open(night+'/'+outfile,'w')
            AlwaysRemove=[]
            AlwaysAccept=[]
            for inpline in inFILE.readlines():
                inplinelist=inpline.rstrip().split()
                ScienceImg=inplinelist[0]
                CalImgs=[imgs for imgs in inplinelist[1:] if imgs not in AlwaysRemove]
                FinalCalImgs=CalImgs[:]
                print('For the science image: '+ScienceImg)
                for img in CalImgs:
                    if img not in AlwaysAccept:
                        iraf.display(night+'/'+img,1)
                        print(night+'/'+img)
                        verdict=''
                        verdict=raw_input('Enter "r" to reject, "ra" to reject always in future, "aa" to always accept in future:')
                        if verdict== 'r' :
                            FinalCalImgs.remove(img)
                            print("Removing this image : "+img)
                        elif verdict== 'ra' :
                            FinalCalImgs.remove(img)
                            AlwaysRemove.append(img)
                            print("Removing this image forever: "+img)
                        elif verdict== 'aa' :
                            AlwaysAccept.append(img)
                            print("Always accept this image forever this night (Dangerous): "+img)
                #Writing the final surviving calibration files to output file
                ouFILE.write(' '.join([ScienceImg]+FinalCalImgs)+'\n')
            ouFILE.close()
            inFILE.close()
    print('All nights over...') 
               
    

def Manual_InspectObj_subrout():
    """ This will display one image after other, and based on user input classify images of each dither position """
    directories=LoadDirectories(CONF=True)
    print("Press _a_ and then _q_ over a good central star for selecting image")
    for night in directories:
        ObjFILE=open(night+'/AllObjects.List','r')
        Obj2CombFILE=open(night+'/AllObjects2Combine.List','w')
        newX=0
        newY=0
        newU_L_Sfilter='(Blah,Blah,Blah)'
        for objline in ObjFILE.readlines():
            img=objline.split()[0]
            iraf.display(night+'/'+img,1)
            print(objline)
            print('\n To discard this image press _q_ without pressing _a_')
            imx=iraf.imexam(Stdout=1)
            if len(imx) <= 1 : #Then discard this image
                print('Discarding image :'+night+'/'+img)
                continue
            #Else, continue below
            if TODO=='P': #If we were doing photometry
                oldX=newX
                oldY=newY
                oldU_L_Sfilter=newU_L_Sfilter
                newU_L_Sfilter=shlex.split(objline)[1]
                FWHM=eval(imx[3].split()[-1])
                newX=eval(imx[2].split()[0])
                newY=eval(imx[2].split()[1])
                #Print blank enter in the output file if the star has shifted
                StarShifted= np.sqrt((newX-oldX)**2 +(newY-oldY)**2) > 1*FWHM
                # or filter wheels have changed.
                FiltersChanged= newU_L_Sfilter != oldU_L_Sfilter 
                if StarShifted or FiltersChanged : Obj2CombFILE.write('\n')

            elif TODO=='S' : #If doing spectroscopy
                Obj2CombFILE.write('\n')   #Enter a blank line anyway. We can ask user to remove it if he/she needs to combine frames.

            #Now, add this img name to dither image list
            Obj2CombFILE.write(img+' '+str(newX)+' '+str(newY)+'\n')
        Obj2CombFILE.close()
        ObjFILE.close()
        print('We have made the selected list of images in '+night+'/AllObjects2Combine.List  \n Add blank lines between file names to prevent them from median combining. \n Remove the blank line between file names, which you want to combine.')
        os.system(TEXTEDITOR+' '+night+'/AllObjects2Combine.List')
    print('All nights over...')
        
             
def SelectionofFrames_subrout():
    """ Selects the images to reduce and create tables of corresponding Flat and Argons """
    directories=LoadDirectories(CONF=True)
    FiltREdic=dict()
    ArgonREdic=dict()
    ObjRE=raw_input("Enter Regular Expression to select the objects from all dirs: ").strip(' ')
    # For RE rules See: http://docs.python.org/2/howto/regex.html#regex-howto
    regexpObj= re.compile(r''+ObjRE)
    #Generating list of objects frames
    for night in directories:
        print("Working on night : "+night)
        slopeimgFILE=open(night+'/SlopeimagesLog.txt')
        slopeimgFILElines=slopeimgFILE.readlines()
        slopeimgFILE.close()
        ObjList=[imgline.rstrip() for imgline in slopeimgFILElines if regexpObj.search(imgline.split()[0]) is not None ]
        FiltList=set()  #Set to store the list of filters needed to find flat/Argon for
        ObjFILE=open(night+'/AllObjects.List','w')
        for Objline in ObjList:
            Name=shlex.split(Objline)[0]
            U_L_Sfilter=tuple([pos.upper() for pos in shlex.split(Objline)[5:8]])  #(UPPER,LOWER,SLIT) filters in Uppercase
            FiltList.add(U_L_Sfilter)
            ObjFILE.write(Name+'    "'+str(U_L_Sfilter)+'"\n')
        ObjFILE.close()
        #Now ask for flats in each filters
        Flatlistdic=dict()
        for filt in FiltList:
            if filt not in FiltREdic.keys() : 
                FiltREdic[filt]=ObjRE  #Setting default to self
            #Ask user again to confirm or change if he/she needs to
            InpfiltRE=raw_input("Enter Regular Expression for the flat of filters %s (default: %s) : "%(str(filt),FiltREdic[filt])).strip(' ')
            if InpfiltRE :
                FiltREdic[filt]=InpfiltRE
            regexpFilt= re.compile(r''+FiltREdic[filt])
            FlatList=[imgline.split()[0] for imgline in slopeimgFILElines if (regexpFilt.search(imgline.split()[0]) is not None) and (filt == tuple([pos.upper() for pos in shlex.split(imgline)[5:8]])) ]
            Flatlistdic[filt]=FlatList  #Saving flat list for this filter set
        
        #What is happening here is very dumb.. I am taking all pattern matching images with same filter combination.
        #It is maximally dumb thing to do for spectra's Continuum and Argon. We should be selecting only contemporary Continuum flats and Argon.
        #I'm feeling brain dead now to implement it. // To Do later...

        #Now if We are doing Spectroscopy, Find the corresponding Argon lamps also
        if TODO=='S':
            Argonlistdic=dict()
            for filt in FiltList:
                if filt not in ArgonREdic.keys() : 
                    ArgonREdic[filt]='.*Argon.*'  #Setting default to *Argon*
                #Ask user again to confirm or change if he/she needs to
                InpfiltRE=raw_input("Enter Regular Expression for the Argon of filters %s (default: %s) : "%(str(filt),ArgonREdic[filt])).strip(' ')
                if InpfiltRE :
                    ArgonREdic[filt]=InpfiltRE
                regexpArg= re.compile(r''+ArgonREdic[filt])
                ArgonList=[imgline.split()[0] for imgline in slopeimgFILElines if (regexpArg.search(imgline.split()[0]) is not None) and (filt == tuple([pos.upper() for pos in shlex.split(imgline)[5:8]])) ]
                Argonlistdic[filt]=ArgonList  #Saving Argon list for this filter set
            
        #Now, load the Object list and write to a file the Obj and corresponding flats/Argons
        ObjFlatFILE=open(night+'/AllObjects-Flat.List','w')
        if TODO=='S': ObjArgonFILE=open(night+'/AllObjects-Argon.List','w')
        for Objline in ObjList:
            Name=shlex.split(Objline)[0]
            U_L_Sfilter=tuple([pos.upper() for pos in shlex.split(Objline)[5:8]])  #(UPPER,LOWER,SLIT) filters in UPPERCASE
            ObjFlatFILE.write(Name+'  '+' '.join(Flatlistdic[U_L_Sfilter])+'\n')
            if TODO=='S' :ObjArgonFILE.write(Name+'  '+' '.join(Argonlistdic[U_L_Sfilter])+'\n')
        ObjFlatFILE.close()
        print('Edit and save the Flat/Argon list associations for this night :'+night)
        os.system(TEXTEDITOR+' '+night+'/AllObjects-Flat.List')
        if TODO=='S': 
            ObjArgonFILE.close()
            os.system(TEXTEDITOR+' '+night+'/AllObjects-Argon.List')
    print('All nights over...')    
    

def LoadDirectories(CONF=False):
    """ Loads the directores and return the list of directories to do analysis """
    try :
        directoriesF=open(MotherDIR+'/directories','r')
    except IOError :
        #Creating a text file containg the directories to visit if it doesn't already exist
        os.system('find . -type d -maxdepth 1 -mindepth 1 | sort > directories ')
        directoriesF=open(MotherDIR+'/directories','r')
    directories=[dirs.rstrip() for dirs in directoriesF.readlines()]
    directoriesF.close()
    if CONF == True :
        #Ask user again to confirm or change if he/she needs to
        InpList=raw_input('Enter the directories to analyse (default: %s) :'%','.join(directories)).strip(' ')
        if InpList : 
            directories=InpList.split(',')
            directoriesF=open(MotherDIR+'/directories','w') #Updateing directories file
            directoriesF.write('\n'.join(directories))
            directoriesF.close()
    if len(directories) == 0 : 
        print('Atleast one directory containing data to be given as input')
        exit(1)
    else :
        return directories
    

def Backup_subrout():
    """ Copies all the files in present directory to the ../BACKUPDIR """
    os.system('mkdir  ../'+BACKUPDIR)
    print("Copying files to ../"+BACKUPDIR)
    os.system('cp -r * ../'+BACKUPDIR)

#-----Main Program Begins here........

try : 
    configfile=open('TIRSPECscript.conf','r')
except IOError :
    print ("Error: Copy the TIRSPECscript.conf into this directory contianing folders of each night data, before running the script.")
    exit(1)
for con in configfile.readlines():
    con=con.rstrip()
    if len(con.split()) >= 2 :
        if con.split()[0] == "VERBOSE=" :
            VER=con.split()[1]
        elif con.split()[0] == "TODO=" :
            TODO=con.split()[1]
        elif con.split()[0] == "TEXTEDITOR=" :
            TEXTEDITOR=shlex.split(con)[1]
        elif con.split()[0] == "IMGCOMB=" :
            IMGCOMBMETHOD=con.split()[1]
        elif con.split()[0] == "DITHERCOMB=" :
            DITHERCOMBMETHOD=con.split()[1]

        elif con.split()[0] == "BPMPHOTO=" :
            PhotBadPixelMaskName=con.split()[1]
        elif con.split()[0] == "BPMSPEC=" :
            SpecBadPixelMaskName=con.split()[1]
        
        elif con.split()[0] == "THRESHOLD=" :
            threshold=con.split()[1]
        elif con.split()[0] == "EPADU=" :
            EPADU=con.split()[1]
        elif con.split()[0] == "READNOISE=" :
            READNOISE=con.split()[1]
        elif con.split()[0] == "DATAMAX=" :
            DATAMAX=con.split()[1]

        elif con.split()[0] == "APPERTURE=" :
            APPERTURE=con.split()[1]
        elif con.split()[0] == "ANNULUS=" :
            ANNULUS=con.split()[1]
        elif con.split()[0] == "DANNULUS=" :
            DANNULUS=con.split()[1]

        elif con.split()[0] == "EXPTIME=" :
            EXPTIMEHDR=con.split()[1]
        elif con.split()[0] == "FILTER=" :
            FILTERHDR=con.split()[1]
        elif con.split()[0] == "UT=" :
            UTHDR=con.split()[1]
        elif con.split()[0] == "OBJECT=" :
            OBJECTHDR=con.split()[1]
        elif con.split()[0] == "COMMENT=" :
            COMMENTHDR=con.split()[1]

        elif con.split()[0] == "OUTPUT=" :
            OUTPUTfile=con.split()[1]
        elif con.split()[0] == "BACKUP=" :
            BACKUPDIR=con.split()[1]

        elif con.split()[0] == "CONVOLVEIMG=" :
            CONVOLVEIMG=con.split()[1]
        elif con.split()[0] == "DOPSF=" :
            DOPSF=con.split()[1]

configfile.close()
MotherDIR=os.getcwd()
#    OUTPUTfilePATH=MotherDIR+'/'+OUTPUTfile
parentdir=MotherDIR.split('/')[-1]

if TODO == 'P' : todoinwords='Photometry'
elif TODO == 'S' : todoinwords='Spectroscopy'
 
print("Very Very Important: Backup your files first. Don't proceed without backup.\n")
print(" ---------------- Welcome to TIRSPEC \033[91m "+todoinwords+" \033[0m Script --------------- \n")
print("Enter the Serial numbers (space seperated if more than one task in succession) \n")
print("0  Backup files in current directory to ../"+BACKUPDIR+"\n")
print("1  Selection of object frames to reduce \n")
print("2  Manually inspect and reject object images by displaying one by one to classify \n")
print("3  Manually inspect and reject Flats/Argons by displaying one by one\n")
print("4  Combine images in a Dither, apply Flat Correction and Bad pixel interpolation\n")
if TODO=='P':print("5  Align and combine combined images of each Dither in Photometry data \n")
if TODO=='P':print("6  Make the list of images, Images4Photo.in to do Photometry \n")
if TODO=='P':print("7  Select Stars and Sky region of the field on first image \n")
#print("7  Remove Cosmic Rays on all the images in Images4Photo.in. IMP:It will OVERWRITE original images.\n")
if TODO=='P':print("8  Create Sextracter config file & coordinate output of first image in this directory \n")
if TODO=='P':print("9  Do Photometry \n")
print("--------------------------------------------------------------- \n")
todo=raw_input('Enter the list : ')
todo=todo.split()
if ("2" in todo) or ("3" in todo) or ("4" in todo) or ("5" in todo) or ("7" in todo) or ("8" in todo) or ("9" in todo) :
    from pyraf import iraf
if ("9" in todo) :
    from astropy.io import ascii
    import astropy.table as table  #Requires astropy version >= 0.3 for vstack function
for task in todo :
    if task == "0" :
        Backup_subrout()
    elif task == "1" :
        SelectionofFrames_subrout()
    elif task == "2" :
        Manual_InspectObj_subrout()
    elif task == "3" :
        Manual_InspectFlat_subrout()
    elif task == "4" :
        CombDith_FlatCorr_subrout(method=IMGCOMBMETHOD)
    elif task == "5" :
        if TODO=='P': AlignNcombine_subrout(method=DITHERCOMBMETHOD)

    elif task == "6" :
        if TODO=='P': Createlist_subrout()

    elif task == "7" :
        if TODO=='P': Star_sky_subrout()

    elif task == "8" :
        if TODO=='P': Sextractor_subrout()

    elif task == "9" : 
        if TODO=='P': Photometry()
print("All tasks over....Enjoy!!!_________indiajoe@gmail.com")
            

