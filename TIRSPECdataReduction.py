#!/usr/bin/env python
#This script is to semi-automate basic data reduction of TIRSPEC data.
# IMP:  Keep ds9 open 
#---------------------------------------------indiajoe
import os
import os.path
import glob
import pyfits
import sys, traceback 
import numpy as np
import warnings
import re
import shlex
import readline

# Other modules inserted inside the code are
# pyraf.iraf

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
print("4  Combine images in a Dither and apply Flat Correction \n")
print("5  Align and combine combined images of each Dither in Photometry data \n")
#print("5  Make the list of images, Images4Photo.in to do Photometry \n")
print("6  Select Stars and Sky region of the field on first image \n")
print("7  Remove Cosmic Rays on all the images in Images4Photo.in. IMP:It will OVERWRITE original images.\n")
print("8  Create Sextracter config file & coordinate output of first image in this directory \n")
print("9  Do Photometry \n")
print("--------------------------------------------------------------- \n")
todo=raw_input('Enter the list : ')
todo=todo.split()
if ("2" in todo) or ("3" in todo) or ("4" in todo) or ("5" in todo) or ("6" in todo) or ("7" in todo) or ("8" in todo) or ("9" in todo) :
    from pyraf import iraf
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
        AlignNcombine_subrout(method=DITHERCOMBMETHOD)
#        Createlist_subrout()

    elif task == "6" :
        print('Not implemented')
        exit(1)
        Star_sky_subrout()
    elif task == "7" :
        print('Not implemented')
        exit(1)
        Cosmicrays_subrout()
    elif task == "8" :
        print('Not implemented')
        exit(1)
        Sextractor_subrout()
    elif task == "9" : 
        print('Not implemented')
        exit(1)
        Photometry()
print("All tasks over....Enjoy!!!_________indiajoe@gmail.com")
            

