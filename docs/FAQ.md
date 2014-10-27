[FAQuesitonLists]: https://github.com/indiajoe/TIRSPEC/wiki/FAQ#frequently-asked-questions-faq "List of FAQs"
## Frequently Asked Questions (FAQ) ##

1. [How to give Regular Expression input for the Flats/Argon?](#1-how-to-give-regular-expression-input-for-the-flatsargon)
2. [Pros and Cons of Cross Dispersed spectrum](#2-pros-and-cons-of-cross-dispersed-spectrum)
3. [Should we edit login.cl before running the pipeline?](#3-should-we-edit-logincl-before-running-the-pipeline)
4. [How do I rename some raw images before reduction?](#4-how-do-i-rename-some-raw-images-before-reduction)
5. [How do I specify my favourite text editor?](#5-how-do-i-specify-my-favourite-text-editor)
6. [Why do we have to give `.*filename.*` instead of `*filename*` in regexp?](#6-why-do-we-have-to-give-filename-instead-of-filename-in-regexp)

[flatregex]: https://github.com/indiajoe/TIRSPEC/wiki/FAQ#1-how-to-give-regular-expression-input-for-the-flatsargon  "Flat Regex"
### 1. How to give Regular Expression input for the Flats/Argon? ###
A typical input request by the pipeline for selecting flat is shown below.
> "Enter regular Expression for the flat of filters ('JOS', 'G','L2')(default: .\*Continu.\*)"

It is asking user for the regular expression to find out the flats for this filter setup (JOS,G,L2).
Typically we name spectroscopic flats with "Continuum" filename. (Argon lamp spectra with "Argon" filename).
Hence when you give `.*Continu.*` to select filenames which match that regular expression.
The code will pick all the files observed in same filter setup which has "Continu" string in it.

But usually we don't want Continuum flats taken for other objects in same filter also to be included. That is why there is an option to give a range of filenumbers. Open SlopeimagesLog.txt of that night and find out the range of filenumbers correspond to this source. After that, you can give the lower limit and upper limit of filenumbers along with `.*Continu.*`

For example   `.*Continu.* 220 330`
Will select all the files which has the string "Continu" in it, and filenumber in the range 220 to 330, and also same filter configuration.

[Back to Top](#frequently-asked-questions-faq)
[XdispSpec]: https://github.com/indiajoe/TIRSPEC/wiki/FAQ#2-pros-and-cons-of-cross-dispersed-spectrum  "XDispersed Mode"
### 2. Pros and Cons of Cross Dispersed spectrum  ###
You have to use short slits for YJ and HK cross dispersed mode. (slit names starting with S).

Advantage : 

1. You get better wavelength coverage. See the wavelength coverage table in the webpage.
2. Simultaneous coverage of both orders of spectrum.

Disadvantage: 

1. Net throughput is less by a factor of 3. So one will have to give more exposure time for good S/N ratio than single order mode. If magnitude is less than ~7, then S/N will not be a problem.

[Back to Top](#frequently-asked-questions-faq)

[logincl]: https://github.com/indiajoe/TIRSPEC/wiki/FAQ#3-should-we-edit-logincl-before-running-the-pipeline  "login.cl file"
### 3. Should we edit login.cl before running the pipeline ###
**Yes**.
login.cl file should be there in iraf directory of your home directory. ie`~/iraf/login.cl`
You should also remove the \# of the following line and set it to imt1024.

`set     stdimage        = imt1024 `

This is due to a bug in ds9 - iraf v2.16 interaction. You get correct coordinates with imexam only when this image size setting is given in login.cl.
 
[Back to Top](#frequently-asked-questions-faq)

[renaming]: https://github.com/indiajoe/TIRSPEC/wiki/FAQ#4-how-do-i-rename-some-raw-images-before-reduction  "Renaming files"
### 4. How do I rename some raw images before reduction? ###
Just renaming the filenames alone will not do. One has to update the SlopeimagesLog.txt file also.
For doing safe renaming in bulk, cd to the night's directory and run `RenameImgs.sh`.

This will open the SlopeimagesLog.txt in a text editor in which you can edit the filenames and save and close. The script will rename the edited filenames accordingly.

Remember to backup your data before doing operations like renaming.

[Back to Top](#frequently-asked-questions-faq)

[texteditor]: https://github.com/indiajoe/TIRSPEC/wiki/FAQ#5-how-do-i-specify-my-favourite-text-editor  "Text Editor"
### 5. How do I specify my favourite text editor?###
It is important that whatever text editor you use, it should be opening in blocking mode. So that the execution of the script will be paused till you save and close the editor.

Following are the ways you can specify some of the commonly used text editors.
* Emacs : `"emacs -nw" `
* vi    : `"vi"`
* gedit : `"gedit -w"`
* nano  : `"nano"`
 
[Back to Top](#frequently-asked-questions-faq)
[pyregex]:https://github.com/indiajoe/TIRSPEC/wiki/FAQ#6-why-do-we-have-to-give-filename-instead-of-filename-in-regexp "Python regex"
### 6. Why do we have to give `.*filename.*` instead of `*filename*` in regexp?###
`.*` is the standard regular expression which python follows for matching any character multiple number of times. It is equivalent to `*` in shell terminal. Hence in python, instead of `*filename*` one should follow the standard regexp, ie. `.*filename.*` to match filenames containing the string `filename` in it.

For more details on other features of python regexp see [PythonReDoc](https://docs.python.org/2/library/re.html), [ReHowTo](https://docs.python.org/2/howto/regex.html).

You can use the power of Regexp to select any complicated pattern of filenames you want in pipeline.

[Back to Top](#frequently-asked-questions-faq)