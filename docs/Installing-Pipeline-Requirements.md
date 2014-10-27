# Installation

TIRSPEC data reduction pipeline heavily depends on standard astronomy python modules and IRAF. This page will help you install all these dependencies required to run pipeline.

If you already have installed the python modules mentioned below. You only have to update them to latest version. Which can be done by the command 

`pip install -U YourModuleNameHere`

# All in one package: Ureka
You could try installing all the required pakages together by installing the Ureka bundle provided by stcsi  http://ssb.stsci.edu/ureka/

But many modules inside it are old. so after installtion, you will have to give the following commands to upgrade all the packages

* `pip install -U numpy`
* `pip install -U scipy`
* `pip install -U astropy`
* `pip install -U matplotlib`
* `pip install -U ipython[all]`

You will still have to install sextractor separately as mentioned below.


# Installing one by one manually.
If you are not using Ureka, you could install every required module manually as explained below. It is not difficult and as a bonus you will also learn how to install software in your gnu/linux machine.

### IRAF (version 2.16 or greater)
IRAF's latest version installation is very easy. (You don't have to create a separate user iraf or anything like that, one has to do for installing 2.14 or older versions)

All you have to do is to download the latest tar.gz from the iraf webpage.
Copy it to `/iraf/iraf/` directory.
Extract and run `./install` script as superuser.

See the page below for detailed example of how to do it.

http://www.astronomy.ohio-state.edu/~khan/iraf/iraf_step_by_step_installation

PS: Beware of old webpages which provides instruction to install older versions of iraf like 2.14. They are very different, and not to be used now.

PPS: copy your login.cl to `~/iraf/` directory. Then you will be able to start pyraf from any directory in your computer.

PPPS: Don't forget to change the following line in login.cl as shown below.

`set     stdimage        = imt1024`

This is because TIRSPEC images are 1024x1024.

### Sextractor
Sextractor is a powerful tool for automated star extraction from the image.
Download the latest version from http://www.astromatic.net/software/sextractor

Installation is straight forward and simply like any other standard gnu/linux packages.

Uncompress the tar.gz.

cd to the uncompressed directory.

`./configure`

`make`

`sudo make install`

If you are getting stuck at an error related to fftw in the `./configure` step, and if you don't want to use some sextractor's new model fitting features, you can skip that feature by giving the following command, instead of simple `./configure`

`./configure --disable-model-fitting`

 This feature is anyway not needed for our pipeline.

PS: Do not install sextractor from ubuntu's repository. Their sextractor's binary has a bug!

## Python modules
Python will the there most likely by default in all gnu/linux computers. So you have to install only the required modules.

pip is the best repository to install latest version of python modules.

Hence, if pip command is missing in your terminal. Install pip first.

`sudo apt-get install python-pip   # On Ubuntu/Debian`

 OR

`sudo yum install python-pip   #On Fedora `

Once you have pip ready. you can install any python modules using it.
You can upgrade pip using pip itself

`sudo pip install -U pip`

PS: Never install python modules using ubuntu/fedora's repository. They are older than the age of the universe!

### numpy
Numpy is the first module we will update

`sudo pip install -U numpy`

### scipy
Scipy is another very useful module to have.

`sudo pip install -U scipy`

This might give some error saying liblas or libatlas is missing.
you can install them using synaptic manager in ubuntu. Remember to install the versions of libblas and libatlas with -dev prefixed to it. (eg: libatlas-base-dev)

After that retrying the scipy installation command.

### ipython
Ipython is a very neat and useful terminal for interactive use of python. It is not needed for the pipeline. But please install it, you will thank me later.

`sudo pip install -U ipython[all]`

### astropy
This is a very important astronomy module used heavily in pipeline. It is important to install the latest version.

`sudo pip install -U astropy`

### pyraf
This module is needed to call IRAF binaries from python.

`sudo pip install -U pyraf`

### matplotlib
This is a very powerful plotting library in python. Before installing if you want to run TIRSPEC's Exposure time calculator, you should also install tcl tk

`sudo yum install tcl tk   # Only for Fedora  Users`

Ubuntu users please install them from your synaptic manager.

Once you have tcl/tk installed. you can install matplotlib using pip. This will combile matplotlib with Tk backend support.

`sudo pip install -U matplotlib`


That is all. now you can start using the scripts in this repository.
