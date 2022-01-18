# SOFA-object-for-Pure-data
SOFA-object-for-Pure-data is a graphical programming language called Pure-data that uses a Spatial Oriented Format (SOFA) file to convolve a user-specified .wav file with IR (impulse response) based on elevation azimuth and distance position information.
It is an object that performs IR (impulse response) convolution according to the user-specified .wav file and elevation azimuth and distance position information, and outputs the result.
With this object, you can manipulate stereophonic sound in real time within Pure-data.

//Sakurai adds points
I added the ability to change the azimuth of the speaker to the existing code.
The existing code could only read one SOFA file, so I made it so that it could read multiple files.
I also made it possible to change multiple SOFA files read by the slider, just as the existing code uses the slider to change the azimuth of the listener.
# Requirement

* Pure-data
* fftw-3.3.10
* libmysofa



# Installation
Pure-data
http://puredata.info/downloads/pure-data
The site to install Puredata.
Please install Puredata according to your OS here.

fftw
http://www.fftw.org/download.html
Install the latest fftw here

> $ tar -zxvf fftw-3.3.9.tar.gz

> $ cd fftw-3.3.9

> $ make

> $ sudo make install

Installation is complete.

Make sure that libfftw3.a, libfftw3.la, libfftw3f.a, libfftw3f.la, and libmysofa.a are in ~/usr/local/lib/.

If you don't have them, move them to ~/usr/local/lib.


libmysofa
https://github.com/hoene/libmysofa
Follow the Compile here.

If you do not have Cunit installed, please install it.
https://sourceforge.net/projects/cunit/files/

Also, if you do not have Node.js installed, please install it.


# Usage
Execute this command
Put the made file in the path of pure-data.
Compilation is complete.

How to use pdsofa~help.pd

> $ cd mysofa~

> $ make
(Don't forget [> $ make clean] if darwin file is in the current directory.)
First, select the .wav file you want to put in the sound file array.

Turn on the DSP.

Press the bang connected to [tabplay~].

Change the speaker azimuth and orientation values, and the sound will be output accordingly.

# Author
* Author Yuki Sakurai
* E-mail s1260177@u-aizu.ac.jp
