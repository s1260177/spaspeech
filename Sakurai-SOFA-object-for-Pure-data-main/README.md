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
If you have not downloaded the following, please do so. (If you get an error when compiling, please make sure you have all of the following items)

Pure-data
http://puredata.info/downloads/pure-data
The site to install Puredata.
Please install Puredata according to your OS here.

Homebrew
https://brew.sh/index_ja

$ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

Brew(cunit)
https://formulae.brew.sh/formula/cunit
$ brew install cunit

Brew(node.js)
https://qiita.com/kyosuke5_20/items/c5f68fc9d89b84c0df09
$ brew install nodebrew
$ nodebrew setup
$ nodebrew install-binary latest
$ nodebrew use ***

Brew(cmake)
$ brew install cmake

fftw
http://www.fftw.org/download.html 
$ tar -zxvf fftw-3.3.10.tar.gz
$ cd fftw-3.3.10

fftw3f
$ ./configure --enable-float --enable-shared
$ make
$ sudo make install

libmysofa
https://github.com/hoene/libmysofa
Follow the Compile here.

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
