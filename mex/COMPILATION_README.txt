

Given a working MEX setup, the mex functions can be compiled from MATLAB using the provided build scripts in the respective projects folder.
The "build_CPlusPlus.m" script builds all components automatically and moves the files to the right folders for TrackNTrace to work.


Compiling the TNT fitter needs a few extra steps to compile the ceres library before compiling the mex file:

- FastPsfFitting -

(x) Execute "git submodule update --init --recursive" which gets the FastPsfFittingCode. 
    For compiling ceres we use Tal Ben-Nun's github repository for windows (we use it for linux as well),
    which is also downloaded with the previous command.
(x) Go to the FastPsfFitting directory
(x) Download the eigen library from
     URL: http://eigen.tuxfamily.org/ 
    and put it in "ceres-windows/eigen"

Windows:
To build ceres with Visual Studio simply open the project files in ceres-windows and 
start the compilation of project "ceres" in "Release" "x64" mode.

Linux (tested with Kubuntu):
With Linux ceres can be compiled with CMake as follows:
- Execute linux_install_forCeresBuild.sh to install neccessary libraries
- Switch to the "ceres-windows/ceres-solver/" directory
- Execute "mkdir build && cd build"
- Execute "cmake-gui ..", press "Configure" once, set "BUILD_SHARED_LIBS" to true, press "Generate"
- Close cmake-gui and execute "make -j4"







