# TEKA-py : Python wrapper to module TEKA 
### Implementation of the TEKA (Time Elastic Kernel Averaging) code defined and used in [1]. The TEKA code is written in C++ and the python wrapper is using the CYTHON compiler

## Requirements
g++ compiler
CYTHON and python3.*

## To install
$ sh install.sh

## To uninstall 
$ sh uninstall.sh

## To test TEKA-py on the Cylinder, Bell, Funnel dataset
$ python3 testCBF.py

![Cylender,Bell,Funnel centroids](CBF_ITEKA_Centroids.jpg?raw=true)


Please cite the following article if needed to.
[1] Marteau, P.F., Times series averaging and denoising from a probabilistic perspective on time-elastic kernels International Journal of Applied Mathematics and Computer Science, De Gruyter pdf
