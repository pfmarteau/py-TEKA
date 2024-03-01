# py-TEKA : Python wrapper to module TEKA 
### Implementation of the TEKA (Time Elastic Kernel Averaging of set of time series) code defined and used in [1]. The TEKA code is written in C++ and the python wrapper is using the CYTHON compiler

## Requirements
- g++ compiler
- CYTHON compiler
- python3.*
- matplotlib

## To install
$ sh install.sh

## To uninstall 
$ sh uninstall.sh

## To test py-TEKA on the Cylinder, Bell, Funnel dataset
$ python3 testCBF.py

<img src="CBF_ITEKA_Centroids.jpg" width="100" height="100">
<img src="CBF_ITEKA_15_c.jpg" width="100" height="100">
<img src="CBF_ITEKA_15_b.jpg" width="100" height="100">
<img src="CBF_ITEKA_15_f.jpg" width="100" height="100">



Please cite the following article if needed to.
[1] Marteau, P.F., Times series averaging and denoising from a probabilistic perspective on time-elastic kernels International Journal of Applied Mathematics and Computer Science, De Gruyter pdf
