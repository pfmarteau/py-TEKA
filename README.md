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
<p float="left">
  <img src="fig/CBF_ITEKA_15_c.jpg" width="200" height="150">
  <img src="fig/CBF_ITEKA_15_b.jpg" width="200" height="150">
  <img src="fig/CBF_ITEKA_15_f.jpg" width="200" height="150">
  <img src="fig/CBF_ITEKA_Centroids.jpg" width="200" height="150">
</p>
<p float="left">
  <img src="fig/CBF_ITEKA_15_c0.jpg" width="200" height="150">
  <img src="fig/CBF_ITEKA_15_b0.jpg" width="200" height="150">
  <img src="fig/CBF_ITEKA_15_f0.jpg" width="200" height="150">
  <img src="fig/CBF_ITEKA_15_b_t0.jpg" width="200" height="150">
</p>

Please cite the following article if needed to.
[1] Marteau, P.F., Times series averaging and denoising from a probabilistic perspective on time-elastic kernels International Journal of Applied Mathematics and Computer Science, Vol 29, num 2, pp 375–392, De Gruyter editor, 2016 [https://arxiv.org/abs/1611.09194], [bibtex](ref/marteau2016.bib)
