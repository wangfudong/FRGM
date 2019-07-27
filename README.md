# A Functional Representation for Graph Matching
The implementation of our work FRGM ([[project page](http://captain.whu.edu.cn/FRGM/)] [[TPAMI](https://ieeexplore.ieee.org/document/8723156)]).

# Introduction
This work presents a functional representation for graph matching (FRGM). From the functional representation perspective, the matching between graphs can be reformulated as a linear functional between the function spaces of graphs for general graph matching. Moreover, the linear functional representation map can be viewed as a new parameterization for Euclidean graph matching, which allows us to estimate the geometric parameters and correspondence matrix simultaneously.

<p align="center">
<img src="images/frgm-1.svg" width="500">
<p>

# Usage
The structure is organized as follows:
```
your_dir/
  -3rd_party
  -data
  -FRGM-D
  -FRGM-E
  -FRGM-G
  -GM_methods
  -PR_methods
```
The 3rd_party consists of some dependent codes (Shape Context, geodesic, linear assignment, etc) and can be downlowed [here](http://captain.whu.edu.cn/FRGM/code/3rd_party.zip). 
The data can be downloaded [here](http://captain.whu.edu.cn/FRGM/data/data.zip).

The [GM_methods](http://captain.whu.edu.cn/FRGM/code/GM_methods.zip) and [PR_methods](http://captain.whu.edu.cn/FRGM/code/PR_methods.zip) consists of the implementations of the compared methods on general graph matching and Euclidean graph matching with geometric deformations, respectively.

# Citation
If you find our work useful in your research, please consider citing:
```
@ARTICLE{8723156, 
author={F. {Wang} and N. {Xue} and Y. {Zhang} and G. {Xia} and M. {Pelillo}}, 
journal={IEEE Transactions on Pattern Analysis and Machine Intelligence}, 
title={A Functional Representation for Graph Matching}, 
year={2019}, 
volume={}, 
number={}, 
pages={1-1}, 
doi={10.1109/TPAMI.2019.2919308}, 
ISSN={0162-8828}, 
month={},}
```
