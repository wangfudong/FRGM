# A Functional Representation for Graph Matching
The implementation of our work FRGM ([[project page](http://captain.whu.edu.cn/FRGM/)] [[arxiv](https://arxiv.org/abs/1901.05179)]).

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
@article{FRGM_2019,
  author    = {Fudong Wang and
               Gui{-}Song Xia and
               Nan Xue and
               Yipeng Zhang and
               Marcello Pelillo},
  title     = {A Functional Representation for Graph Matching},
  year      = {2019},
  url       = {http://arxiv.org/abs/1901.05179},
  archivePrefix = {arXiv},
  eprint    = {1901.05179},
}
```
