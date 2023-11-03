# Fast Ellipsoidal Conformal and Quasi-Conformal Map

<img src = "https://github.com/garyptchoi/ellipsoidal-map/blob/main/cover.png" height="360" />

* **Fast Ellipsoidal Conformal Map (FECM)**: Compute an ellipsoidal conformal parameterization of a genus-0 closed surface using the method in [1].

* **Fast Ellipsoidal Quasi-Conformal Map (FEQCM)**: Compute an ellipsoidal quasi-conformal parameterization of a genus-0 closed surface with prescribed landmark constraints using the method in [1].

Any comments and suggestions are welcome. 

If you use this code in your own work, please cite the following papers:

[1] G. P. T. Choi, 
    "[Fast ellipsoidal conformal and quasi-conformal parameterization of genus-0 closed surfaces.](https://arxiv.org)"
    Preprint, 2023.

Copyright (c) 2023, Gary P. T. Choi

https://www.math.cuhk.edu.hk/~ptchoi/

===============================================================

Usage:
* `map = ellipsoidal_conformal_map(v,f,a,b,c)`
* `map = ellipsoidal_quasiconformal_map(v,f,a,b,c,landmark,target,lambda)`

Input:
* `v`: nv x 3 vertex coordinates of a genus-0 triangle mesh
* `f`: nf x 3 triangulations of a genus-0 triangle mesh
* `a,b,c`: the elliptic radii of the target ellipsoid
* `landmark`: k x 1 vertex indices of the landmarks
* `target`: k x 3 coordinates of the target positions on the ellipsoid
* `lambda`: a positive balancing factor for landmark matching

Output:
* `map`: nv x 3 vertex coordinates of the ellipsoidal conformal/quasi-conformal parameterization

================================================================

Remarks:
* The input surface mesh `(v,f)` is assumed to be optimally aligned with the x,y,z-axes beforehand.
* For `ellipsoidal_quasiconformal_map`, the `target` coordinates must lie on the prescribed ellipsoid. 
