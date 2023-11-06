% Fast ellipsoidal conformal map (FECM) for genus-0 closed surfaces
%
% Compute an ellipsoidal conformal parameterization of a genus-0 closed
% surface using the method in [1].
%
% Usage:
% map = ellipsoidal_conformal_map(v,f,a,b,c)
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 triangle mesh
% f: nf x 3 triangulations of a genus-0 triangle mesh
% a,b,c: the elliptic radii of the target ellipsoid
%
% Output:
% map: nv x 3 vertex coordinates of the ellipsoidal conformal parameterization
%
% Remark:
% - The input surface mesh is assumed to be optimally aligned with the 
%   x,y,z-axes beforehand.
%
% If you use this code in your own work, please cite the following papers:
%
% [1] G. P. T. Choi, 
%     "Fast ellipsoidal conformal and quasi-conformal parameterization of genus-0 closed surfaces".
%     Preprint, arXiv:2311.01788, 2023.
% 
% Copyright (c) 2023, Gary P. T. Choi

addpath('code');
addpath('data');

%% Example 1: Hippocampus
load('hippocampus.mat');

plot_mesh(v,f);
view([-90 10]);
title('Input surface');

map = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map,f);
view([-30 5]);
title('Ellipsoidal conformal parameterization');

%% Example 2: Buddha
load('buddha.mat');

plot_mesh(v,f);
view([-70 10]);
title('Input surface');

map = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map,f);
view([-70 10]);
title('Ellipsoidal conformal parameterization');

%% Example 3: Bimba
load('bimba.mat');

plot_mesh(v,f);
view([90 0]);
title('Input surface');

map = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map,f);
view([90 0]);
title('Ellipsoidal conformal parameterization');

%% Example 4: Lion Vase
load('vaselion.mat');

plot_mesh(v,f);
view([90 0]);
title('Input surface');

map = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map,f);
view([90 0]);
title('Ellipsoidal conformal parameterization');

%% Example 5: Moai
load('moai.mat');

plot_mesh(v,f);
view([80 10]);
title('Input surface');

map = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map,f);
view([80 10]);
title('Ellipsoidal conformal parameterization');

%% Example 6: Bulldog
load('bulldog.mat');

plot_mesh(v,f);
view([10 10]);
title('Input surface');

map = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map,f);
view([-25 10]);
title('Ellipsoidal conformal parameterization');


%% Example 7: Nefertiti
load('nefertiti.mat');

plot_mesh(v,f);
view([10 10]);
title('Input surface');

map = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map,f);
view([-25 10]);
title('Ellipsoidal conformal parameterization');