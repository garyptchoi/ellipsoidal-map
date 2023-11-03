% Fast ellipsoidal quasi-conformal map (FEQCM) for genus-0 closed surfaces
%
% Compute an ellipsoidal quasi-conformal parameterization of a genus-0
% closed surface with prescribed landmark constraints using the method in [1].
%
% Usage:
% map = ellipsoidal_quasiconformal_map(v,f,a,b,c,landmark,target,lambda)
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 triangle mesh
% f: nf x 3 triangulations of a genus-0 triangle mesh
% a,b,c: the elliptic radii of the target ellipsoid
% landmark: k x 1 vertex indices of the landmarks
% target: k x 3 coordinates of the target positions on the ellipsoid
% lambda: a positive balancing factor for landmark matching
%
% Output:
% map: nv x 3 vertex coordinates of the ellipsoidal quasi-conformal parameterization
%
% Remarks:
% - The input surface mesh is assumed to be optimally aligned with the 
%   x,y,z-axes beforehand.
% - the "target" coordinates must lie on the prescribed ellipsoid.
%
% If you use this code in your own work, please cite the following papers:
%
% [1] G. P. T. Choi, 
%     "Fast ellipsoidal conformal and quasi-conformal parameterization of genus-0 closed surfaces".
%     Preprint, 2023.
% 
% Copyright (c) 2023, Gary P. T. Choi

addpath('code');
addpath('data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 1: Buddha
load('buddha_with_landmark.mat');

plot_mesh(v,f);
hold on;
plot3(v(landmark,1),v(landmark,2),v(landmark,3),...
    'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',15);
view([-90 15]);
title('Input surface');

%% ellipsoidal conformal map (landmarks not aligned)
map_fecm = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map_fecm,f);
hold on;
plot3(map_fecm(landmark,1),map_fecm(landmark,2),map_fecm(landmark,3),...
    'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',15);
plot3(target(:,1),target(:,2),target(:,3),...
    'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',15);
view([-90 15]);
title('Ellipsoidal conformal parameterization');

%% ellipsoidal quasi-conformal map (landmarks aligned)
map_feqcm = ellipsoidal_quasiconformal_map(v,f,a,b,c,landmark,target,5);

plot_mesh(map_feqcm,f);
hold on;
plot3(map_feqcm(landmark,1),map_feqcm(landmark,2),map_feqcm(landmark,3),...
    'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',15);
plot3(target(:,1),target(:,2),target(:,3),...
    'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',15);
view([-90 15]);
title('Ellipsoidal quasi-conformal parameterization');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 2: Bulldog
load('bulldog_with_landmark.mat');

plot_mesh(v,f);
hold on;
plot3(v(landmark,1),v(landmark,2),v(landmark,3),...
    'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',15);
view([35 15]);
title('Input surface');

%% ellipsoidal conformal map (landmarks not aligned)
map_fecm = ellipsoidal_conformal_map(v,f,a,b,c);

plot_mesh(map_fecm,f);
hold on;
plot3(target(:,1),target(:,2),target(:,3),...
    'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',15);
plot3(map_fecm(landmark,1),map_fecm(landmark,2),map_fecm(landmark,3),...
    'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',15);
view([15 10]);
title('Ellipsoidal conformal parameterization');

%% ellipsoidal quasi-conformal map (landmarks aligned)
map_feqcm = ellipsoidal_quasiconformal_map(v,f,a,b,c,landmark,target,5);

plot_mesh(map_feqcm,f);
hold on;
plot3(target(:,1),target(:,2),target(:,3),...
    'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',15);
plot3(map_feqcm(landmark,1),map_feqcm(landmark,2),map_feqcm(landmark,3),...
    'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',15);
view([15 10]);
title('Ellipsoidal quasi-conformal parameterization');
