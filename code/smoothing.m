function smooth_mu = smoothing(mu,Smooth_Operator,Operator)
% Smoothing the Beltrami coefficient.
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2014-2022, Gary Pui-Tung Choi
% https://math.mit.edu/~ptchoi/

vmu = Operator.f2v*mu;
nvmu = Smooth_Operator\abs(vmu);
vmu = nvmu.*(cos(angle(vmu))+sqrt(-1)*sin(angle(vmu)));
smooth_mu = Operator.v2f*vmu;

end