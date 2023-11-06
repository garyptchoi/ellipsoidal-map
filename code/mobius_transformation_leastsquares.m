function fz = mobius_transformation_leastsquares(z, p, q)

% If you use this code in your own work, please cite the following papers:
%
% [1] G. P. T. Choi, 
%     "Fast ellipsoidal conformal and quasi-conformal parameterization of genus-0 closed surfaces".
%     Preprint, arXiv:2311.01788, 2023.
% 
% Copyright (c) 2023, Gary P. T. Choi

if size(p,2) == 2
    p = p(:,1)+1i*p(:,2);
end

if size(q,2) == 2
    q = q(:,1)+1i*q(:,2);
end

% only consider r*exp(1i*t)*z
d = @(x) sum(abs(x(1)*exp(1i*x(2))*p - q).^2);

% Optimization setup
x0 = [1,0]; % initial guess
lb = [-100,-4*pi]; % lower bound for the parameters
ub = [100,4*pi]; % upper bound for the parameters
options = optimoptions('fmincon','Display','off');

% Optimization (may further supply gradients for better result, not yet implemented)
x = fmincon(d,x0,[],[],[],[],lb,ub,[],options);

% obtain the conformal parameterization with area distortion corrected
fz = x(1)*exp(1i*x(2))*z;
