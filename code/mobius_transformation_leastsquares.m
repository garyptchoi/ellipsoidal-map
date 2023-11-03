function fz = mobius_transformation_leastsquares(z, p, q)

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

% Most general:
% (az+b)/(cz+d)
%
% d = @(x) sum(abs(((x(1)+x(2)*1i)*p(isfinite(p)&isfinite(q))+(x(3)+x(4)*1i))./...
%     ((x(5)+x(6)*1i)*p(isfinite(p)&isfinite(q))+(x(7)+x(8)*1i))...
%     - q((isfinite(p)&isfinite(q)))).^2);
% 
% % Optimization setup
% x0 = [1,0,0,0,0,0,1,0]; % initial guess
% lb = [-1,-1,-1,-1,-1,-1,-1,-1]*100; % lower bound for the parameters
% ub = [1,1,1,1,1,1,1,1]*100; % upper bound for the parameters
% options = optimoptions('fmincon','Display','off');
% 
% % Optimization (may further supply gradients for better result, not yet implemented)
% x = fmincon(d,x0,[],[],[],[],lb,ub,[],options);
% 
% % obtain the conformal parameterization with area distortion corrected
% fz = ((x(1)+x(2)*1i)*z+(x(3)+x(4)*1i))./((x(5)+x(6)*1i)*z+(x(7)+x(8)*1i));
% 
