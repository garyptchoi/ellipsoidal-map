function map = ellipsoidal_quasiconformal_map(v,f,a,b,c,landmark,target,lambda)
% Fast ellipsoidal quasi-conformal map (FEQCM) for genus-0 closed surfaces
%
% Compute an ellipsoidal quasi-conformal parameterization of a genus-0
% closed surface with prescribed landmark constraints using the method in [1].
%
% Input:
% v: nv x 3 vertex coordinates of a genus-0 triangle mesh
% f: nf x 3 triangulations of a genus-0 triangle mesh
% a,b,c: the radii of the target ellipsoid
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

% check whether the target landmark positions lie on the ellipsoid

if max(abs(target(:,1).^2/a^2+target(:,2).^2/b^2+target(:,3).^2/c^2-1)) > 1e-5
    error('The target landmark positions do not lie on the ellipsoid!');
end

nv = length(v);

% Assume the surface is already optimally rotated
% Define north and south poles using the top and bottom points
[~,id_top] = max(v(:,3)); 
[~,id_bottom] = min(v(:,3));
northpole = v(id_top,:);
southpole = v(id_bottom,:);

% Find the faces closest to the prescribed poles
v_centroid = (v(f(:,1),:) + v(f(:,2),:) + v(f(:,3),:))/3;
[~,f_id_north] = min((v_centroid(:,1)-northpole(:,1)).^2+...
                 (v_centroid(:,2)-northpole(:,2)).^2+...
                 (v_centroid(:,3)-northpole(:,3)).^2);
[~,f_id_south] = min((v_centroid(:,1)-southpole(:,1)).^2+...
                    (v_centroid(:,2)-southpole(:,2)).^2+...
                    (v_centroid(:,3)-southpole(:,3)).^2);

%% Spherical conformal map for the input surface

% Spherical conformal map (using Choi et al., SIAM J. Imaging Sci. 2015)
S = spherical_conformal_map(v,f);

% Project the sphere onto the plane
p = stereographic(S);
z = complex(p(:,1),p(:,2));

% find the centroid of the two polar triangles
S_north = (S(f(f_id_north,1),:)+S(f(f_id_north,2),:)+S(f(f_id_north,3),:))/3;
S_south = (S(f(f_id_south,1),:)+S(f(f_id_south,2),:)+S(f(f_id_south,3),:))/3;

% project onto sphere
S_north = S_north/sqrt(S_north(:,1)^2+S_north(:,2)^2+S_north(:,3)^2);
S_south = S_south/sqrt(S_south(:,1)^2+S_south(:,2)^2+S_south(:,3)^2);

% stereographic porjections
p_north = stereographic(S_north);
p_south = stereographic(S_south);

% apply mobius transformation
z_north = complex(p_north(:,1),p_north(:,2));
z_south = complex(p_south(:,1),p_south(:,2));
% compute a conformal map that maps z(id_bottom) to 0 and z(id_top) to Inf
z = (z-z_south)./(z-z_north);

% Fixing the rotation arbitrarily
[~,id_right] = max(v(:,1));
z = z*exp(-1i*angle(z(id_right)));

%% balancing scheme

w = stereographic_ellipsoid_south_pole(stereographic_ellipsoid(z,a,b,c),a,b,c);
w = complex(w(:,1),w(:,2));

% Compute the perimeters
perimeter_north = (abs(z(f(f_id_north,1))-z(f(f_id_north,2))) + ...
    abs(z(f(f_id_north,2))-z(f(f_id_north,3))) + ...
    abs(z(f(f_id_north,3))-z(f(f_id_north,1))));
perimeter_south = (abs(w(f(f_id_south,1))-w(f(f_id_south,2))) + ...
    abs(w(f(f_id_south,2))-w(f(f_id_south,3))) + ...
    abs(w(f(f_id_south,3))-w(f(f_id_south,1))));
% rescale to get the best distribution
r_perimeter = (sqrt(perimeter_north*perimeter_south))/(perimeter_north);

z = z*r_perimeter;
p = [real(z), imag(z)];

%% Construct conformal projection from plane to ellispoid 

% apply inverse ellipsoidal stereographic projection 
Eabc = stereographic_ellipsoid(p,a,b,c);

f_punctured = f;
f_punctured(f_id_north,:) = [];

% compute Beltrami coefficient
[~,id] = sort(-abs(z));

% compute Beltrami coefficient for inverse ellipsoidal stereographic
% projection using explicit formula
x = p(:,1);
y = p(:,2);
Px = [2*a*(-x.^2+y.^2+1)./(1+x.^2+y.^2).^2, ...
    -4*b*x.*y./(1+x.^2+y.^2).^2, ...
    4*c*x./(1+x.^2+y.^2).^2];
Py = [-4*a*x.*y./(1+x.^2+y.^2).^2, ...
    2*b*(x.^2-y.^2+1)./(1+x.^2+y.^2).^2, ...
    4*c*y./(1+x.^2+y.^2).^2];
E = sum(Px.*Px,2);
F = sum(Px.*Py,2);
G = sum(Py.*Py,2);
% vertex-based BC
mu_v = (E-G+2*1i*F)./(E+G+2*sqrt(E.*G-F.^2));
% face-based BC
mu = (mu_v(f_punctured(:,1))+mu_v(f_punctured(:,2))+mu_v(f_punctured(:,3)))/3;

% Construct QC map
id_fixed = 1:min(50,nv/100);
p_lbs = linear_beltrami_solver(p,f_punctured,mu,id(id_fixed),p(id(id_fixed),:));

% Composition: conformal map from plane to Eabc
Fx = scatteredInterpolant(p_lbs(:,1),p_lbs(:,2),Eabc(:,1));
Fy = scatteredInterpolant(p_lbs(:,1),p_lbs(:,2),Eabc(:,2));
Fz = scatteredInterpolant(p_lbs(:,1),p_lbs(:,2),Eabc(:,3));

%% Landmark-aligned quasi-conformal map

% find the target landmark position on the plane
Fpx = scatteredInterpolant(p(:,1),p(:,2),p_lbs(:,1));
Fpy = scatteredInterpolant(p(:,1),p(:,2),p_lbs(:,2));
target_plane = stereographic_ellipsoid(target,a,b,c);
target_plane = [Fpx(target_plane(:,1),target_plane(:,2)), ...
    Fpy(target_plane(:,1),target_plane(:,2))];

% Apply Mobius transformation to roughly align the landmarks
zm = mobius_transformation_leastsquares(z,z(landmark,:),target_plane);
x1 = real(zm);
y1 = imag(zm);

% Landmark-constrained optimized harmonic map on the plane
L = cotangent_laplacian([x1,y1],f)/4;
M = L - sparse(landmark,landmark,lambda*ones(1,length(landmark)),nv,nv);
c = zeros(nv,1); 
d = c;
c(landmark) = -lambda*target_plane(:,1);
d(landmark) = -lambda*target_plane(:,1);

p1 = f(f_id_north,1);
p2 = f(f_id_north,2);
p3 = f(f_id_north,3);
fixed = [p1,p2,p3,find(~inpolygon(x1,y1,x1([p1 p2 p3]),y1([p1 p2 p3])))']';
[mrow,mcol,mval] = find(M(fixed,:));
M = M - sparse(fixed(mrow),mcol,mval,length(x1),length(x1)) + ...
        sparse(fixed,fixed,ones(size(fixed)),length(x1),length(x1));
c(fixed) = x1(fixed);
d(fixed) = y1(fixed);

z = M \ complex(c,d);

% Overlap correction
vv = [x1,y1];
initial_map = [real(z),imag(z)];
bdy_target = initial_map(fixed,:);
Operator = createOperator(vv',f');
map_flash = initial_map;
alpha = 1;
beta = 1;
penalty = 0.01;
sigmaIncrease = 0.5;
mu_bound = 0.99;
balancing_ratio = exp(-sqrt(lambda)/4); 

update_mu = beltrami_coefficient(vv, f, initial_map);
overlap_count = sum(abs(update_mu)>=1);
mu_diff = max(abs(update_mu - zeros(length(f),1)));
landmark_error = mean(sqrt(sum((initial_map(landmark,1:2)-target_plane).^2,2)));
fprintf('  Iteration  Mu_Difference  Overlap   Landmark error\n');
IterationNumber = 0;

fprintf('%7.0f %15f %8.0f %17f \n',[IterationNumber,mu_diff,overlap_count,landmark_error]);

while overlap_count ~= 0

    mu = update_mu;

    % Smooth BC
    penalty = penalty + sigmaIncrease;
    Smooth_Operator = speye(nv) + ...
        1/penalty*(alpha*speye(nv) - beta*L/2);
    smooth_mu = smoothing(update_mu,Smooth_Operator,Operator);
    smooth_mu(abs(smooth_mu)>=mu_bound) = ...
        smooth_mu(abs(smooth_mu)>=mu_bound)./...
        abs(smooth_mu(abs(smooth_mu)>=mu_bound))*mu_bound;

    % find BC direction for non overlap exact landmark matching
    map_flash = linear_beltrami_solver(vv,f,smooth_mu,[fixed;landmark],[bdy_target;target_plane]); 
    exact_mu = beltrami_coefficient(vv, f, map_flash);

    % combine both direction, update BC
    combined_mu = exact_mu + balancing_ratio*(smooth_mu - exact_mu);
    combined_mu(abs(combined_mu)>=mu_bound) = ...
        combined_mu(abs(combined_mu)>=mu_bound)./...
        abs(combined_mu(abs(combined_mu)>=mu_bound))*mu_bound;

    % reconstruct without any landmark constraints
    map_flash = linear_beltrami_solver(vv,f,combined_mu,fixed,bdy_target); 
    update_mu = beltrami_coefficient(vv, f, map_flash);

    % evaluate the result
    landmark_error = mean(sqrt(sum((map_flash(landmark,1:2)-target_plane).^2,2)));
    mu_diff = max(abs(update_mu - mu));
    overlap_count = sum(abs(update_mu)>=1);
    IterationNumber = IterationNumber + 1;
    fprintf('%7.0f %15f %8.0f %17f \n',[IterationNumber,mu_diff,overlap_count,landmark_error]);

    if IterationNumber > 500
        warning('Iteration exceeds 500. Automatically terminated.');
        % consider changing the parameters to improve the performance
        break;
    end
end

x = map_flash(:,1);
y = map_flash(:,2);

%% inverse mapping

map = zeros(nv,3);
map(:,1) = Fx(x,y);
map(:,2) = Fy(x,y);
map(:,3) = Fz(x,y);

end

function v = stereographic(u)
% Stereographic projection and its inverse
    if size(u, 2) == 1
      u = [real(u), imag(u)];
    end
    x = u(:,1);
    y = u(:,2);
    if size(u,2) < 3
        z = 1 + x.^2 + y.^2;
        v = [2*x ./ z, 2*y ./ z, (-1 + x.^2 + y.^2) ./ z];
        v(isnan(z)|(~isfinite(z)),1) = 0;
        v(isnan(z)|(~isfinite(z)),2) = 0;
        v(isnan(z)|(~isfinite(z)),3) = 1;
    else
        z = u(:,3);
        v = [x ./ (1-z), y ./ (1-z)];
        v(isnan(v)) = Inf;
    end
end

function P = stereographic_ellipsoid(p,a,b,c)
% Ellipsoidal stereographic projection and its inverse
    if size(p, 2) == 1
      p = [real(p), imag(p)];
    end
    u = p(:,1);
    v = p(:,2);
    if size(p,2) < 3
      z = 1 + u.^2 + v.^2;
      P = [2*a*u./z, 2*b*v./z, c*(-1+u.^2+v.^2)./z];
      
      P(isnan(z)|(~isfinite(z)),1) = 0;
      P(isnan(z)|(~isfinite(z)),2) = 0;
      P(isnan(z)|(~isfinite(z)),3) = c;
        
    else
      z = p(:,3);
      P = [(u/a)./(1-z/c), (v/b)./(1-z/c)];
      P(isnan(P)) = Inf;
    end
end

function P = stereographic_ellipsoid_south_pole(p,a,b,c)
% South-pole ellipsoidal stereographic projection and its inverse
    if size(p, 2) == 1
      p = [real(p), imag(p)];
    end
    u = p(:,1);
    v = p(:,2);
    if size(p,2) < 3
      z = 1 + u.^2 + v.^2;
      P = [2*a*u./z, 2*b*v./z, -c*(-1+u.^2+v.^2)./z];
      P(isnan(z)|(~isfinite(z)),1) = 0;
      P(isnan(z)|(~isfinite(z)),2) = 0;
      P(isnan(z)|(~isfinite(z)),3) = -c;
        
    else
      z = p(:,3);
      P = [(u/a)./(1+z/c), (v/b)./(1+z/c)];
      P(isnan(P)) = Inf;
    end
end