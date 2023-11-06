function map = ellipsoidal_conformal_map(v,f,a,b,c)
% Fast ellipsoidal conformal map (FECM) for genus-0 closed surfaces
%
% Compute an ellipsoidal conformal parameterization of a genus-0 closed
% surface using the method in [1].
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


nv = length(v);

% Assume the surface is already optimally rotated 
% Define north and south poles using the top and bottom points
[~,id_top] = max(v(:,3)); 
[~,id_bottom] = min(v(:,3));
northpole = v(id_top,:);
southpole = v(id_bottom,:);

% Find the faces closest to the poles
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

% apply inverse ellipsoidal stereographic projection 
Eabc = stereographic_ellipsoid(p,a,b,c);

%% Construct conformal projection from plane to ellispoid 

f_punctured = f;
f_punctured(f_id_north,:) = [];

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
