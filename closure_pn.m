function [Mx,My,moment_order] = closure_pn(n_mom)
%STARMAP_CLOSURE_PN
%   Creates P_N moment system matrices, to be used by
%   STARMAP_SOLVER, a second order staggered grid finite
%   difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 2D geometry.
%
%   Version 1.9
%   Copyright (c) 03/13/2015 Benjamin Seibold and Martin Frank
%   http://www.math.temple.edu/~seibold
%   http://www.mathcces.rwth-aachen.de/5people/frank/start
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5)
%
%   StaRMAP project website:
%   http://math.temple.edu/~seibold/research/starmap/

%   For license, see file starmap_solver.m, as published on
%   http://www.math.temple.edu/~seibold/research/starmap

%========================================================================
n_sys = (n_mom+1)*(n_mom+2)/2; % number of system components
Mx = sparse(zeros(n_sys)); My = Mx; s = size(Mx);
for m = 1:n_mom
    i = 1:m; p = m*(m-1)/2+i;
    v = D(m,-m+2*(ceil(i/2)-1));
    Mx(sub2ind(s,p,p+m)) = v;
    My(sub2ind(s,p,p+m-(-1).^i)) = -(-1).^i.*v;
    i = 1:m-1; p = m*(m-1)/2+i;
    v = F(m,-m+2+2*(ceil(i/2)-1));
    Mx(sub2ind(s,p,p+m+2)) = -v;
    My(sub2ind(s,p-(-1).^i,p+m+2)) = (-1).^i.*v;
end
m = 1:2:n_mom; i = m.*(m+1)/2;
Mx(i,:) = sqrt(2)*Mx(i,:); My(i,:) = sqrt(2)*My(i,:);
m = 2:2:n_mom; i = (m+1).*(m+2)/2;
Mx(:,i) = sqrt(2)*Mx(:,i); My(:,i) = sqrt(2)*My(:,i);
Mx = full((Mx+Mx')/2); My = full((My+My')/2);
% Order of moments
moment_order = ceil(sqrt(2*(1:n_sys)+1/4)-3/2);

%========================================================================
function y = D(l,m)
y = sqrt((l-m).*(l-m-1)/(2*l+1)/(2*l-1));

function y = F(l,m)
y = sqrt((l+m).*(l+m-1)/(2*l+1)/(2*l-1));
