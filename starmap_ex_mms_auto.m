function starmap_ex_mms_auto
%STARMAP_EX_MMS_AUTO
%   Example case for STARMAP_SOLVER, a second order staggered
%   grid finite difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 2D slab geometry.
%
%   Created by the file starmap_create_mms.m
%
%   Version 1.9
%   Copyright (c) 03/13/2014 Benjamin Seibold and Martin Frank
%   http://www.math.temple.edu/~seibold
%   http://www.mathcces.rwth-aachen.de/5people/frank/start
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5)

%   For license, see file starmap_solver.m, as published on
%   http://www.math.temple.edu/~seibold/research/starmap

%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
'name','Manufactured Solution Test',... % name of example
'closure','P',... % type of closure (can be 'P' or 'SP')
'n_mom',3,... % order of moment approximation
'mom_output',1,... % order of moment approximation
'diff_order',4,... % difference order for the finite difference scheme
'CFL',2,... 
'num_bc',4,... % number of bc
't_end',0.3,... % final time (initial time is t=0)
'sigma_a',@sigma_a,... % absorption coefficient (defined below)
'sigma_s0',@sigma_s0,... % isotropic scattering coefficient (def. below)
'sigma_sm',@sigma_sm,... % aniso. scattering coefficient (defined below)
'source',@source,... % source term (defined below)
'ic',@initial,... % initial condition
'ax',[0 3 0 3],... % coordinates of computational domain%'t_plot',[0 0.0204082 0.0408163 0.0612245 0.0816327  0.102041  0.122449  0.142857  0.163265  0.183673  0.204082   0.22449  0.244898  0.265306  0.285714  0.306122  0.326531  0.346939  0.367347  0.387755  0.408163  0.428571   0.44898  0.469388  0.489796  0.510204  0.530612   0.55102  0.571429  0.591837  0.612245  0.632653  0.653061  0.673469  0.693878  0.714286  0.734694  0.755102   0.77551  0.795918  0.816327  0.836735  0.857143  0.877551  0.897959  0.918367  0.938776  0.959184  0.979592         1],... % coordinates of computational domain%'output',@output... % output routine (none here)
't_plot',[],...
'output',[]...
);

%========================================================================
% Moment System Setup and Solver Execution
%========================================================================
switch par.closure                      % Define closure matrix function.
    case  'P', closurefun = 'closure_pn';
    case 'SP', closurefun = 'starmap_closure_spn';
end                                                      % Compute moment
[par.system_data.Ax,par.system_data.Ay,par.mom_order] = feval(closurefun,par.n_mom);  % matrices.
par.n_eqn = length(par.system_data.Ax);  % size of system.

%========================================================================
% Run solver
%========================================================================
points = [16  32  64 128];

E_2 = zeros(par.n_eqn,length(points));
for k = 1:length(points)                       % Loop over various grid resolutions.
    par.n = [1 1]*points(k);                   % Numbers of grid cells.
    sol = solver(par);                         % Run solver.
    for j = 1:1                                % Loop over solution components.
        [X,Y] = ndgrid(sol(j).x,sol(j).y);     % Grid on which solution lives.
        Utrue = solution(X,Y,j,par.t_end);     % Evaluate true solution.
        D = abs(sol(j).U-Utrue);               % Difference between num. and true sol.
        int_x = dot(transpose(D),transpose(sol(j).Px * D),2);  % integral along x.
        int_xy = sum(sol(j).Py*int_x);  % integral along xy.
        E_2(j,k) = int_xy;%Sc. L2 error.
    end
end

loglog(sqrt(points),E_2(1,:));
% subplot(1,2,1);
% contourf(X,Y,sol(1).U);
% colorbar;
% 
% subplot(1,2,2);
% contourf(X,Y,Utrue);
% colorbar;






%========================================================================
% Problem Specific Functions
%========================================================================
function f = initial(x,y,k)
switch k
case 1, f = feval(@(x,y,t)exp((t-x+3.0./2.0).^2.*-2.0e1-(t-y+3.0./2.0).^2.*2.0e1),x,y,0);
case 2, f = feval(@(x,y,t)0.0,x,y,0);
case 3, f = feval(@(x,y,t)0.0,x,y,0);
case 4, f = feval(@(x,y,t)0.0,x,y,0);
case 5, f = feval(@(x,y,t)0.0,x,y,0);
case 6, f = feval(@(x,y,t)0.0,x,y,0);
case 7, f = feval(@(x,y,t)0.0,x,y,0);
case 8, f = feval(@(x,y,t)0.0,x,y,0);
case 9, f = feval(@(x,y,t)0.0,x,y,0);
case 10, f = feval(@(x,y,t)0.0,x,y,0);
end

function f = source(x,y,t,k)
switch k
case 1, f = feval(@(x,y,t)exp((t-x+3.0./2.0).^2.*-2.0e1-(t-y+3.0./2.0).^2.*2.0e1).*(t.*-8.0e1+x.*4.0e1+y.*4.0e1+t.*cos(y.*pi.*2.0)-1.2e2),x,y,t);
case 2, f = feval(@(x,y,t)sqrt(3.0).*exp((t-x+3.0./2.0).^2.*-2.0e1-(t-y+3.0./2.0).^2.*2.0e1).*(t.*4.0e1-x.*4.0e1+6.0e1).*(1.0./3.0),x,y,t);
case 3, f = feval(@(x,y,t)sqrt(3.0).*exp((t-x+3.0./2.0).^2.*-2.0e1-(t-y+3.0./2.0).^2.*2.0e1).*(t.*4.0e1-y.*4.0e1+6.0e1).*(1.0./3.0),x,y,t);
case 4, f = feval(@(x,y,t)0.0,x,y,t);
case 5, f = feval(@(x,y,t)0.0,x,y,t);
case 6, f = feval(@(x,y,t)0.0,x,y,t);
case 7, f = feval(@(x,y,t)0.0,x,y,t);
case 8, f = feval(@(x,y,t)0.0,x,y,t);
case 9, f = feval(@(x,y,t)0.0,x,y,t);
case 10, f = feval(@(x,y,t)0.0,x,y,t);
end

function f = solution(x,y,k,t)
switch k
case 1, f = feval(@(x,y,t)exp((t-x+3.0./2.0).^2.*-2.0e1-(t-y+3.0./2.0).^2.*2.0e1),x,y,t);
case 2, f = feval(@(x,y,t)0.0,x,y,t);
case 3, f = feval(@(x,y,t)0.0,x,y,t);
case 4, f = feval(@(x,y,t)0.0,x,y,t);
case 5, f = feval(@(x,y,t)0.0,x,y,t);
case 6, f = feval(@(x,y,t)0.0,x,y,t);
case 7, f = feval(@(x,y,t)0.0,x,y,t);
case 8, f = feval(@(x,y,t)0.0,x,y,t);
case 9, f = feval(@(x,y,t)0.0,x,y,t);
case 10, f = feval(@(x,y,t)0.0,x,y,t);
end

function f = sigma_a(x,y,t)
f = feval(@(x,y,t)t.*cos(y.*pi.*2.0),x,y,t);

function f = sigma_s0(x,y,t)
f = feval(@(x,y,t)1.0,x,y,t);

function f = sigma_sm(x,y,m,t)
switch m
case 0, f = feval(@(x,y,t)1.0,x,y,t);
case 1, f = feval(@(x,y,t)9.0./1.0e1,x,y,t);
case 2, f = feval(@(x,y,t)8.1e1./1.0e2,x,y,t);
case 3, f = feval(@(x,y,t)7.29e2./1.0e3,x,y,t);
end

function output(par,x,y,U,step)

imagesc(x,y,U'), axis xy equal tight;
title(sprintf('t = %0.2f',par.t_plot(step)));
colorbar;
xlabel('x'), ylabel('y')

drawnow