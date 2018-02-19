function output = solver(par)

% check whether we have anisoptropic scattering or not
flag_ani_scatter = isfield(par,'sigma_sm')&&isfield(par,'mom_order');

% if no anisotropic scattering then we initialize the function by a zero
if ~flag_ani_scatter % If no anisotropic part of scattering kernel given,
    par.sigma_sm = @zero; % structures for isotropic scattering, possibly
end     

if ~isfield(par,'ic'),       par.ic = @zero; end% Default: no init. cond.
if ~isfield(par,'sigma_a'),  par.sigma_a = @zero; end%Default: no absorp.
if ~isfield(par,'sigma_s0'), par.sigma_s0 = @zero; end%Defa.: no scatter.
if ~isfield(par,'source'),  par.source = @zero; end % Default: no source.
if par.num_bc ~=4
    assert(1 == 0, 'not valid num bc'); 
end        

% corresponding to every row in Ax, stores the non-zero indices
Ix = cellfun(@find,num2cell(par.system_data.Ax',1),'Un',0);         
Iy = cellfun(@find,num2cell(par.system_data.Ay',1),'Un',0);                

% in total we have four grids.x1 and x2 are the two grids in the x
% direction.
% x1 n y1
c11 = 1; 
% x1y2
c12 = [];
c21 = [];
c22 = [];
while length([c11,c12,c21,c22])<par.n_eqn                         
    c11 = unique(vertcat(c11',Ix{c21},Iy{c12})');            
    c12 = unique(vertcat(c12',Ix{c22},Iy{c11})');
    c21 = unique(vertcat(c21',Ix{c11},Iy{c22})');            
    c22 = unique(vertcat(c22',Ix{c12},Iy{c21})');           
end

assert(isempty(intersect(c11,c12)), 'incorrect ordering');
assert(isempty(intersect(c12,c22)), 'incorrect ordering');
assert(isempty(intersect(c21,c22)), 'incorrect ordering');
assert(isempty(intersect(c21,c11)), 'incorrect ordering');

% incase we do not have enough variables
if isempty(c12), c12 = []; end 
if isempty(c21), c21 = []; end 
if isempty(c22), c22 = []; end 

% 
gtx = ones(1,par.n_eqn); gty = gtx;                              % Grid types
gtx([c21,c22]) = 2; gty([c12,c22]) = 2;                     % in x and y.

% size of the grid. h1 for the x direction and h2 for the y direction
h = (par.ax([2 4])-par.ax([1 3]))./par.n;

% a crude approximation for delta_t
par.dt = min(h)/abs(eigs(par.system_data.Ax,1,'lm'))/par.CFL;

% the first and second grid in the x direction
x{1} = par.ax(1):h(1):par.ax(2);
x{2} = [par.ax(1) par.ax(1) + h(1)/2:h(1):par.ax(2)-h(1)/2 par.ax(2)];

% the first and second grid it the y direction
y{1} = par.ax(3):h(2):par.ax(4);
y{2} = [par.ax(3) par.ax(3) + h(2)/2:h(2):par.ax(4)-h(2)/2 par.ax(4)];

% develop all the four grids. X contains the x coordinates and Y contains
% the y coordinates.
[X,Y] = cellfun(@ndgrid,[x;x]',[y;y],'Un',0);


% we need to know which elemtns are coupled with which one at the
% boundaries. The 
% ID = 1
% a loop over all the boundaries
bc_coupling = cell(par.num_bc,1);
bc_coupling_g = cell(par.num_bc,1);

% % we need to do this fixing for the machine error
% for i = 1 : par.num_bc
%     bc_coupling{i} = cellfun(@(a) find(abs(a) > 1e-14) ,num2cell(par.system_data.penalty_B{i},2),'Un',0);
%     bc_coupling_g{i} = cellfun(@(a) find(abs(a)> 1e-14) ,num2cell(par.system_data.penalty{i},2),'Un',0);
% end

switch par.diff_order 
    case 2
            % differential opeartor for the X grid. There are two types of
            % x grids. DX{1} acts on something defined on 2. same for all
            % the other operators. 
           [DX{1},DX{2},PX{1},PX{2}] = sbp_staggered_2(par.n(1),h(1));
           
           % differential operators for the y grid
           [DY{1},DY{2},PY{1},PY{2}] = sbp_staggered_2(par.n(2),h(2));

           
    case 4
            % differential opeartor for the X grid. There are two types of
            % x grids. DX{1} acts on something defined on 2. same for all
            % the other operators. 
           [DX{1},DX{2},PX{1},PX{2}] = sbp_staggered_4(par.n(1),h(1));
           
           % differential operators for the y grid
           [DY{1},DY{2},PY{1},PY{2}] = sbp_staggered_4(par.n(2),h(2));

       
end

% we transpose DY for the data structure
DY = cellfun(@transpose,DY,'Un',0);

% scaling for the boundary conditions
bc_scaling = [1/PX{1}(1,1) 1/PY{1}(1,1) 1/PX{1}(1,1) 1/PY{1}(1,1)];

% initialize the solution variables
% n_sys is the total number of equations in the system
% dxU and dyU are the derivatives in the x and y direction.
% force is the forcing term 
% set the initial conditions in case they are needed else initialize with
% zero
U = cell(1,par.n_eqn); dxU = U; dyU = U; force = U; UTemp = U;

% data structure for storing the values at the boundaries.
bc_values = cell(1,par.n_eqn);
k_RK = cell(4,1);
for j = 1:par.n_eqn                          
    U{j} = X{gtx(j),gty(j)}*0+...                   
        capargs(par.ic,X{gtx(j),gty(j)},Y{gtx(j),gty(j)},j);     
   
    force{j} = X{gtx(j),gty(j)}*0;
    bc_values{j} = X{gtx(j),gty(j)}*0;
end


%% Time Loop
par.t_plot = [par.t_plot par.t_end inf]; plot_count = 1;
cputime = zeros(1,3);
t = 0; step_count = 0;

% compute the material properties at t = 0
% stores which moment order sits on which grid
sg = [gtx;gty;par.mom_order]';   
% we assume that nothing is being prescribed to the zeroth moment
sg = unique(sg(logical(par.mom_order),:),'rows')';           

s0 = cellfun(@(x,y)capargs(...            
                par.sigma_a,x,y,0),X,Y,'Un',0); 
            
sA = cellfun(@(x,y)capargs(...           
                par.sigma_a,x,y,0),X,Y,'Un',0); 
           
% sM depends upon the moment order and hence should be computed differently
sM = cell(2,2,max(par.mom_order));
for j = 1:size(sg,2)       
       sM{sg(1,j),sg(2,j),sg(3,j)} = capargs(par.sigma_sm,...
                                    X{sg(1,j),sg(2,j)},Y{sg(1,j),sg(2,j)},sg(3,j),0);
end

for j = par.source_ind     % Source: evaluate only active components.
        force{j} = capargs(par.source,X{gtx(j),gty(j)},...  % Evaluate source
                           Y{gtx(j),gty(j)},0,j); 
end

while t < par.t_end
    
    
    if t+par.dt > par.t_end
        par.dt = par.t_end-t;
    end
    
    
    % the ode sytem can be written as U_t = Op.
    % RK = 2 implementation
     tic
            t_temp = [t t + par.dt/2 t + par.dt/2 t + par.dt];
            dt_temp = [0 par.dt/2 par.dt/2 par.dt];
            weight = [1/6 2/6 2/6 1/6];
     
     UTemp = U; 
     for RK = 1 : 4
            evaluate = par.time_dep & (t_temp(RK) > 0);
            if evaluate(1)
                s0 = cellfun(@(x,y)capargs(...            
                        par.sigma_a,x,y,0),X,Y,'Un',t_temp(RK));
            end
            
            if evaluate(2)
                sA = cellfun(@(x,y)capargs(...           
                    par.sigma_a,x,y,0),X,Y,'Un',t_temp(RK)); 
            end
            
            if evaluate(3)
                for j = 1:size(sg,2)       
                    sM{sg(1,j),sg(2,j),sg(3,j)} = capargs(par.sigma_sm,...
                                                X{sg(1,j),sg(2,j)},Y{sg(1,j),sg(2,j)},sg(3,j),t_temp(RK));
                 end
            end
            
            if evaluate(4)
                for j = par.source_ind    
                         force{j} = capargs(par.source,X{gtx(j),gty(j)},... 
                                           Y{gtx(j),gty(j)},t_temp(RK),j); 
                end
            end
            
            if evaluate(5)
                % add boundary inhomogeneity here
            end
            
            for i = 1:par.n_eqn
                % we store the derivatives on the neighbouring grid. That's
                % why Ix{i}(1). Same for the y direction. 
                dxU{i} = DX{gtx(Ix{i}(1))} * UTemp{i};
                
                % we have stored the transpose of DY
                dyU{i} = UTemp{i} * DY{gty(Iy{i}(1))};
                
            end

%             % extract all the values at x = 1, last row of the matrix.
%             bc_ID = 1;
%             values = cellfun(@(a) a(end,:),U,'Un',0);
%             for j = 1 : par.n_eqn
%                 bc_values{j}(end,:) = bc_scaling(bc_ID) * sumcell(values(bc_coupling{bc_ID}{j}),...
%                                       par.system_data.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j}));
%             end
%             
%             bc_ID = 2;
%             values = cellfun(@(a) a(:,end),U,'Un',0);
%             for j = 1 : par.n_eqn
%                 bc_values{j}(:,end) = bc_scaling(bc_ID) * sumcell(values(bc_coupling{bc_ID}{j}), ...
%                                       par.system_data.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j}));
%             end
%             
%             bc_ID = 3;
%             values = cellfun(@(a) a(1,:),U,'Un',0);
%             for j = 1 : par.n_eqn
%                 bc_values{j}(1,:) = bc_scaling(bc_ID) * sumcell(values(bc_coupling{bc_ID}{j}), ...
%                                     par.system_data.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j}));
%             end
%             
%             
%             bc_ID = 4;
%             values = cellfun(@(a) a(:,1),U,'Un',0);
%             for j = 1 : par.n_eqn
%                 bc_values{j}(:,1) = bc_scaling(bc_ID) * sumcell(values(bc_coupling{bc_ID}{j}), ...
%                                     par.system_data.penalty_B{bc_ID}(j,bc_coupling{bc_ID}{j}));
%             end
            
            
            for i = 1 : par.n_eqn
                % multiplication by the system matrices
                W = -sumcell([dxU(Ix{i}),dyU(Iy{i})],...
                    [par.system_data.Ax(i,Ix{i}),par.system_data.Ay(i,Iy{i})]);
                                    
                k_RK{RK}{i} = (W - (s0{gtx(i),gty(i)}+sA{gtx(i),gty(i)}).*UTemp{i} + force{i} + ...
                                        bc_values{i});
                % if not the zeroth order moment then we also need to find
                % the contribution from the anistropic scattering. 
                if par.mom_order(i) ~= 0
                    k_RK{RK}{i} = k_RK{RK}{i} +  sM{gtx(i),gty(i),par.mom_order(i)}.*UTemp{i};
                end
                
                if RK ~= 4
                    UTemp{i} = U{i} + k_RK{RK}{i} * dt_temp(RK + 1);
                end
            end
            
            
     end
    
    for RK = 1 : 4
        for i = 1 : par.n_eqn
            U{i} = U{i} + weight(RK) * k_RK{RK}{i} * par.dt;
        end
    end
    
    step_count = step_count + 1;
    t = t + par.dt;
    cputime(1) = cputime(1) + toc;
    
    tic
    while t>=par.t_plot(plot_count)-1e-14  % If current time has exceeded
        lambda = (par.t_plot(plot_count)-t+par.dt)/par.dt; %plotting time, define
        Uplot = U{par.mom_output}; % linear interpolation in time.
        if plot_count==length(par.t_plot)-1      % Is final time reached?
            if nargout                   % If yes, save solution at final
                output = struct('x',x(gtx),'y',y(gty),'U',U);     % time.
            end
        else                           % If not, invoke plotting routine.
            xplot = x{gtx(par.mom_output)};             % Assign grids at
            yplot = y{gty(par.mom_output)};          % outputted moments.
%             if length(par.mom_output)==1        % If only a single moment
%                 xplot = xplot{:}; yplot = yplot{:}; Uplot = Uplot{par.mom_output};% is
%             end                         % plotted, remove cell structure.
            if nargout(par.output)                     % Call output rou,
                par = par.output(par,xplot,yplot,Uplot,plot_count);%tine-
            else                          % allowing for it to modify the
                par.output(par,xplot,yplot,Uplot,plot_count)% struct par.
            end
        end
        plot_count = plot_count+1;
    end
    cputime(2) = cputime(2) + toc;
end

fprintf('%0.0f time steps\n',step_count)           % Display test
cputime = reshape([cputime;cputime/sum(cputime)*1e2],1,[]);   % case info
fprintf(['CPU-times\n advection:%15.2fs%5.0f%%\n',... % and CPU times.
    'plotting:%16.2fs%5.0f%%\n'],cputime)

end


function z = capargs(fct,varargin)
% Call function fct with as many arguments as it requires (at least 1),
% and ignore further arguments.
narg = max(nargin(fct),1);
z = fct(varargin{1:narg});

end

function f = zero(varargin)
% Zero function.
f = zeros(size(varargin{1}));
end

function S = sumcell(A,w)
% Add vector of cells A, weighted with vector w.
S = A{1}*w(1); for j = 2:length(w), S = S+A{j}*w(j); end
end

