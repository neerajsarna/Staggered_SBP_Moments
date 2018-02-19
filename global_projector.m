function[projector] = global_projector(MomOrder,alpha,Perm,InvPerm)

l_values = 0:1:MomOrder;

% coefficients for a particular l
num_l = l_values + ones(1,length(l_values));

% total number of equations
nEqn = ((MomOrder + 1 ) * (MomOrder + 2))/2;

% position of each block of l
pos_l = zeros(MomOrder+1,1);

projector = zeros(nEqn,nEqn);

pos_l(1) = 1;

for i = 2:MomOrder+1
    for j = 1:i-1
        % add all the previous lengths 
        pos_l(i) = pos_l(i) + num_l(j);
    end
        pos_l(i) = pos_l(i) + 1;
end

for l = 0 : 1 : MomOrder
    global_id = pos_l(l+1):num_l(l+1) + pos_l(l+1)-1;
    projector(global_id,global_id) = local_projector(l,alpha);
end

% develop projector for u_{perm}. See pdf
projector = Perm*projector*InvPerm;
projector = sparse(projector);
end