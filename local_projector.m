% returns the local projector
function[proj] = local_projector(l,alpha)

listk = get_listk(l);
nEqn = l + 1;


% first we set the main diagonal
ia = 1:1:nEqn;
j = 1:1:nEqn;
v = cos(abs(listk) * alpha);

proj = sparse(ia,j,v,nEqn,nEqn);

% now we set the antidiagonal
j = nEqn:-1:1;
v = sin(listk * alpha);

proj = proj + sparse(ia,j,v,nEqn,nEqn);
end