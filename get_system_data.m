function[system_data] = get_system_data(filenames)

%% boundary matrices for x = 1
%[system_data.BInflow] = create_sp_mat(filenames.BInflow);

%% Permutation matrix for the projector
[system_data.Perm] = create_sp_mat(filenames.Perm);
[system_data.InvPerm] = create_sp_mat(filenames.InvPerm);

%% System matrix
[system_data.Ax] = create_sp_mat(filenames.Ax);
[system_data.Ay] = create_sp_mat(filenames.Ay);
%[system_data.Sigma] = create_sp_mat(filenames.Sigma);

system_data.Ax = full(system_data.Ax);
system_data.Ay = full(system_data.Ay);
end