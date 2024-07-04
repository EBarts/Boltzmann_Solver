%% FIG4(c,d): single spectrum cut at E_F = -15 meV
% for next figures, firstly, use 
%load('single_run_kp.mat')
%load('paoflowTB_modelsingle_run_01-29-2024 23-52.mat');
Plot_fkmu_spectrum_for_paper_kp

%% FIG4(e): contours with E along z
%load('paoflowTB_modelsingle_run_01-31-2024 12-03.mat')
delta_fkmu_init_dec = zeros(size(V_eig_Boltzmann));
delta_fkmu_init_dec = squeeze(delta_fkmu_init_dec(1,:));
set_of_eigs = D_eig_Boltzmann_array > 1e-5 &...
    D_eig_Boltzmann_array < 1e1;
set_of_eigs_ind = find(set_of_eigs);
set_of_eigs_ind = set_of_eigs_ind';
for isite = set_of_eigs_ind
    delta_fkmu_init_dec(:) = delta_fkmu_init_dec(:) + ...
        delta_fkmu_init_decomp_coef(isite)*...
        V_eig_Boltzmann(:, isite);
end
Nz = 101;
c_max0 = 0;
c_max0 = Plot_func_fkmu_kp(delta_fkmu_init_dec,0,...
    Energy_window,Boltzmann_ind,Energy_E_F_all,Nz, c_max0, "total");

%% FIG4(f): fast modes, contours with E along z
delta_fkmu_init_dec = zeros(size(V_eig_Boltzmann));
delta_fkmu_init_dec = squeeze(delta_fkmu_init_dec(1,:));
D_eig_cutoff = 0.58; % for Fig3
% D_eig_cutoff = 0.05; % for Fig4

set_of_eigs = D_eig_Boltzmann_array>D_eig_cutoff &...
    D_eig_Boltzmann_array < 1e1;
set_of_eigs_ind = find(set_of_eigs);
set_of_eigs_ind = set_of_eigs_ind';
for isite = set_of_eigs_ind
    delta_fkmu_init_dec(:) = delta_fkmu_init_dec(:) + ...
        delta_fkmu_init_decomp_coef(isite)*...
        V_eig_Boltzmann(:, isite);
end
%Nz = 61;
Plot_func_fkmu_kp(delta_fkmu_init_dec,0,...
    Energy_window,Boltzmann_ind,Energy_E_F_all,Nz,c_max0,"fast");

%% FIG4(g): slow modes, contours with E along z
delta_fkmu_init_dec = zeros(size(V_eig_Boltzmann));
delta_fkmu_init_dec = squeeze(delta_fkmu_init_dec(1,:));
set_of_eigs = D_eig_Boltzmann_array>1e-5 &...
    D_eig_Boltzmann_array < D_eig_cutoff;
set_of_eigs_ind = find(set_of_eigs);
set_of_eigs_ind = set_of_eigs_ind';
for isite = set_of_eigs_ind
    delta_fkmu_init_dec(:) = delta_fkmu_init_dec(:) + ...
        delta_fkmu_init_decomp_coef(isite)*...
        V_eig_Boltzmann(:, isite);
end
%Nz = 61;
Plot_func_fkmu_kp(delta_fkmu_init_dec,0,...
    Energy_window,Boltzmann_ind,Energy_E_F_all,Nz,c_max0,"slow");