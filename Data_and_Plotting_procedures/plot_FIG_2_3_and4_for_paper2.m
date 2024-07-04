%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for Ploting figures 2 and 3 for paper
% run first to get the file (or use preloaded):
% find_EF_and_T_dependence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIG2: efficiency and Mz
% load file first

% zz response
% 81x81x81% T=10K
%load('paoflowTB_model_run_01-26-2024 17-34.mat') 
%eps_efficiency = observables_vs_EF(5,:);
%M_Edelstein = observables_vs_EF(6,:);

load('paoflowTB_model_run_02-12-2024 06-48.mat')
%load('paoflowTB_model_run_02-22-2024 04-28.mat')
eps_efficiency_zz = observables_vs_EF(3,:);
M_Edelstein_zz = observables_vs_EF(6,:);
Spin_momentum_zz = observables_vs_EF(9,:);
Orbital_momentum_zz = observables_vs_EF(12,:);
%observables_vs_EF(7:12) = total_SL_at_t0;
load('paoflowTB_model_run_02-13-2024 01-36.mat')
eps_efficiency_xx = observables_vs_EF(3,:);
M_Edelstein_xx = observables_vs_EF(4,:);
%Spin_momentum = observables_vs_EF(9,:);
%Orbital_momentum = observables_vs_EF(12,:);
Plot_observable_for_paper3(E_Fermi_array,...
    eps_efficiency_zz, eps_efficiency_xx,...
    M_Edelstein_zz, M_Edelstein_xx,...
    Spin_momentum_zz, Orbital_momentum_zz);

%% FIG3(a): spectral decomposition with E field along z
load('paoflowTB_model_run_01-09-2024 17-15.mat')% 81x81x81% T=10K
%load('paoflowTB_model_run_01-26-2024 17-34.mat')% 101x101x101% T=10K
Sz_relax_time = squeeze(observables_vs_EF(2,:));
Plot_Ef_and_T_dep_for_paper2(fkmu_init_vs_EF,...
    E_Fermi_array,Sz_relax_time,N_Deig_points,"EFdep",TBoltz_K)

%% FIG3(b): spectral decomposition with E field along x
load('paoflowTB_model_run_01-09-2024 18-53.mat') %81x81x81% T=10K
Plot_Ef_and_T_dep_for_paper2(fkmu_init_vs_EF,...
    E_Fermi_array,0,N_Deig_points,"EFdep",TBoltz_K)

%% FIG3(c,d): single spectrum cut at E_F = -15 meV
% for next figures, firstly, use 
%load('single_run_paoflow.mat')
load('paoflowTB_modelsingle_run_01-26-2024 19-42.mat')
Plot_fkmu_spectrum_for_paper2
%Plot_fkmu_spectrum_for_paper1pst % for FIG 4

%% FIG3(e): contours with E along z
% load('paoflowTB_modelsingle_run_01-31-2024 12-03.mat')
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
Nz = 61;
c_max0 = 0;
c_max0 = Plot_func_fkmu_PAOFLOW4(delta_fkmu_init_dec,0,...
    Energy_window,Boltzmann_ind,Energy_E_F_all,Nz, c_max0, "total");

%% FIG3(f): fast modes, contours with E along z
delta_fkmu_init_dec = zeros(size(V_eig_Boltzmann));
delta_fkmu_init_dec = squeeze(delta_fkmu_init_dec(1,:));
%D_eig_cutoff = 0.33; % for Fig3
D_eig_cutoff = 0.05; % for Fig4

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
Plot_func_fkmu_PAOFLOW4(delta_fkmu_init_dec,0,...
    Energy_window,Boltzmann_ind,Energy_E_F_all,Nz,c_max0,"fast");

%% FIG3(g): slow modes, contours with E along z
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
Plot_func_fkmu_PAOFLOW4(delta_fkmu_init_dec,0,...
    Energy_window,Boltzmann_ind,Energy_E_F_all,Nz,c_max0,"slow");