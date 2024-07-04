%% Single E_F run

%E_Fermi = -20*1e-3;
E_Fermi = -20*1e-3;
TBoltz_K = 10;
N_Deig_points = 50;%50; %25;
%model_ind = "paoflowTB_model";
model_ind = "kp_model";

is_single_run = 1;

addpath('core_codes_kp')  

tic
Boltzmann_Solver;
find_fkmu;
toc

observables_vs_EF = zeros(1,9);
observables_vs_EF(1) = exp_fit2_coef(1);
observables_vs_EF(2) = -exp_fit2_coef(2);
observables_vs_EF(5) = spin_conversion_efficiency(1);
%observables_vs_EF(6) = RE_coeff;
%observables_vs_EF(7) = SHE_coeff;
%observables_vs_EF(8) = total_el_current_z;
%observables_vs_EF(9) = total_spin_current_zz;

%
abs_fkmu = delta_fkmu_init_decomp_coef1;
Filename_to_save = [char(model_ind),'single_run_kp_' ,datestr(now,'mm-dd-yyyy HH-MM')];

make_plots_FIG4

%save(Filename_to_save,...
%    "E_Fermi","TBoltz_K","D_eig_Boltzmann_array1",...
%    'D_eig_Boltzmann_array','observables_vs_EF',...
%    "delta_fkmu_init_decomp_coef1","Energy_E_F",...
%    'J_polariz_matrix_el','Energy_window','spin_conversion_efficiency',...
%    'Boltzmann_ind','Energy_E_F_all','time_array',...
%    'delta_fkmu_init_decomp_coef','V_eig_Boltzmann','-v7.3');
