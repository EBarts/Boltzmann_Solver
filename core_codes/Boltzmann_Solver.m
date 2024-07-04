%% The script solves Boltzmann equations using exact diagonalization

Energy_window = 1.*1e-3;%1e-3;
%Energy_window = 2*1e-3;%1e-3;

if is_single_run == 1
     Energy_window = 1.*1e-3;
end

%
if model_ind == "paoflowTB_model"
    Introduce_paoflowTB_model
    E_Fermi0 = -0.020;
elseif model_ind == "kp_model"
    Introduce_kp_model
    E_Fermi0 = -0.005;%-0.002;
elseif model_ind == "toyTB_model_3D"
    Introduce_toyTB_model 
    E_Fermi0 = -0.005;%-0.002;
elseif model_ind == "kp_model_free_2DEG"  
    Introduce_kp_model_no_SOC; 
    E_Fermi0 = -0.002;
end


%% Solve Boltzmann equations
% energy difference cutoff (the window function), in eV
% delta_E_cutoff = 0.2*1e-3;
% % Kronecker Delta symbol cutoff, in eV
% E_width = 0.1*1e-3;

delta_E_cutoff = 1*1e-3;%0.2*1e-3;
% Kronecker Delta symbol cutoff, in eV
E_width = 0.2*1e-3;


% [~, ~, Density_of_States0]...
%     = BoltzmannKernel_calc_rel_time...
%     (Energy_E_F(:) + E_Fermi - E_Fermi0, wave_functions, E_width,delta_E_cutoff, 1);

Density_of_States0 = Density_of_States_Ef_calc...
    (energy_list, E_Fermi0, Energy_window ,E_width,delta_E_cutoff);

[BoltzmannKernel, Gamma_at_Ef, Density_of_States_Ef]...
    = BoltzmannKernel_calc_rel_time...
    (Energy_E_F(:), wave_functions, E_width,delta_E_cutoff, Density_of_States0);


BoltzmannKernel = (BoltzmannKernel + BoltzmannKernel')/2;

[V_Boltzmann, D_Boltzmann] = eig(BoltzmannKernel);

%D_Boltzmann(1:19)

% store the eigenvalues and normalize the eigenvectors
D_eig_Boltzmann_array = zeros(Nsites,1); % eigen energies
V_eig_Boltzmann = zeros(Nsites,Nsites);
V_eig_Boltzmann_square = 0*D_eig_Boltzmann_array;

for isite = 1:Nsites
    D_eig_Boltzmann_array(isite) = D_Boltzmann(isite,isite);
    V_eig_Boltzmann_square(isite) = ...
        sqrt(V_Boltzmann(:,isite)'*V_Boltzmann(:,isite));
    V_eig_Boltzmann(:, isite) = V_Boltzmann(:,isite)...
        /sqrt(V_Boltzmann(:,isite)'*V_Boltzmann(:,isite));
end

%D_eig_Boltzmann_array(1:19)

clear V_Boltzmann D_Boltzmann
%toc
return
