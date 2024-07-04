function [abs_fkmu, observables_vs_EF,...
    Gamma_at_Ef_binned, Density_of_States_Ef] = ...
    calc_Boltzmann_observables(E_Fermi,TBoltz_K,N_Deig_points,model_ind)
    
    %% The script solves Boltzmann equations using exact diagonalization

Energy_window = 1.*1e-3;

%% Introduce wavefunctions and energy bands (Introduce_paoflowTB_model)
% Load path
% filename_paoflow =
% 'model_PAOFLOW_TB/Te_data_61x61x61_kxy=0_06_kz=0_12.mat';
filename_paoflow = 'model_PAOFLOW_TB/Te_data_41x41x41_kxy=0_06_kz=0_12.mat';
load(filename_paoflow,'ek','vk','pksp');
energy_list = ek;
clear ek

% Energy is counted near E_Fermi
Energy_E_F_all = energy_list - max(energy_list(:)) - E_Fermi;

% Select states near Ef
Boltzmann_ind = abs(Energy_E_F_all) < Energy_window;

% Total number of selected states
Nsites = max(size(find(Boltzmann_ind(:))));

% The energy array of the selected states
Energy_E_F = Energy_E_F_all(Boltzmann_ind);

wave_functions = vk.';
wave_functions = wave_functions(:,Boltzmann_ind);
clear vk

Velocity0x = squeeze(pksp(:,1));
Velocity0y = squeeze(pksp(:,2));
Velocity0z = squeeze(pksp(:,3));
Velocity0 = sqrt(Velocity0x.^2 + Velocity0y.^2 + Velocity0z.^2 );
Velocity =  Velocity0(Boltzmann_ind);
Velocity0x = Velocity0x(Boltzmann_ind);
Velocity0y = Velocity0y(Boltzmann_ind);
Velocity0z = Velocity0z(Boltzmann_ind);
clear pksp

E_Fermi0 = -0.020;

%% Solve Boltzmann equations (Boltzmann_Solver)
% energy difference cutoff (the window function), in eV
delta_E_cutoff = 1e-3;
% Kronecker Delta symbol cutoff, in eV
E_width = 0.2*1e-3;

% Density of states at E0 to calculate tao0
Density_of_States0 = Density_of_States_Ef_calc...
    (energy_list, E_Fermi0, Energy_window ,E_width,delta_E_cutoff);

% Construct the relaxation matrix
[BoltzmannKernel, Gamma_at_Ef, Density_of_States_Ef]...
    = BoltzmannKernel_calc_rel_time...
    (Energy_E_F(:), wave_functions, E_width,delta_E_cutoff, Density_of_States0);

% Make sure that it does not have large imaginary component
BoltzmannKernel = (BoltzmannKernel + BoltzmannKernel')/2;

% Diagonalize the Kernel
[V_Boltzmann, D_Boltzmann] = eig(BoltzmannKernel);

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
clear V_Boltzmann D_Boltzmann

%% Decompose excited population over eigestates \delta f_kmu (find_fkmu)
TBoltz = TBoltz_K/11605; % in eV (11605 K = 1 eV)
Fermi_step_weight = (1./(4 *cosh(Energy_E_F/2/TBoltz)));

angle_E = pi/2;%pi/2; % 0 if along x;
delta_fkmu_init_decomp_coef = zeros(1,Nsites);

%% EDELSTEIN
delta_fkmu_init = Fermi_step_weight.*...
    (cos(angle_E)*Velocity0x + sin(angle_E)*Velocity0z);

sumall = @(X) sum(X(:));
set_of_eigs = D_eig_Boltzmann_array>1e-5 &...
    D_eig_Boltzmann_array < 1e1;

%
set_of_eigs_ind = find(set_of_eigs);
set_of_eigs_ind = set_of_eigs_ind';

%set_of_eigs_ind = 1:Nsites;
for isite = set_of_eigs_ind
    delta_fkmu_init_decomp_coef(isite) = sumall(...
        conj(V_eig_Boltzmann(:, isite)).*delta_fkmu_init(:))...
    /D_eig_Boltzmann_array(isite);
end

% calculate spin decay time of this population 
% and fit it with two exponential functions (slow and fast)
Sx_matrix_el = matrix_el_calc_list(wave_functions, "Sx",model_ind);
%size(Sx_matrix_el)
%Nsites 
Lx_matrix_el = real(matrix_el_calc_list(wave_functions, "Lx",model_ind));
Sy_matrix_el = matrix_el_calc_list(wave_functions, "Sy",model_ind);
Ly_matrix_el = real(matrix_el_calc_list(wave_functions, "Ly",model_ind));
Sz_matrix_el = matrix_el_calc_list(wave_functions, "Sz",model_ind);
Lz_matrix_el = real(matrix_el_calc_list(wave_functions, "Lz",model_ind));
% Henceforward, Sz implies Mz, mu_Bohr
if floor(angle_E) == 1
    J_polariz_matrix_el = (Sz_matrix_el + Lz_matrix_el)/(3/2);
elseif floor(angle_E) == 0
    J_polariz_matrix_el = (Sx_matrix_el + Lx_matrix_el)/(3/2);
end

SL_matrix_el = zeros(6,Nsites);
SL_matrix_el(1,:) = Sx_matrix_el;
SL_matrix_el(2,:) = Sy_matrix_el;
SL_matrix_el(3,:) = Sz_matrix_el;
SL_matrix_el(4,:) = Lx_matrix_el;
SL_matrix_el(5,:) = Ly_matrix_el;
SL_matrix_el(6,:) = Lz_matrix_el;

%SL_matrix_el
delta_fkmu_spin_polarization = zeros(1,Nsites);
for isite = set_of_eigs_ind
    delta_fkmu_spin_polarization(isite) = sumall(...
        V_eig_Boltzmann(:, isite).*J_polariz_matrix_el(:));
end

%% Spin time dependence
time_array = linspace(0,1e2,1e3);
total_spin_in_time = 0*time_array;
population_in_time = 0*time_array;

for it = 1:numel(time_array)
    for isite = set_of_eigs_ind
        total_spin_in_time(it) = total_spin_in_time(it) +...
            delta_fkmu_spin_polarization(isite)*...
            delta_fkmu_init_decomp_coef(isite)*...
            exp(-time_array(it)*D_eig_Boltzmann_array(isite));
%
        population_in_time(it) = population_in_time(it) +...
            abs(delta_fkmu_init_decomp_coef(isite)*...
            exp(-time_array(it)*D_eig_Boltzmann_array(isite)));
    end
end

delta_fkmu_init_in_time = zeros(numel(total_spin_in_time),Nsites);

for it = [1, 10, 50, 100, 200]%1:10%numel(time_array)
    for isite = set_of_eigs_ind%1:Nsites
        for istates=1:Nsites
        delta_fkmu_init_in_time(it,istates) = delta_fkmu_init_in_time(it,istates) +...
            delta_fkmu_init_decomp_coef(isite)*V_eig_Boltzmann(istates, isite)...
            *exp(-time_array(it)*D_eig_Boltzmann_array(isite));
        end
    end
end

% spin_conversion_efficiency = abs(total_spin_in_time)...
%     ./population_in_time(1);

%% Calculate charge and spin currents
population_at_t0 = zeros(size(delta_fkmu_init));
el_current_x = zeros(size(delta_fkmu_init));
el_current_y = zeros(size(delta_fkmu_init));
el_current_z = zeros(size(delta_fkmu_init));
spin_current_zz = zeros(size(delta_fkmu_init));
SL_at_t0 = zeros(6, Nsites);
%for isite = set_of_eigs_ind
for istate = 1:Nsites
    population_at_t0(istate) = sumall(...
        delta_fkmu_init_decomp_coef(set_of_eigs)...
        .*V_eig_Boltzmann(istate, set_of_eigs));
    el_current_x(istate) = sumall(...
        Velocity0x(istate)...
        .*delta_fkmu_init_decomp_coef(set_of_eigs)...
        .*V_eig_Boltzmann(istate, set_of_eigs));
    el_current_y(istate) = sumall(...
        Velocity0y(istate)...
        .*delta_fkmu_init_decomp_coef(set_of_eigs)...
        .*V_eig_Boltzmann(istate, set_of_eigs));
    el_current_z(istate) = sumall(...
        Velocity0z(istate)...
        .*delta_fkmu_init_decomp_coef(set_of_eigs)...
        .*V_eig_Boltzmann(istate, set_of_eigs));
    spin_current_zz(istate) = sumall(...
        Velocity0z(istate)*J_polariz_matrix_el(istate)...
        .*delta_fkmu_init_decomp_coef(set_of_eigs)...
        .*V_eig_Boltzmann(istate, set_of_eigs));
    %SL_matrix_el()
    for iSL = 1:6
        SL_at_t0(iSL,istate) = sumall(...
            SL_matrix_el(iSL,istate)...
            .*delta_fkmu_init_decomp_coef(set_of_eigs)...
            .*V_eig_Boltzmann(istate, set_of_eigs));
    end
end
%*delta_fkmu_spin_polarization(set_of_eigs)...
%% Rashba-Edelstein coefficicent, SHE coef. and S-to-CC eff.
%population_at_t0 = sumall(abs(population_at_t0));
population_at_t0 = sum(abs(population_at_t0),"all");
%size(population_at_t0)
spin_conversion_efficiency = abs(total_spin_in_time)/population_at_t0;
 
total_el_current_z = sum(el_current_z,"all");
total_el_current_y = sum(el_current_y,"all");
total_el_current_x = sum(el_current_x,"all");
total_spin_current_zz = sumall(spin_current_zz);

%Mz_matrix_el = 2*Sz_matrix_el + Lz_matrix_el
total_Mz_at_t0 = sumall(2*SL_at_t0(3,:) + SL_at_t0(6,:));
total_Mx_at_t0 = sumall(2*SL_at_t0(1,:) + SL_at_t0(4,:));
total_My_at_t0 = sumall(2*SL_at_t0(2,:) + SL_at_t0(5,:));

%
total_SL_at_t0 = sum(SL_at_t0,2);

% convert RE_coeff to Mz, mu_Bohr in j_z = 100 A/cm^2
% forgetting 1e-8 in front!

RE_coeff_conv_factor_to_M = 1.923;


RE_coeff_zz = RE_coeff_conv_factor_to_M*...
    abs(total_Mz_at_t0/total_el_current_z);
RE_coeff_xx = RE_coeff_conv_factor_to_M*...
    abs(total_Mx_at_t0/total_el_current_x);
RE_coeff_yy = RE_coeff_conv_factor_to_M*...
    abs(total_My_at_t0/total_el_current_y);

SHE_coeff = total_spin_current_zz/total_el_current_z;

%% SHE efficiency coefficient
%return % to supress plots

% Make plots
%N_Deig_points = 50
D_eig_min = 0.;
D_eig_max = 1.0;%1.5;%1.0;%1.5;
D_eig_Boltzmann_array1 = linspace(0,D_eig_max,N_Deig_points);
D_eig_Boltzmann_hstep = D_eig_max/(N_Deig_points-1)/2;


delta_fkmu_init_decomp_coef1 = 0*D_eig_Boltzmann_array1;
for iDeig = 1:N_Deig_points
    delta_fkmu_init_decomp_coef1(iDeig) = sumall(abs(...
        delta_fkmu_init_decomp_coef(...
        D_eig_Boltzmann_array1(iDeig)-D_eig_Boltzmann_hstep<D_eig_Boltzmann_array &...
        D_eig_Boltzmann_array <= D_eig_Boltzmann_array1(iDeig)+D_eig_Boltzmann_hstep)));
end

figure
subplot(2,2,1)
hplots(1) = plot(D_eig_Boltzmann_array1',...
     abs(delta_fkmu_init_decomp_coef1)','-*');
xlim([0, 1.5])

subplot(2,2,2)
hplots(2) = plot(time_array,spin_conversion_efficiency);
ylim([0,0.5])

subplot(2,2,3)
exp_fit1 = fit(time_array',population_in_time','exp2');
plot(exp_fit1,time_array',population_in_time','--');
%val(x) = a*exp(b*x) + c*exp(d*x)
exp_fit1_coef = coeffvalues(exp_fit1);
slocal_legend = findobj('type','legend');
delete(slocal_legend);

subplot(2,2,4)
exp_fit2 = fit(time_array',total_spin_in_time','exp1');
plot(exp_fit2,time_array',total_spin_in_time','--');

%exp_fit2 = fit(time_array',total_spin_in_time','exp2');
%plot(exp_fit2,time_array',total_spin_in_time','--');
exp_fit2_coef = coeffvalues(exp_fit2);
%
%plot(Gamma_at_Ef,'*')

slocal_legend=findobj('type','legend');
delete(slocal_legend);
clear slocal_legend
% Tune plots
for ih = 1:2
    hplots(ih).LineWidth = 2;
end

myfntsize_label = 10;
myfntsize_gca = 10;

myfig = gcf;
set(gcf,'Renderer','painters')

for isubplt = 1:4
    subplot(2,2,isubplt)
%     ylabel('meV','FontSize',myfntsize_label);
    xlabel('$$t/\tau_0$$', ...
        'Interpreter', 'latex','FontSize',myfntsize_label);
    set(gca,'Fontsize',myfntsize_gca);
    set(gca,'LineWidth',1.5)
    set(gca,'box','on')
end%

subplot(2,2,1)
xlabel('$$\Gamma_i = {\tau_0}/{\tau_i}$$', ...
    'Interpreter', 'latex','FontSize',myfntsize_label);
ylabel('$$|{c_i}|,$$ at $t=0$', ...
    'Interpreter', 'latex','FontSize',myfntsize_label);

subplot(2,2,2)

ylabel('pseudospin $$\epsilon(t)$$', ...
    'Interpreter', 'latex','FontSize',myfntsize_label);

sprtinf_table_T = sprintf('T = %d K',TBoltz_K);
sprtinf_table_E = sprintf('E_F = %0.3f eV',E_Fermi);

text(0.4,0.8, sprtinf_table_T,'Units','normalized',...
    'FontSize',11)
text(0.4,0.5, sprtinf_table_E,'Units','normalized',...
        'FontSize',11)
%    'FontWeight','bold','FontSize',11)
xlim([0, 100])

subplot(2,2,3)
ylabel('$$\sum_{i}|{c_i}(t)|$$', ...
    'Interpreter', 'latex','FontSize',myfntsize_label);
xlim([0, 10])

%val(x) = a*exp(-b*x) + c*exp(-d*x)
sprtinf_exp1A = sprintf('a = %0.3f, c = %0.3f',...
    exp_fit1_coef(1),exp_fit1_coef(3));
sprtinf_exp1Gamma = sprintf('b = %0.3f, d = %0.3f',...
    -exp_fit1_coef(2),-exp_fit1_coef(4));
%fprintf
text(0.45,0.4, sprtinf_exp1A,'Units','normalized',...
    'FontSize',9)
text(0.45,0.25, sprtinf_exp1Gamma,'Units','normalized',...
    'FontSize',9)

subplot(2,2,4)
ylabel('$$<S_z>$$, a.u.', ...
    'Interpreter', 'latex','FontSize',myfntsize_label);

sprtinf_exp1A = sprintf('a = %0.3f',...
    exp_fit2_coef(1));
sprtinf_exp1Gamma = sprintf('b = %0.3f',...
    -exp_fit2_coef(2));
%fprintf
text(0.45,0.35, sprtinf_exp1A,'Units','normalized',...
    'FontSize',9)
text(0.45,0.15, sprtinf_exp1Gamma,'Units','normalized',...
    'FontSize',9)

x0 = 25;
y0 = 10;
width = 17.0;
height = 8;%7.5;%15.0;

set(gcf,'units','centimeters','position',[x0,y0,width,height])
set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
graphname = sprintf('data_temp_observables/T%d_K_and_EF%0.5f',...
    TBoltz_K, E_Fermi);
exportgraphics(myfig,[graphname,'.png']);

    %% Save the results
    observables_vs_EF = zeros(1,15);
    observables_vs_EF(1) = exp_fit2_coef(1);
    observables_vs_EF(2) = -exp_fit2_coef(2);
    observables_vs_EF(3) = spin_conversion_efficiency(1);
    observables_vs_EF(4) = RE_coeff_xx;
    observables_vs_EF(5) = RE_coeff_yy;
    observables_vs_EF(6) = RE_coeff_zz;
    observables_vs_EF(7:12) = total_SL_at_t0;
    observables_vs_EF(13) = total_el_current_x;
    observables_vs_EF(14) = total_el_current_y;
    observables_vs_EF(15) = total_el_current_z;
    %
    abs_fkmu = delta_fkmu_init_decomp_coef1;
    Gamma_at_Ef_binned = sum(Gamma_at_Ef(:));
end