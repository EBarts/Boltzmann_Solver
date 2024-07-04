%% Decompose excited population over eigestates \delta f_kmu
%E_Fermi = -0.02;
%TBoltz_K = 2.0 %2.5 % in Kelvin
TBoltz = TBoltz_K/11605; % in eV (11605 K = 1 eV)
%Fermi_step_weight = (1./(4 *cosh(Energy_E_F/2/TBoltz)));
Fermi_step_weight = (1./cosh(Energy_E_F/2/TBoltz)).^2/4;

%Velocity0x = zeros(Nz, Nx, dNbands);
%Velocity0z = zeros(Nz, Nx, dNbands);
angle_E = pi/2;%pi/2; % 0 if along x;
%delta_fkmu_init = zeros(Nz, Nx, dNbands);
delta_fkmu_init_decomp_coef = zeros(1,Nsites);

%% EDELSTEIN
delta_fkmu_init = Fermi_step_weight.*...
    (cos(angle_E)*Velocity0x + sin(angle_E)*Velocity0z);

%% Optically PST:
%delta_fkmu_init = Fermi_step_weight.*abs(Velocity0z)...
%     .*sign(Energy_E_F);


% delta_fkmu_init = delta_fkmu_init/...
%     max(abs(delta_fkmu_init(:)));

sumall = @(X) sum(X(:));
set_of_eigs = D_eig_Boltzmann_array>1e-5 &...
    D_eig_Boltzmann_array < 1e1;

% set_of_eigs = D_eig_Boltzmann_array>1e-3 &...
%     D_eig_Boltzmann_array < 1e1;

% set_of_eigs = D_eig_Boltzmann_array>1e-10 &...
%     D_eig_Boltzmann_array < 1e5;

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

if model_ind == "kp_model"||...
        model_ind == "kp_model_free_2DEG"
% for Iz pseudospin expectation value
%Sz_matrix_el = ...
J_polariz_matrix_el = ...
    wave_functions_one_component_up.*wave_functions_one_component_up-...
    wave_functions_one_component_down.*wave_functions_one_component_down;
Mz_matrix_el = J_polariz_matrix_el;
% elseif model_ind == "toyTB_model_3D" 
% Jz_polariz_matrix_el = ...
%         (matrix_el_calc(V_eig_array, Ny, Nx, Nz, Nbands, "Sz")...
%         + real(matrix_el_calc(V_eig_array, Ny, Nx, Nz, Nbands, "Lz")))/(3/2);
%     Mz_matrix_el = ...
%         2*matrix_el_calc(V_eig_array, Ny, Nx, Nz, Nbands, "Sz")...
%         + real(matrix_el_calc(V_eig_array, Ny, Nx, Nz, Nbands, "Lz"));
%     Jz_polariz_matrix_el = Jz_polariz_matrix_el(Boltzmann_ind);
%     Mz_matrix_el = Mz_matrix_el(Boltzmann_ind);
%     Sz_matrix_el = matrix_el_calc(V_eig_array, Ny, Nx, Nz, Nbands, "Sz");
% Sz_matrix_el = Sz_matrix_el(Boltzmann_ind);

% Sz_matrix_el = matrix_el_calc_list...
%     (V_eig_array(:,Boltzmann_ind), Nsites, Nbands, "Sz");
elseif model_ind == "paoflowTB_model" ||...
        model_ind == "toyTB_model_3D" 
    Sx_matrix_el = matrix_el_calc_list(wave_functions, "Sx",model_ind);
    %size(Sx_matrix_el)
    %Nsites 
    Lx_matrix_el = real(matrix_el_calc_list(wave_functions, "Lx",model_ind));
    Sy_matrix_el = matrix_el_calc_list(wave_functions, "Sy",model_ind);
    Ly_matrix_el = real(matrix_el_calc_list(wave_functions, "Ly",model_ind));
    Sz_matrix_el = matrix_el_calc_list(wave_functions, "Sz",model_ind);
    Lz_matrix_el = real(matrix_el_calc_list(wave_functions, "Lz",model_ind));
    %size(wave_functions)
    %wave_functions(1:5,1:5)
    %save('Matrices_results_MATLAB.mat','Lz_matrix_el','wave_functions')
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
    %Mz_matrix_el = 2*Sz_matrix_el + Lz_matrix_el;
end
%Sz_matrix_el(1:5)
%Lz_matrix_el(1:5)
%SL_matrix_el
delta_fkmu_spin_polarization = zeros(1,Nsites);
for isite = set_of_eigs_ind
    delta_fkmu_spin_polarization(isite) = sumall(...
        V_eig_Boltzmann(:, isite).*J_polariz_matrix_el(:));
end

%% Spin time dependence
% time_array = linspace(0,5,1e4);
% time_array = linspace(10,1e3,1e4); %exp(-10) = 4.5400e-05
%time_array = linspace(0,1e2,1e4);
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
if model_ind == "paoflowTB_model"
    RE_coeff_conv_factor_to_M = 1.923;
elseif model_ind == "kp_model" || ...
        model_ind == "toyTB_model_3D"||...
        model_ind == "kp_model_free_2DEG"
    RE_coeff_conv_factor_to_M = 1.923*0.529;
end

RE_coeff_zz = RE_coeff_conv_factor_to_M*...
    abs(total_Mz_at_t0/total_el_current_z);
RE_coeff_xx = RE_coeff_conv_factor_to_M*...
    abs(total_Mx_at_t0/total_el_current_x);
RE_coeff_yy = RE_coeff_conv_factor_to_M*...
    abs(total_My_at_t0/total_el_current_y);

%real([total_Mx_at_t0, total_My_at_t0, total_Mz_at_t0])
%[total_el_current_x, total_el_current_y, total_el_current_z]
%real([RE_coeff_xx, RE_coeff_yy, RE_coeff_zz])

SHE_coeff = total_spin_current_zz/total_el_current_z;

%% SHE efficiency coefficient
%return % to supress plots

% Make plots
%N_Deig_points = 50
D_eig_min = 0.;
D_eig_max = 1.0;%1.5;%1.0;%1.5;
D_eig_Boltzmann_array1 = linspace(0,D_eig_max,N_Deig_points);
D_eig_Boltzmann_hstep = D_eig_max/(N_Deig_points-1)/2;
    %(D_eig_Boltzmann_array1(2) - D_eig_Boltzmann_array1(1))/2

%Gamma_at_Ef1 = linspace(0,D_eig_max,N_Deig_points);

delta_fkmu_init_decomp_coef1 = 0*D_eig_Boltzmann_array1;
for iDeig = 1:N_Deig_points
    delta_fkmu_init_decomp_coef1(iDeig) = sumall(abs(...
        delta_fkmu_init_decomp_coef(...
        D_eig_Boltzmann_array1(iDeig)-D_eig_Boltzmann_hstep<D_eig_Boltzmann_array &...
        D_eig_Boltzmann_array <= D_eig_Boltzmann_array1(iDeig)+D_eig_Boltzmann_hstep)));
end

%to suppress plots
%return

% delta_fkmu_init_decomp_coef1(1) = sumall(...
%         delta_fkmu_init_decomp_coef(...
%         D_eig_Boltzmann_array <= D_eig_Boltzmann_array1(1)));
% for iDeig = 2:N_Deig_points
%     delta_fkmu_init_decomp_coef1(iDeig) = sumall(...
%         delta_fkmu_init_decomp_coef(...
%         D_eig_Boltzmann_array1(iDeig-1)<D_eig_Boltzmann_array &...
%         D_eig_Boltzmann_array <= D_eig_Boltzmann_array1(iDeig)));
% end
% shift by the half of the bin
% D_eig_Boltzmann_array1 = ...
%     D_eig_Boltzmann_array1 - D_eig_max/N_Deig_points/2;
figure
subplot(2,2,1)
%exp_fit1 = fit(time_array',population_in_time','exp2');
%plot(exp_fit1,time_array',population_in_time','--');
% gauss_fit = fit(D_eig_Boltzmann_array1',...
%     abs(delta_fkmu_init_decomp_coef1)','gauss2');
hplots(1) = plot(D_eig_Boltzmann_array1',...
     abs(delta_fkmu_init_decomp_coef1)','-*');
% hold on
% % ft_gauss = fittype('a*exp(-((x-b)/c)^2)+d');
% [curve_fit,~] = fit(D_eig_Boltzmann_array1'...
%     ,abs(delta_fkmu_init_decomp_coef1)','gauss2');
% plot(curve_fit,'m')
% hold off
% slocal_legend=findobj('type','legend');
% delete(slocal_legend);
%xlim([0, 1.])
xlim([0, 1.5])

subplot(2,2,2)
hplots(2) = plot(time_array,spin_conversion_efficiency);
ylim([0,0.5])

subplot(2,2,3)
%plot(time_array,population_in_time);
%xlim([10, time_array(end)])
% exp_fit1 = fit(time_array',total_spin_in_time','exp1');
% plot(exp_fit1,time_array',population_in_time','--');
exp_fit1 = fit(time_array',population_in_time','exp2');
plot(exp_fit1,time_array',population_in_time','--');
%size(plot(exp_fit1,time_array',population_in_time','--'))

%val(x) = a*exp(b*x) + c*exp(d*x)
exp_fit1_coef = coeffvalues(exp_fit1);
slocal_legend = findobj('type','legend');
delete(slocal_legend);

subplot(2,2,4)
%plot(time_array,total_spin_in_time);
%xlim([10, time_array(end)])
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
%text(-0.2,1., 'b','Units','normalized',...
% text(0.1,0.95, 'b','Units','normalized',...
%     'FontWeight','bold','FontSize',11)
%ylim([1.05*min(scaleNNN(NNN_J/100)), ...
%    1.05*max(scaleNNN(NNN_K))]);
%

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
%val(x) = a*exp(-b*x) + c*exp(-d*x)
% sprtinf_exp1A = sprintf('a = %0.3f, c = %0.3f',...
%     exp_fit2_coef(1),exp_fit2_coef(3));
% sprtinf_exp1Gamma = sprintf('b = %0.3f, d = %0.3f',...
%     -exp_fit2_coef(2),-exp_fit2_coef(4));
sprtinf_exp1A = sprintf('a = %0.3f',...
    exp_fit2_coef(1));
sprtinf_exp1Gamma = sprintf('b = %0.3f',...
    -exp_fit2_coef(2));
%fprintf
text(0.45,0.35, sprtinf_exp1A,'Units','normalized',...
    'FontSize',9)
text(0.45,0.15, sprtinf_exp1Gamma,'Units','normalized',...
    'FontSize',9)

%subplot(2,2,1)
%text(-0.2,1., 'a','Units','normalized',...
% text(0.1,0.95, 'a','Units','normalized',...
%     'FontWeight','bold','FontSize',11)
%ylim([1.05*min(Ac_array),1.05*max(NN_K)]);

%subplot(2,2,2)
%text(-0.2,1., 'b','Units','normalized',...
% text(0.1,0.95, 'b','Units','normalized',...
%     'FontWeight','bold','FontSize',11)
%ylim([1.05*min(scaleNNN(NNN_J/100)), ...
%    1.05*max(scaleNNN(NNN_K))]);
%


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
%exportgraphics(myfig,[graphname,'.pdf']);
%exportgraphics(myfig,'Spin_relaxation.pdf')

