% filename_paoflow = 'model_PAOFLOW_TB/output_Te_data2.mat';
% filename_paoflow = 'model_PAOFLOW_TB/output_Te_data.mat';
%filename_paoflow = 'model_PAOFLOW_TB/Te_data_401x1x401_kxy=0_06_kz=0_12.mat';
filename_paoflow = 'model_PAOFLOW_TB/Te_data_41x41x41_kxy=0_06_kz=0_12.mat';

%filename_paoflow = 'model_PAOFLOW_TB/Te_data_61x61x61_kxy=0_06_kz=0_12.mat';

%filename_paoflow = 'model_PAOFLOW_TB/Te_data_81x81x81_kxy=0_06_kz=0_12.mat';

%filename_paoflow = 'model_PAOFLOW_TB/Te_data_101x101x101_kxy=0_06_kz=0_12.mat';

load(filename_paoflow,'ek','vk','pksp');
energy_list = ek;
clear ek

Energy_E_F_all = energy_list - max(energy_list(:)) - E_Fermi;
% Select states near Ef
Boltzmann_ind = abs(Energy_E_F_all) < Energy_window;
% total number of selected states
Nsites = max(size(find(Boltzmann_ind(:))));
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