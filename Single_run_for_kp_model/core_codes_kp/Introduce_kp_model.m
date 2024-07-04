Nx = 101;%251;%51; %201 -> 354 sec
Ny = Nx;  Nz = Nx;
kx_window = 0.04;%0.05; % in inverse Angstrem
ky_window = 0.04;%0.05; 
kz_window = 0.06;%0.07; 
%
slice_to_2D_ind = 0;
if slice_to_2D_ind == 1
    ky_window = 0;
    Ny = 1;
end

[kx_mesh,ky_mesh, kz_mesh] = ndgrid(...
    linspace(-kx_window,+kx_window,Nx),...
    linspace(-ky_window,+ky_window,Ny),...
    linspace(-kz_window,+kz_window,Nz));

Akp = -32.6; Bkp = -36.4;
betakp = 2.47; 
Deltakp = 63*1e-3;
%Deltakp = 0.6*63*1e-3;
%Deltakp = 0.2*63*1e-3;

% the upper branch in eV is
energy_list = ...
    Akp*(kx_mesh.^2 + ky_mesh.^2) + Bkp*kz_mesh.^2 + ...
    sqrt(betakp^2*kz_mesh.^2 + Deltakp^2);

Energy_E_F_all = energy_list - max(energy_list(:)) - E_Fermi;

% Select states near Ef
Boltzmann_ind = abs(Energy_E_F_all)< Energy_window;
% total number of selected states
Nsites = max(size(find(Boltzmann_ind(:))));
Energy_E_F = Energy_E_F_all(Boltzmann_ind);

kx_mesh_E_F = kx_mesh(Boltzmann_ind);
ky_mesh_E_F = ky_mesh(Boltzmann_ind);
kz_mesh_E_F = kz_mesh(Boltzmann_ind);

% pseudospin up state
wave_functions_one_component_up =...
    sqrt(1 + betakp*kz_mesh_E_F./...
    sqrt(betakp^2*kz_mesh_E_F.^2 + Deltakp^2))/sqrt(2); 
% pseudospin down state
wave_functions_one_component_down =...
    sqrt(1 - betakp*kz_mesh_E_F./...
    sqrt(betakp^2*kz_mesh_E_F.^2 + Deltakp^2))/sqrt(2); 

% The doublet of the J = \pm 3/2 states 
wave_functions = zeros(2, Nsites);
wave_functions(1, :) = wave_functions_one_component_up(:);
wave_functions(2, :) = wave_functions_one_component_down(:);

%% Calculate velocities
% for the upper band ('plus')
Velocity0x = 2*Akp*kx_mesh;
Velocity0y = 2*Akp*ky_mesh;
Velocity0z = (2*Bkp + betakp^2./sqrt(betakp^2*kz_mesh.^2 + Deltakp^2)).*kz_mesh;
Velocity0 = sqrt(Velocity0x.^2 + Velocity0y.^2 + Velocity0z.^2);
Velocity0x = Velocity0x(Boltzmann_ind);
Velocity0y = Velocity0y(Boltzmann_ind);
Velocity0z = Velocity0z(Boltzmann_ind);

Velocity = Velocity0(Boltzmann_ind);
clear Akp Bkp Deltakp betakp

return

%%
E_Fermi = -0.004; Energy_window = 0.001; Introduce_kp_model
figure
velocity_convers_to_cm_per_s = 2.418; % keep in mind that 10 ^6 is neglected
surf(velocity_convers_to_cm_per_s*squeeze(Velocity0))
shading flat
colorbar
figure
surf(velocity_convers_to_cm_per_s*squeeze(energy_list - max(energy_list(:))))

