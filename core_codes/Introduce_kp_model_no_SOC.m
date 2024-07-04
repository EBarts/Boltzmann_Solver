Nx = 81;%351;%251;%51; %201 -> 354 sec
Ny = Nx;  Nz = Nx;
kx_window = 0.05;%0.05; % in inverse Angstrem
ky_window = 0.05;%0.05; 
kz_window = 0.05;%0.07; 
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

%Akp = -32.6; Bkp = -32.6;
Akp = -32.6; Bkp = Akp;
betakp = 0; Deltakp = 0;
% 
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
    1 + 0*kz_mesh_E_F; 
% pseudospin down state
wave_functions_one_component_down =...
    0*kz_mesh_E_F; 


% The doublet of the J = \pm 3/2 states 
wave_functions = zeros(2, Nsites);
wave_functions(1, :) = wave_functions_one_component_up(:);
wave_functions(2, :) = wave_functions_one_component_down(:);
%wave_functions = wave_functions(:,Boltzmann_ind);

%% Calculate velocities
Velocity0x = 2*Akp*kx_mesh;
Velocity0y = 2*Akp*ky_mesh;
Velocity0z = 2*Bkp*kz_mesh;
Velocity0 = sqrt(Velocity0x.^2 + Velocity0y.^2 + Velocity0z.^2);
Velocity0x = Velocity0x(Boltzmann_ind);
Velocity0y = Velocity0y(Boltzmann_ind);
Velocity0z = Velocity0z(Boltzmann_ind);

Velocity = Velocity0(Boltzmann_ind);
%clear Akp Bkp Deltakp betakp

return