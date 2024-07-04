%%
%filename_toyTB = 'model_toyTB_3D/bands_toyTB_3Dmodel.mat';
% filename_toyTB = ['model_toyTB_3D/',filename];
filename_toyTB = 'model_toyTB_3D/bands_toyTB_3Dmodel_81_81_81.mat';
%    'model_toyTB_3D/bands_toyTB_3Dmodel_71_71_71.mat';
%    'model_toyTB_3D/bands_toyTB_3Dmodel_81_81_81.mat';
%    'model_toyTB_3D/bands_toyTB_3Dmodel_71_71_71.mat';
%    'model_toyTB_3D/bands_toyTB_3Dmodel_61_61_61';
load(filename_toyTB,'D_eig_array','V_eig_array',...
    'DeltaLambda','a_in_triangle','c_lat','ab_lat',...
    'kx_lattice','ky_lattice','kz_lattice');
c_lat_in_Ang = c_lat*a_in_triangle;
ab_lat_in_Ang = ab_lat*a_in_triangle;
iband = 12;
% or iband = 11:12 for two bands
%%
[Ny, Nz, Nx, Nbands] = size(D_eig_array);
%D_eig_array = zeros(Ny, Nz, Nx, Nbands); % eigen energies
%V_eig_array = zeros(Ny, Nz, Nx, Nbands,Nbands); % eigen states
energy_list = zeros(Ny,Nz, Nx);
energy_list(:,:,:) = DeltaLambda*1e-3*D_eig_array(:,:,:, iband);
% energy_list = ek;
% clear ek

Energy_E_F_all = energy_list - max(energy_list(:)) - E_Fermi;
% Select states near Ef
Boltzmann_ind = abs(Energy_E_F_all)< Energy_window;
% total number of selected states
Nsites = max(size(find(Boltzmann_ind(:))));
Energy_E_F = Energy_E_F_all(Boltzmann_ind);

wave_functions_one_component = zeros(Ny,Nz, Nx);
%Nsites = Ny*Nz*Nx;
wave_functions = zeros(Nbands, Ny*Nz*Nx);
%
for icomp = 1:Nbands
    wave_functions_one_component(:,:,:) = ...
        V_eig_array(:,:,:,icomp,iband);
    wave_functions(icomp, :) = wave_functions_one_component(:);
end

wave_functions = wave_functions(:,Boltzmann_ind);


%%
%if Ny == 1
% Calculate velocities
%Velocity0 = NaN + zeros(Nz, 6);
Velocity0a = zeros(Ny,Nz, Nx);
Velocity0z = zeros(Ny,Nz, Nx);
Velocity0b = zeros(Ny,Nz, Nx);

%Velocity0 = zeros(Ny, Nz, Nx, dNbands);
kz_lattice_unit = kz_lattice(10) - kz_lattice(9);
kx_lattice_unit = kx_lattice(10) - kx_lattice(9);
if Ny ~= 1
    ky_lattice_unit = ky_lattice(10) - ky_lattice(9);
end
    
% z velocity component
for ikz = 2:Nz-1
    Velocity0z(:,ikz,:) = ...
        (energy_list(:,ikz+1,:) - energy_list(:,ikz-1,:) )/(2*kz_lattice_unit);
end
% a velocity component
for ikx = 2:Nx-1
    Velocity0a(:,:,ikx) = ...
        (energy_list(:,:,ikx+1) - energy_list(:,:,ikx-1) )/(2*kx_lattice_unit);
end

if Ny ~= 1
    for iky = 2:Ny-1
        Velocity0b(iky,:,:) = ...
            (energy_list(iky+1,:,:) - energy_list(iky - 1,:,:) )/(2*ky_lattice_unit);
    end
end

% edges:
% z velocity component
ikz = 1;
Velocity0z(:,ikz,:) = ...
    (energy_list(:,ikz+1,:) - energy_list(:,ikz,:) )/kz_lattice_unit;
ikz = Nz;
Velocity0z(:,ikz,:) = ...
    (energy_list(:,ikz,:) - energy_list(:,ikz-1,:) )/kz_lattice_unit;
% a velocity component
ikx = 1;
Velocity0a(:,:,ikx) = ...
    (energy_list(:,:,ikx+1) - energy_list(:,:,ikx) )/kx_lattice_unit;
ikx = Nx;
Velocity0a(:,:,ikx) = ...
    (energy_list(:,:,ikx) - energy_list(:,:,ikx-1) )/kx_lattice_unit;
% b velocity component

if Ny~=1
iky = 1;
Velocity0b(iky,:,:) = ...
    (energy_list(iky+1,:,:) - energy_list(iky,:,:) )/ky_lattice_unit;
iky = Ny;
Velocity0b(iky,:,:) = ...
    (energy_list(iky,:,:) - energy_list(iky-1,:,:) )/ky_lattice_unit;
end


Velocity0a = Velocity0a*ab_lat_in_Ang*sqrt(3)/(4*pi);
Velocity0b = Velocity0b*ab_lat_in_Ang*sqrt(3)/(4*pi);

Velocity0z = Velocity0z*c_lat_in_Ang/(2*pi);
Velocity0x = Velocity0a - Velocity0b;
Velocity0y = (Velocity0a + Velocity0b)/sqrt(3);

Velocity0 = sqrt(Velocity0x.^2 + Velocity0y.^2 + Velocity0z.^2 );


Velocity = Velocity0(Boltzmann_ind);
Velocity0x = Velocity0x(Boltzmann_ind);
Velocity0y = Velocity0y(Boltzmann_ind);
Velocity0z = Velocity0z(Boltzmann_ind);
% 
% Velocity0x = squeeze(pksp(Boltzmann_ind,1));
% Velocity0y = squeeze(pksp(Boltzmann_ind,2));
% Velocity0z = squeeze(pksp(Boltzmann_ind,3));
% Velocity0 = sqrt(Velocity0x.^2 + Velocity0y.^2 + Velocity0z.^2 );
% clear pksp kq