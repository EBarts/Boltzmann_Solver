function Plot_func_fkmu(fkmu_list, isblt)

if isblt ~=0
    %subplot(2,3,isblt);
    subplot(1,3,isblt);
end
set(gcf,'Renderer','painters')
%% Plotting
% load file first
%load('data_kp_model/bandsKP2D.mat')
load('PAOFLOW_TB_model/bands_temp2D.mat')

Nx = 201;
Nz = Nx;
Ny = 1;
kx_window = 0.03;%0.05; % in inverse Angstrem
ky_window = 0.03;%0.05; 
kz_window = 0.1;%0.07; 
if Ny == 1
    ky_window = 0;
end
[kx_mesh,ky_mesh, kz_mesh] = ndgrid(...
    linspace(-kx_window,+kx_window,Nx),...
    linspace(-ky_window,+ky_window,Ny),...
    linspace(-kz_window,+kz_window,Nz));

myfntsize_label = 10;
myfntsize_gca = 10;

%fkmu_list = V_eig_Boltzmann(:, Eig_fkmu_ind);   
% for 2D:
%fkmu = reshape(fkmu_list, Nx, Ny, Nz);
fkmu = zeros(Nx, Ny, Nz);
fkmu(Boltzmann_ind) = fkmu_list;

Boltzmann_ind1 = abs(Energy_E_F_all)< 0.3*Energy_window;
fkmu1 = zeros(Nx, Ny, Nz);
fkmu1(Boltzmann_ind1) = fkmu(Boltzmann_ind1);

iNy = 1;
fkmu2D = squeeze(fkmu1(:,iNy,:));

X = squeeze(kx_mesh(:,iNy,:));
Z = squeeze(kz_mesh(:,iNy,:));
%Energy_E_F_all2D = squeeze(Energy_E_F_all(:,iNy,:));
Energy_E_F_all2D = squeeze(reshape(Energy_E_F_all,[Nx,iNy,Nz]));
%pcolor
% transpose because we need (Nz, Nx) array for plotting
surf(X,Z,Energy_E_F_all2D,fkmu2D);
%shading flat
colormap('jet')
colorbar
zlim([-Energy_window +Energy_window])
%E_Fermi = -0.005;%-0.02;
%Energy_E_F_all = energy_kp_plus - E_Fermi;
%Energy_window = 3e-3;%5e-3; % 5*1e-3;
return

%% spin projected
%calc_ind = "Lz";
Sz_matrix_el = fkmu;%matrix_el_calc(V_eig_array, Nz, Nbands, calc_ind);
%mat_el_fin = zeros(Nz, Nx, Nbands);

%Sz_matrix_el_b = zeros(Nz, Nx);
for iband = 11:12 %1:Nbands %7:12
    Sz_matrix_el_b = Sz_matrix_el(:,:, iband-10);
    D_eig = DeltaLambda*1e-3*...
        (D_eig_array(:,:, iband) - max(max(D_eig_array(:,:, 12)))); % in eV
    hold on
    
    surf(X,Z,D_eig,Sz_matrix_el_b);
    shading flat
    zlim([-0.15 -0.0])
    hold off
%     hold on
%     for ikzplt = 1: Nz
%         if arrayE_test(ikzplt,iband) ~= 0
%         scatter(kz_lattice(ikzplt),DeltaLambda*1e-3*...
%                 (arrayE_test(ikzplt,iband) - max(D_eig_array(:, 12))),'xk')
%         end
%     end
%     hold off
    %patch(kz_lattice, D_eig, EdgeColor ='interp')
    %plot(kz_lattice, D_eig, Sz_matrix_el_b, 'k','LineWidth',1.5)
end

if max(imag(Sz_matrix_el(:))) >1e-4
    fprintf('ERROR: the matrix elements are imaginary: %d\n',...
        max(imag(Sz_matrix_el(:))));
end
colormap('jet')
colorbar


% hold on
% scatter(DeltaLambda*1e-3*...
%         (arrayE_test(:, iband) - max(D_eig_array(:, 12))),'LineWidth',5.2)
% hold off

yaxis_point = -0.1;%-0.045;
xGamma_point = 0.0;
xZ_point = 0.5 - 0.01;
% text(xGamma_point,yaxis_point, '$$ \Gamma $$ ','Units','normalized',...
%     'FontSize',11,'FontWeight','bold','Interpreter', 'latex');
% text(xZ_point,yaxis_point, '$$ Z $$','Units','normalized',...
%     'FontSize',11,'FontWeight','bold', 'Interpreter', 'latex');



xlabel('k_x   ',...
    'Fontsize',myfntsize_label);

ylabel('k_z  ',...
    'Fontsize',myfntsize_label);

set(gca,'Fontsize',myfntsize_gca);
set(gca,'LineWidth',1.5)
set(gca,'box','on')

x0 = 10;
y0 = 5;
width = 35.0;
height = 20.5;%5.5;%9.5;%7.5;%15.0;

set(gcf,'units','centimeters','position',[x0,y0,width,height])
myfig = gcf;
set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%exportgraphics(myfig,'BZ_all_bands.pdf')
%exportgraphics(myfig,'BZ_all_bands.png')
end