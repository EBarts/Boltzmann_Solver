function c_max = Plot_func_fkmu_kp(fkmu_list,isblt,...
    Energy_window,Boltzmann_ind,Energy_E_F_all,Nz, c_max0,plot_ind)
figure
if isblt ~=0
    %subplot(2,3,isblt);
    subplot(1,3,isblt);
end
set(gcf,'Renderer','painters')
%% Plotting
% load file with energy bands first
%load('model_PAOFLOW_TB/Te_data_61x61x61_kxy=0_06_kz=0_12.mat','ek','kq');
%load('model_PAOFLOW_TB/Te_data_101x101x101_kxy=0_06_kz=0_12.mat','ek','kq');
% ek size (Nx*Ny*Nz,1)
% kq size (Nx*Ny*Nz,3)
%kq1 = reshape(kq,[Nz,Ny,Nx,3]);
%Energy_E_F_all = Energy_E_F_all*1e3; % in meV

E_Fermi = -20*1e-3;
addpath('core_codes_kp')  

Introduce_kp_model


energy_list = energy_list - max(energy_list(:));
%ek = ek - max(ek(:));



E_Fermi = -20;
% Energy_E_F_all = ek*1e3 - E_Fermi; % in meV
% Energy_window = Energy_window *1e3;

Energy_E_F_all = energy_list*1e3 - E_Fermi; % in meV
Energy_window = Energy_window *1e3;


Nx = Nz;
Ny = Nz;
% kx_window = 0.05;%0.05; % in inverse Angstrem
% ky_window = 0.05;%0.05;
% kz_window = 0.1;%0.07;
% if Ny == 1
%     ky_window = 0;
% end

myfntsize = 9;

%fkmu_list = V_eig_Boltzmann(:, Eig_fkmu_ind);   
% for 2D:
%fkmu = reshape(fkmu_list, Nx, Ny, Nz);
%fkmu = zeros(Ny,Nz, Nx);
fkmu = zeros(Nx, Ny, Nz);
fkmu(Boltzmann_ind) = fkmu_list;

% energy_list = zeros(Ny,Nz, Nx);
%Boltzmann_ind1 = abs(Energy_E_F_all)< 0.3*Energy_window;
Boltzmann_ind1 = abs(Energy_E_F_all)< ...
    Energy_window;
    %0.5*Energy_window;
fkmu1 = zeros(Nx, Ny, Nz);
fkmu1(Boltzmann_ind1) = fkmu(Boltzmann_ind1);

%Energy_E_F_all = reshape(Energy_E_F_all,[Nz,Ny,Nx]);
iNy = round(Ny/2); % 1;


%kq_all = reshape(kq,[Nz,Ny,Nx,3]);
%Xk = squeeze(kq_all(:,iNy,:,1));
%Zk = squeeze(kq_all(:,iNy,:,3));
Xk = squeeze(kx_mesh(:,iNy,:));
Zk = squeeze(kz_mesh(:,iNy,:));



%fkmu2D = squeeze(fkmu1(iNy,:,:));
%fkmu2D = squeeze(fkmu1(:,iNy,:));
fkmu2D = squeeze(...
    (fkmu1(:,iNy,:) + fkmu1(:,iNy+1,:) + fkmu1(:,iNy-1,:))/3 ...
    );

%Ny,Nz, Nx
%[X,Z] = meshgrid(1:Nx,1:Nz);
[X,Z] = meshgrid(-(Nx-1)/2:(Nx-1)/2,...
    -(Nz-1)/2:(Nz-1)/2);
%Energy_E_F_all2D = squeeze(Energy_E_F_all(:,iNy,:));
%Energy_E_F_all2D = squeeze(reshape(Energy_E_F_all,[Nx,iNy,Nz]));
%Energy_E_F_all2D = reshape(Energy_E_F_all,[Nx,Ny,Nz]);
%Energy_E_F_all2D = squeeze(Energy_E_F_all2D(:,iNy,:));
Energy_E_F_all2D = squeeze(Energy_E_F_all(:,iNy,:));
% pcolor
% transpose because we need (Nz, Nx) array for plotting

%surf(X,Z,Energy_E_F_all2D,fkmu2D);

%surf(Xk,Zk,Energy_E_F_all2D,fkmu2D);
%contour(X,Z,Energy_E_F_all2D,'k');

%hold on
%pcolor(X,Z,fkmu2D);
%hold off

hold on

% pcolor(X,Z,fkmu2D);
% contour(X,Z,Energy_E_F_all2D,'k');
pcolor(Xk,Zk,fkmu2D);
contour(Xk,Zk,Energy_E_F_all2D,4,'k');
%contour(Xk,Zk,Energy_E_F_all2D,'k');
hold off

%shading interp
shading flat

%colormap('jet')
%colormap(jet(1000))
% zmin = min(fkmu2D(:));
% zmax = max(fkmu2D(:));
% colormap(redwhiteblue(zmin, zmax));

%colormap(redwhiteblue)
%colormap('parula')
%colormap(seismic(128))
%colormap(redblue)
colormap(seismic)
%caxis([-0.5 0.5])
% 
c_max = max(abs(fkmu2D(:)));
%caxis([-1 1])

if c_max0 == 0
    caxis([-c_max c_max])
else
    caxis([-c_max0 c_max0])
end
caxis([-c_max c_max])
%c_max = 0;
colorbar;


x0 = 10;
y0 = 5;
width = 17.0/3;
height = 2;%10/3;%12/3;%6/2;%8;%7.5;%15.0;

% xlim([100, 300])
% ylim([50, 400])

% xlabel('k_x   ',...
%     'Fontsize',myfntsize_label);
% 
% ylabel('k_z  ',...
%     'Fontsize',myfntsize_label);


% if plot_ind == "slow"
% 
% end
text(0.4,0.7,char(plot_ind),...
    'Interpreter', 'latex',...
    'Units','normalized',...
    'FontSize',myfntsize);
%text(0.1,0.92,'total',...
% text(0.4,0.7,'PST',...
%     'Interpreter', 'latex',...
%     'Units','normalized',...
%     'FontSize',1.1*myfntsize)

% text(0.45,0.5,'+',...
%     'Interpreter', 'latex',...
%     'Units','normalized',...
%     'FontSize',1.1*myfntsize_label)
% 
% text(0.4,0.3,'slow',...
%     'Interpreter', 'latex',...
%     'Units','normalized',...
%     'FontSize',1.1*myfntsize_label)

if plot_ind == "total"
    label_fkmu = '(e)'; 
elseif plot_ind == "fast"
    label_fkmu = '(f)';
elseif plot_ind == "slow"
    label_fkmu = '(g)';
end
text(-0.2,1.0,label_fkmu,...
    'Units','normalized',...
    'FontSize',1.*myfntsize)

set(gca,'Fontsize',myfntsize);
%set(gca,'LineWidth',1.)
set(gca,'box','on')

% width = 20.0;
% height = 20.;%5.5;%9.5;%7.5;%15.0;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
myfig = gcf;
set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
pdf_name = ['Figures_for_paper_kp/fkmu_contour_',char(plot_ind),'.pdf'];
%exportgraphics(gcf,'Figures_for_paper/fkmu_contour.pdf');

%mat_name = ['Figures_for_paper_kp/fkmu_contour_',char(plot_ind)];
%save(mat_name,'myfig');

exportgraphics(gcf,pdf_name);
end