function Plot_Ef_and_T_dep...
    (fkmu_init,Y_mesh,N_Deig_points,plot_ind,second_par)
%% Plotting
% N_Deig_points = 50;
D_eig_min = 0.;
D_eig_max = 1.5;

Ny = numel(Y_mesh);
if plot_ind == "EFdep"
    filename_to_save = 'Gamma_vs_T.png';
    Y_mesh = 1e3*Y_mesh; % eV to meV
end

D_eig_Boltzmann_mesh = ...
    linspace(D_eig_min,D_eig_max,N_Deig_points);
D_eig_Boltzmann_mesh = repmat(D_eig_Boltzmann_mesh',[1 Ny]);
Y_mesh = repmat(Y_mesh,[N_Deig_points 1]);

pcolor(D_eig_Boltzmann_mesh,Y_mesh,fkmu_init)
shading flat
%shading interp
colormap('jet')
colorbar
%zlim([-Energy_window +Energy_window])
%xlim([0.5,1.5])
xlim([0.0,1.5])
%xlim([0.0,1.])


myfntsize_label = 10;
myfntsize_gca = 10;

myfig = gcf;
set(gcf,'Renderer','painters')

xlabel('\Gamma_i/\Gamma_0  ',...
    'Fontsize',myfntsize_label);


if plot_ind == "Tdep"
    ylabel('T, K   ',...
        'Fontsize',myfntsize_label);
    sprtinf_sec_par = sprintf(...
        'E_F = %0.1f meV',1e3*second_par);
elseif plot_ind == "EFdep"
    ylabel('E_F, meV   ',...
        'Fontsize',myfntsize_label);
    sprtinf_sec_par = sprintf(...
        'T = %0.1f K',second_par);
end
text(0.05,0.1,sprtinf_sec_par,...
    'Units','normalized','Color', 'w',...
    'FontSize',1.5*myfntsize_label)
text(0.8,0.7,'$$|{c_i}(t=0)|$$', 'Color', 'w',...
    'Interpreter', 'latex',...
    'Units','normalized','FontSize',1.5*myfntsize_label);

set(gca,'Fontsize',myfntsize_gca);
set(gca,'LineWidth',1.5)
set(gca,'box','on')

x0 = 10;
y0 = 5;
width = 17.0;%35.0;
height = 12.5;%20.5;%5.5;%9.5;%7.5;%15.0;


set(gcf,'units','centimeters','position',[x0,y0,width,height])
myfig = gcf;
set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


if plot_ind == "Tdep"
    filename_to_save = 'Gamma_vs_T.png';
elseif plot_ind == "EFdep"
    filename_to_save = 'Gamma_vs_EF.png';
end
exportgraphics(myfig,filename_to_save)
%exportgraphics(myfig,'BZ_all_bands.pdf')
end