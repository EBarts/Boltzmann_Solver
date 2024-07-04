function Plot_observable_for_paper3...
    (E_Fermi_array, eps_efficiency_zz,eps_efficiency_xx,...
    M_Edelstein_zz, M_Edelstein_xx,...
    Spin_momentum, Orbital_momentum)

myfntsize = 9;% in points
%4.233 in mm = 0.4233 cm
myfig = gcf;
x0 = 10; y0 = 5;
width = 16.2/2;%16.8/2;%17.0/2; 
height = 10;%9;%6;%6;%8;%7.5;%15.0;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Points','PaperSize',[pos(3), pos(4)])
%set(myfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%set(myfig,'Units','inches');
%set(myfig,'Units','centimeters');


%https://designwizard.com/blog/colour-combination/
%https://www.htmlcsscolor.com/hex/101820
%if eps_ind == "eps_zz"
%color_plt = [16, 24, 32]/255; % 101820 in rgb
color_plt_zz = [40, 2, 116]/255;  % FE7A36  in rgb
%color_plt = [54, 82, 173]/255;  % FE7A36  in rgb
%str_color = '#FF0000';

%color_plt = [95, 67, 30]/255;
color_plt_xx =  [254, 122, 54]/255;% FE7A36 in rgb

%F95700FF
% Convert color code to 1-by-3 RGB array (0~1 each)
%color_plt = sscanf(str_color(2:end),'%2x%2x%2x',[1 3])/255;
%color_plt 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1)

hplt(1) = plot(1e3*E_Fermi_array,M_Edelstein_zz,...
    'Color',color_plt_zz);
hold on
hplt(2) = plot(1e3*E_Fermi_array,M_Edelstein_xx,...
    'Color',color_plt_xx);
NMR_exp = yline(1.5,'--','NMR experiment');
NMR_exp.Color = [.80 0 .40];
NMR_exp.LabelHorizontalAlignment = 'center';
%NMR_exp.LabelVerticalAlignment = 'middle';
hold off

yl1 = ylabel('$$M, 10^{-8} \mu_{\rm B}$$ per Te', ...
    'Interpreter', 'latex');
pos_yl1 = get(yl1,'Pos');

%ylim([0.5,1.5])
ylim([0.5,2.0])
%text(0.45,0.75,'$$ M_{z}$$ ',...
text(0.45,0.4,'$$ M_{z}$$ ',...
    'Interpreter', 'latex', ...
    'Units','normalized','Color',color_plt_zz);
text(0.05,0.2,'$$ M_{x}$$ ',...
    'Interpreter', 'latex',  ...
    'Units','normalized','Color',color_plt_xx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2)
hplt(3) =plot(1e3*E_Fermi_array,eps_efficiency_zz,...
    'Color',color_plt_zz);
hold on
hplt(4) = plot(1e3*E_Fermi_array,eps_efficiency_xx,...
    'Color',color_plt_xx);
hold off

%ylabel('$$\varepsilon^{\rm EE}$$',...
yl2 = ylabel('Efficiency',...
    'Interpreter', 'latex');
pos_yl2 = get(yl2,'Pos');

text(0.45,0.65,'$$ \varepsilon_{zz}$$ ',...
    'Interpreter', 'latex', ...
    'Units','normalized','Color',color_plt_zz);
text(0.05,0.25,'$$ \varepsilon_{xx}$$ ',...
    'Interpreter', 'latex', ...
    'Units','normalized','Color',color_plt_xx);
%    'Units','pixels',...
% xlabel('$$E_F$$, {\rm meV}',...
%     'Interpreter', 'latex');
%ylim([0,0.6])
ylim([0,1.0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,3)
arb_momentum_units = 1e4;

hplt(5)=plot(1e3*E_Fermi_array,...
    Spin_momentum/arb_momentum_units,'Color',color_plt_zz);
%ylim([-1.5,0.5])
hold on
hplt(6)=plot(1e3*E_Fermi_array,...
    Orbital_momentum/arb_momentum_units,'Color','b');
hold off
%ylim([-1.2,-0.05])

%ylabel('$$\langle S_z, L_z \rangle$$',...
yl3 = ylabel('Accumulation, a.u.',...
    'Interpreter', 'latex');
pos_yl3 = get(yl3,'Pos');

set(yl1,'Pos',[pos_yl3(1) pos_yl1(2) pos_yl1(3)]);

%set(yl2,'Pos',[pos_yl3(1) pos_yl2(2) pos_yl2(3)]);

%    'Units','pixels',...
xlabel('$$E_F$$, {\rm meV}',...
    'Interpreter', 'latex');
% ylim([0,0.6])
% ylim([0,1.0])
text(0.15,0.85,'$$\langle S_z \rangle$$ ',...
    'Interpreter', 'latex',...
    'Units','normalized','Color',color_plt_zz);

text(0.5,0.25,'$$\langle L_z \rangle$$ ',...
    'Interpreter', 'latex',...
    'Units','normalized','Color','b');

for isbplt = 1:6
    hplt(isbplt).LineWidth = 1.7;
end

subplot(3,1,1)
text(0.6,0.15,'$$j_{\rm el}=100 \, {\rm A} \, {\rm cm}^{-2}$$ ',...
    'Interpreter', 'latex',...
    'Units','normalized');


subplot(3,1,1)
text(0.02,0.9,'(a)',...
    'Units','normalized',...
    'FontSize',1.1*myfntsize)

subplot(3,1,2)
text(0.02,0.9,'(b)',...
    'Units','normalized',...
    'FontSize',1.1*myfntsize)

subplot(3,1,3)
text(0.02,0.9,'(c)',...
    'Units','normalized',...
    'FontSize',1.1*myfntsize)

for isbplt = 1:3
    subplot(3,1,isbplt)
    % get handle to current axes
    a = gca;
    %
    set(a,'Fontsize',1.0*myfntsize);
    %set(a,'LineWidth',1.)
    % set box property to off and remove background color
    set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
    b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
    axes(a)
    % link axes in case of zooming
    linkaxes([a b])


end
ylim([-1.4,0.1])

%     xpos = -33; % (find this out from get(yl,'pos') on the desired label x-location)
%     %yl = ylabel('Label Here');
%     yl = get(a,'Ylabel');
%     pos=get(yl,'Pos');
%     set(yl,'Pos',[xpos pos(2) pos(3)]);

set(gcf,'Renderer','painters')

filename_to_save = ...
    'Figures_for_paper/Observables_vs_EF_for_paper';

exportgraphics(myfig,[filename_to_save,'.eps'],...
    'BackgroundColor','none','ContentType','vector')
%exportgraphics(myfig,[filename_to_save,'.pdf'])

end