figure
delta_arb_units = max(abs(delta_fkmu_init_decomp_coef1(:)));%2e3;
% hplot = plot(D_eig_Boltzmann_array1',...
%      abs(delta_fkmu_init_decomp_coef1)','.-');
%hplot = plot(D_eig_Boltzmann_array1',...
%         abs(delta_fkmu_init_decomp_coef1)'/delta_arb_units,'*-');
%     abs(delta_fkmu_init_decomp_coef1)','.-');
data_X = D_eig_Boltzmann_array1';
data_Y = abs(delta_fkmu_init_decomp_coef1)'/delta_arb_units;
%histogram(data_Y, 'BinEdges',D_eig_Boltzmann_array1);
ybars = [0 0.8];
%plot(data_X,data_Y,'LineWidth',1.7,'LineStyle','*')

area(data_X,data_Y,'FaceColor',[0.8500 0.3250 0.0980])
hold on
hplot = scatter(data_X,data_Y,'.k');
%patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8])
%plot(x, y, 'bp')
hold off
%axis([0  10    0  10])

%hplot = plot(D_eig_Boltzmann_array1',...
%         abs(delta_fkmu_init_decomp_coef1)'/delta_arb_units,'*-');


%xlim([0, 1.5])
%xlim([0, 1.0])
%ylim([0., 1.2])
% xlim([0.15, 0.5]) % for Fig.3
xlim([0., 1.0]) % for Fig.4
hplot.LineWidth = 1.7;

%return
myfntsize = 9;

myfig = gcf;
set(gcf,'Renderer','painters')

ylabel('$$A_i^{\rm EE}$$  ',...
        'Interpreter', 'latex',...
    'Fontsize',myfntsize);

xlabel('$$\Gamma_i$$  ',...
        'Interpreter', 'latex',...
    'Fontsize',myfntsize);


% text(0.03,0.85,'$$E_F = -20$$ meV',...
%     'Interpreter', 'latex',...
%     'Units','normalized',...
%     'FontSize',1.*myfntsize)
% 
% text(0.75,0.8,'fast',...
%     'Interpreter', 'latex',...
%     'Units','normalized',...
%     'FontSize',1.*myfntsize)
% 
% text(0.2,0.47,'slow',...
%     'Interpreter', 'latex',...
%     'Units','normalized',...
%     'FontSize',1.*myfntsize)
% 
% text(-0.25,1.05,'(c)',...
%     'Units','normalized',...
%     'FontSize',1.*myfntsize)


%set(gca,'Fontsize',myfntsize);
%set(gca,'LineWidth',1.)
%set(gca,'box','on')



x0 = 10;
y0 = 5;
width = 17.0/3;%17.0/2/2;
height = 3;%8;%7.5;%15.0;

% width = 17.0;%35.0;
% height = 12.5;%20.5;%5.5;%9.5;%7.5;%15.0;


set(gcf,'units','centimeters','position',[x0,y0,width,height])
myfig = gcf;
%set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

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

%xlim([0.15, 0.7]) % for Fig.3
xlim([0.0, 1.0]) % for Fig.4

filename_to_save = 'Figures_for_paper_kp/spectrum_for_paper';
exportgraphics(myfig,[filename_to_save,'.pdf'])

save(filename_to_save,'myfig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
% hplot = plot(D_eig_Boltzmann_array1',...
%      abs(delta_fkmu_init_decomp_coef1)','.-');
hplot = plot(time_array,spin_conversion_efficiency);
%ylim([0,0.5])
ylim([0,1.0])

xlim([0, 10])
hplot.LineWidth = 1.7;

myfig = gcf;
set(gcf,'Renderer','painters')

ylabel('$$\varepsilon^{\rm EE}(t)$$  ',...
        'Interpreter', 'latex',...
    'Fontsize',myfntsize);

xlabel('$$ t, \tau_0$$  ',...
        'Interpreter', 'latex',...
    'Fontsize',myfntsize);

text(0.25,0.65,'$$\sim {\rm exp}\left(-\Gamma_{\rm spin}\frac{t}{\tau_0}\right)$$',...
    'Interpreter', 'latex',...
    'Color','k','Units','normalized',...
    'FontSize',1.*myfntsize)

text(-0.25,1.05,'(d)',...
    'Units','normalized',...
    'FontSize',1.*myfntsize)




set(gcf,'units','centimeters','position',[x0,y0,width,height])
myfig = gcf;
%set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

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

filename_to_save = 'Figures_for_paper_kp/spinvstime_for_paper';
exportgraphics(myfig,[filename_to_save,'.pdf'])

%save(filename_to_save,'myfig');

