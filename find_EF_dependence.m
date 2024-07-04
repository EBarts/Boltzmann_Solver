%% 
addpath('core_codes')  
N_Deig_points = 50;

E_FN = 24;

E_Fermi_array = linspace(-0.0001,-0.03,E_FN);

fkmu_init_vs_EF = zeros(N_Deig_points,E_FN);
observables_vs_EF = zeros(15,E_FN);
Gamma_at_Ef_total = zeros(1,E_FN);
Density_of_States_Ef = zeros(1,E_FN);
TBoltz_K = 10;
%
model_ind = "paoflowTB_model";
tic
parfor iEF = 1:E_FN
    [fkmu_init_vs_EF(:,iEF),observables_vs_EF(:,iEF),...
        Gamma_at_Ef_total(iEF), Density_of_States_Ef(iEF)] = ...
        calc_Boltzmann_observables(E_Fermi_array(iEF),...
        TBoltz_K,N_Deig_points,model_ind);
iEF 
end
toc

Filename_to_save = [char(model_ind),'_run_' ,datestr(now,'mm-dd-yyyy HH-MM')];
save(Filename_to_save);

%%
figure
Plot_Ef_and_T_dep(fkmu_init_vs_EF,...
    E_Fermi_array,N_Deig_points,"EFdep",TBoltz_K)

%%
figure
subplot(6,1,1)

hplt(1)=plot(1e3*E_Fermi_array,observables_vs_EF(3,:));
ylabel('Efficiency');

subplot(6,1,2)
%1./
hplt(2) = plot(1e3*E_Fermi_array,observables_vs_EF(2,:));
ylabel('S_z Relax. rate/\Gamma_0');

subplot(6,1,3)
hplt(3) = plot(1e3*E_Fermi_array,observables_vs_EF(1,:));
ylabel('Amplitude, a.u.');
xlabel('E_F, meV');

subplot(6,1,4)
hplt(4)=plot(1e3*E_Fermi_array,observables_vs_EF(4,:));
ylabel('<Mx>, 1e-8 \mu_B');

subplot(6,1,5)
hplt(5) = plot(1e3*E_Fermi_array,observables_vs_EF(5,:));
ylabel('<My>, 1e-8 \mu_B');

subplot(6,1,6)
hplt(6) = plot(1e3*E_Fermi_array,observables_vs_EF(6,:));
ylabel('<Mz>, 1e-8 \mu_B');

for isbplt = 1:6
    subplot(6,1,isbplt)
    set(gca,'Fontsize',10);
    set(gca,'LineWidth',1.5)
    set(gca,'box','on') 
end

for iplt = 1:6
    hplt(iplt).LineWidth = 2.0;
end

x0 = 10;
y0 = 5;
width = 20.0;
height = 20.;%5.5;%9.5;%7.5;%15.0;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
myfig = gcf;
set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

exportgraphics(gcf,'Observables_vs_EF.png')

%%
figure
subplot(6,1,1)

hplt(1)=plot(1e3*E_Fermi_array,observables_vs_EF(7,:));
ylabel('Sx');

subplot(6,1,2)
%1./
hplt(2) = plot(1e3*E_Fermi_array,observables_vs_EF(8,:));
ylabel('Sy');

subplot(6,1,3)
hplt(3) = plot(1e3*E_Fermi_array,observables_vs_EF(9,:));
ylabel('Sz');
xlabel('E_F, meV');

subplot(6,1,4)
hplt(4)=plot(1e3*E_Fermi_array,observables_vs_EF(10,:));
ylabel('Lx');

subplot(6,1,5)
hplt(5) = plot(1e3*E_Fermi_array,observables_vs_EF(11,:));
ylabel('Ly');

subplot(6,1,6)
hplt(6) = plot(1e3*E_Fermi_array,observables_vs_EF(12,:));
ylabel('Lz');

for isbplt = 1:6
    subplot(6,1,isbplt)
    set(gca,'Fontsize',10);
    set(gca,'LineWidth',1.5)
    set(gca,'box','on') 
end

for iplt = 1:6
    hplt(iplt).LineWidth = 2.0;
end

x0 = 10;
y0 = 5;
width = 20.0;
height = 20.;%5.5;%9.5;%7.5;%15.0;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
myfig = gcf;
set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

exportgraphics(gcf,'Observables_vs_EF_SL.png')

%%
figure
subplot(3,1,1)

hplt(1)=plot(1e3*E_Fermi_array,observables_vs_EF(13,:));
ylabel('jx');

subplot(3,1,2)
%1./
hplt(2) = plot(1e3*E_Fermi_array,observables_vs_EF(14,:));
ylabel('jy');

subplot(3,1,3)
hplt(3) = plot(1e3*E_Fermi_array,observables_vs_EF(15,:));
ylabel('jz');
xlabel('E_F, meV');


for isbplt = 1:3
    subplot(3,1,isbplt)
    set(gca,'Fontsize',10);
    set(gca,'LineWidth',1.5)
    set(gca,'box','on') 
end

for iplt = 1:3
    hplt(iplt).LineWidth = 2.0;
end

x0 = 10;
y0 = 5;
width = 20.0;
height = 20.;%5.5;%9.5;%7.5;%15.0;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
myfig = gcf;
set(myfig,'Units','Inches');
pos = get(myfig,'Position');
set(myfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

exportgraphics(gcf,'Observables_vs_EF_currents.png')