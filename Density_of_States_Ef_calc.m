function Density_of_States_Ef = Density_of_States_Ef_calc...
    (energy_list, E_Fermi, Energy_window ,E_width,delta_E_cutoff)

Energy_E_F_all = energy_list - max(energy_list(:)) - E_Fermi;

% Number of states near Ef
Boltzmann_ind = abs(Energy_E_F_all)< Energy_window;
energy_E_F = Energy_E_F_all(Boltzmann_ind);

Nsites = max(size(energy_E_F));

DeltaDirac_func = @(X) E_width/(E_width^2 + X^2)/pi;

% Calculate the density of states at E Fermi
Density_of_States_Ef = 0.0;

for n = 1:Nsites
    absdeltaE = abs(energy_E_F(n));
    if absdeltaE < delta_E_cutoff
        Density_of_States_Ef = Density_of_States_Ef +...
             DeltaDirac_func(absdeltaE);
    end
end

% test 2:
if Density_of_States_Ef < 1e-12
    disp(['Hello test 2 = %f',num2str(Density_of_States_Ef)]);
end


end