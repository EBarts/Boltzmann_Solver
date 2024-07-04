function [BoltzmannKernel, Gamma_at_Ef, Density_of_States_Ef]...
    = BoltzmannKernel_calc_rel_time...
    (energy_E_F, wave_functions, E_width,delta_E_cutoff, Density_of_States0)

[Nbands,Nsites] = size(wave_functions);
tao_W = zeros(Nsites,Nsites);
wave_function_in = zeros(Nbands,1);
wave_function_out = zeros(Nbands,1);

% Calculate the (tao * W) part of the kernel
% F_n = tao_W(n,m) * F_m
% weight the energy delta function with an envelope function
DeltaDirac_func = @(X) E_width/(E_width^2 + X^2)/pi;


% and Calculate analytically
% single-band energy dependent mean relaxation time
Gamma_at_Ef = zeros(Nsites,1);


for n = 1:Nsites
    Gamma_et_Ef_n = 0;
    for m = 1:Nsites
        absdeltaE = abs(energy_E_F(n) - energy_E_F(m));
        if absdeltaE < delta_E_cutoff && n~=m
            wave_function_in(:) = wave_functions(:,m);
            wave_function_out(:) = wave_functions(:,n);
            mat_element = abs(wave_function_out'*wave_function_in)^2;
            tao_W(n,m) = DeltaDirac_func(absdeltaE)*mat_element;
            Gamma_et_Ef_n = Gamma_et_Ef_n + tao_W(n,m);
        end
    end
    Gamma_at_Ef(n) = Gamma_et_Ef_n;
end


% Calculate the density of states at E Fermi
Density_of_States_Ef = 0.0;

for n = 1:Nsites
    absdeltaE = abs(energy_E_F(n));
    if absdeltaE < delta_E_cutoff
        Density_of_States_Ef = Density_of_States_Ef +...
             DeltaDirac_func(absdeltaE);
    end
end

time_scale = "absolute";

if time_scale == "relative"
    Density_of_States_scale = Density_of_States_Ef;
elseif time_scale == "absolute"
    Density_of_States_scale = Density_of_States0;
end

tao_W = tao_W/Density_of_States_scale;
Gamma_at_Ef = Gamma_at_Ef/Density_of_States_scale;

% Construct the complete Boltzmann kernel
BoltzmannKernel = - tao_W;
sum_tao_W = sum(tao_W,2);
for n = 1:Nsites
    BoltzmannKernel(n,n) = BoltzmannKernel(n,n) + sum_tao_W(n);
end


%Density_of_States_Ef = Density_of_States_Ef/Nsites;

% test 2:
if Density_of_States_Ef < 1e-12
    disp(['Hello test 2 = %f',num2str(Density_of_States_Ef)]);
end

end