function [abs_fkmu, observables_vs_EF,...
    Gamma_at_Ef_binned, Density_of_States_Ef] = ...
    gimme_observables(E_Fermi,TBoltz_K,N_Deig_points,model_ind)
    is_single_run = 0;
    Boltzmann_Solver;
    find_fkmu;
    observables_vs_EF = zeros(1,9);
    
    % a b c d coefficients
    % val(x) = a*exp(b*x) + c*exp(d*x)
%     observables_vs_EF(1:4) = -exp_fit2_coef(:);
%     if abs(observables_vs_EF(2)) <...
%         abs(observables_vs_EF(4))
%         observables_vs_EF(1:4)=...
%             observables_vs_EF([3,4,1,2]);
%     end
    %observables_vs_EF(1:2) = exp_fit2_coef(:);
    observables_vs_EF(1) = exp_fit2_coef(1);
    observables_vs_EF(2) = -exp_fit2_coef(2);
    observables_vs_EF(3) = spin_conversion_efficiency(1);
    observables_vs_EF(4) = RE_coeff_xx;
    observables_vs_EF(5) = RE_coeff_yy;
    observables_vs_EF(6) = RE_coeff_zz;
    observables_vs_EF(7:12) = total_SL_at_t0;
    observables_vs_EF(13) = total_el_current_x;
    observables_vs_EF(14) = total_el_current_y;
    observables_vs_EF(15) = total_el_current_z;
%     observables_vs_EF(5) = spin_conversion_efficiency(1);
%     observables_vs_EF(6) = RE_coeff;
%     observables_vs_EF(7) = SHE_coeff;
%     observables_vs_EF(8) = total_el_current_z;
%     observables_vs_EF(9) = total_spin_current_zz;
    %
    abs_fkmu = delta_fkmu_init_decomp_coef1;
    Gamma_at_Ef_binned = sum(Gamma_at_Ef(:));
    %abs(delta_fkmu_init_decomp_coef(set_of_eigs));
end