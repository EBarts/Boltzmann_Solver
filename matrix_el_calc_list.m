% Calculation of single electron energy bands of tellurium states
% over the 1st Brillouin zone
% for 1D helical chain of elemental Te
% based on the codes for calculation of iodine p-states in CrI_3 compound
% Evgenii Barts, University of Groningen, The Netherlands
% 20 July, 2023

function mat_el_fin = matrix_el_calc_list...
    (V_eig_array, calc_ind,model_ind)
[Nstates, Npoints] = size(V_eig_array);
if model_ind == "paoflowTB_model" % reshuffle wave functions
    for isublat = 1:3
        array_orb = 6*(isublat - 1) + (1:6);
        array_paoflow = 6*(isublat - 1) + [6,5,4,3,2,1];
        V_eig_array(array_orb,:) =...
            V_eig_array(array_paoflow,:);
    end
end
V_eigk = zeros(Nstates,1);
%V_eigk = zeros(Nbands,Nbands); % eigen states for a single k point
mat_el_fin = zeros(Npoints, 1);
%char_calc_ind = char(calc_ind);
for ikt = 1:Npoints
    V_eigk(:) = V_eig_array(:,ikt);
    %Transformation matrix to |J,Jz> basis
    %%
    U_matrix = [-sqrt(1/2),0,sqrt(1/6),        0,0,-sqrt(1/3); ...
                -1i*sqrt(1/2),0,-1i*sqrt(1/6), 0,0,1i*sqrt(1/3); ...
                0,sqrt(2/3),0,                 0,-sqrt(1/3),0; ...
                ...
                0,-sqrt(1/6),    0,            sqrt(1/2),-sqrt(1/3),0; ...
                0,-1i*sqrt(1/6), 0,            -1i*sqrt(1/2),-1i*sqrt(1/3),0; ...
                0,0,     sqrt(2/3),            0, 0,sqrt(1/3)   ];

    %Tensor product
    U_matrix_all_ligands = kron(eye(3), U_matrix);
    %
    V_iegkib = zeros(Nstates,1);

    if calc_ind == "Sz"
        sigma_z_matrix = [eye(3),  0*eye(3);...
                    0*eye(3), -eye(3)];
        sigma_z_all_ligands = kron(eye(3), sigma_z_matrix);
        matrix_to_calc = sigma_z_all_ligands/2; %Sz = sigma_z/2
    elseif calc_ind == "Sx"
        sigma_x_matrix = [0*eye(3),  1*eye(3);...
                    1*eye(3), 0*eye(3)];
        sigma_x_all_ligands = kron(eye(3), sigma_x_matrix);
        matrix_to_calc = sigma_x_all_ligands/2;
    elseif calc_ind == "Sy"
        sigma_y_matrix = [0*eye(3),  -1i*eye(3);...
                    1i*eye(3), 0*eye(3)];
        sigma_y_all_ligands = kron(eye(3), sigma_y_matrix);
        matrix_to_calc = sigma_y_all_ligands/2;
    elseif calc_ind == "Lz"
        %%
        Lz_plus1 = -1/sqrt(2)*[1, 1i, 0].'; %column vector
        Lz_min1 = 1/sqrt(2)*[1, -1i, 0].';
        Lz = Lz_plus1*Lz_plus1' - Lz_min1*Lz_min1';
        Lz_matrix = kron(eye(2), Lz);
        Lz_all_ligands = kron(eye(3), Lz_matrix);
        matrix_to_calc = Lz_all_ligands; %
    elseif calc_ind == "Lx"
        %%
        Lx_plus1 = -1/sqrt(2)*[0, 1, 1i].'; %column vector
        Lx_min1 = 1/sqrt(2)*[0, 1, -1i].';
        Lx = Lx_plus1*Lx_plus1' - Lx_min1*Lx_min1';
        Lx_matrix = kron(eye(2), Lx);
        Lx_all_ligands = kron(eye(3), Lx_matrix);
        matrix_to_calc = Lx_all_ligands; %
    elseif calc_ind == "Ly"
        %%
        Ly_plus1 = -1/sqrt(2)*[1i, 0, 1].'; %column vector
        Ly_min1 = 1/sqrt(2)*[-1i, 0, 1].';
        Ly = Ly_plus1*Ly_plus1' - Ly_min1*Ly_min1';
        Ly_matrix = kron(eye(2), Ly);
        Ly_all_ligands = kron(eye(3), Ly_matrix);
        matrix_to_calc = Ly_all_ligands; %
    elseif calc_ind == "pz"
        %%
        Pz_vec = [0, 0, 1].'; %column vector
        Pz = Pz_vec*Pz_vec';
        Pz_matrix = kron(eye(2), Pz);
        Pz_all_ligands = kron(eye(3), Pz_matrix);
        matrix_to_calc = Pz_all_ligands; %
    elseif calc_ind == "Triangle_plus" ||...
            calc_ind == "Triangle_min" || calc_ind == "Triangle_s"
        %
        if calc_ind == "Triangle_plus"
            Triangle_vec = 1/sqrt(3)*[1, exp(2*pi*1i/3), exp(-2*pi*1i/3)].'; %column vector
        elseif calc_ind == "Triangle_min"
            Triangle_vec = 1/sqrt(3)*[1, exp(-2*pi*1i/3), exp(2*pi*1i/3)].';
        elseif calc_ind == "Triangle_s"
            Triangle_vec = 1/sqrt(3)*[1, 1, 1].';
        end
        Triangle_projector = Triangle_vec*Triangle_vec';
        Triangle_all_ligands = kron(Triangle_projector, eye(6));
        matrix_to_calc = Triangle_all_ligands;
    end

    V_iegkib(:) = V_eigk(:);
    V_iegkib_transformed = U_matrix_all_ligands * V_iegkib;
    mat_el_fin(ikt) = V_iegkib_transformed'*matrix_to_calc*V_iegkib_transformed;

end

end