function [orbit_matrix,cost,delta_v_mag]= orbit_solver(r1_vec,r2_vec)

    [char_star] = load_characteristic_values();

    orbit_matrix = lambert_solver(r1_vec,r2_vec);

    %WORKTASK 1 - VISUALIZE ORBIT MATRIX SOLUTIONS TO MAKE SURE THEY MAKE
    %SENSE
    
    %WORKTASK 3A - FUNCTIION HEADERS AND LOGIC TREE: EACH FUNCTION USED IN
    %THIS SOLVER SHOULD HAVE A FUNCTION HEADER AND BE LABELLED IN THE
    %FUNCTION TREE

    %WORKTASK 3B - VARIABLE TITLES: EACH INPUT AND OUTPUT VARIABLE SHOULD
    %HAVE ITS INFORMATION AND UNIT IN THE FUNCTION HEADER

    %WORKTASK 2 - UNDERSTANDING THE LAMBERT PROBLEM

    %magnitudes of position vectors

    r1 = norm(r1_vec);
    r2 = norm(r2_vec);

    load_characteristic_values();
    VE = 29.78./char_star.v;
    %write delta v fn separately
    [RA, RB, VA, VB, delVAmag, delVBmag, delVA, delVB] = deltaV_ver2(r1_vec, r2_vec, r1,...
        orbit_matrix(:,2), orbit_matrix(:,7), VE) %orbit_matrix(..2) = e_solutions, orbit_matrix(..7) = p_solutions


    orbital_period = sqrt((orbit_matrix(:,1).^3)./((char_star.mu)./(4*pi^2)));

    delta_v_mag = sqrt(sum(delta_v.^2,2));

    cost = orbital_period .* delta_v_mag;

end