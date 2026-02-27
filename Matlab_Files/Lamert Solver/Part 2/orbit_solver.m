function [orbit_matrix, cost]= orbit_solver(r1_vec,r2_vec)

    [char_star] = load_charecteristic_values();

    [delta_v,orbit_matrix]=lambert_solver(r1_vec,r2_vec);
   
    orbital_period = sqrt((orbit_matrix(:,1).^3)./((char_star.mu)./(4*pi^2)));

    delta_v_mag = sqrt(sum(delta_v.^2,2));
    
    cost = orbital_period .* delta_v_mag;

end