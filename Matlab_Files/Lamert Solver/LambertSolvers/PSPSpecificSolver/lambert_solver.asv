function orbit_matrix = lambert_solver(r1_vec,r2_vec)

char_values = load_characteristic_values();
mu_S  = char_values.mu;

n = 100; %number of revolutions for lambertInitial Solver, hardcoded 

%coding chain of funciton to generate state vectors and delta v
[a_trans, c, s, r1, r2] = lambertInitial(r1_vec,r2_vec,mu_S,n); 

% Solving for P and E solutions
for index = 1:100
    p_sol = lambertSolverP(a_trans(index), c, s, r1, r2);
    pAB_1 = p_sol{2, 1};
    pAB_2 = p_sol{2, 2};
    pAB_1other = p_sol{2, 3};
    pAB_2other = p_sol{2, 4};
    p_solutions(index,:) = [pAB_1, pAB_2, pAB_1other, pAB_2other];
    e_AB_1 = sqrt(1 - (pAB_1 / a_trans(index)));
    e_AB_2 = sqrt(1 - (pAB_2 / a_trans(index)));
    e_AB_1other = sqrt(1 - (pAB_1other / a_trans(index)));
    e_AB_2other = sqrt(1 - (pAB_2other / a_trans(index)));
    e_solutions(index, :) = [e_AB_1, e_AB_2, e_AB_1other, e_AB_2other];
end

%constructing the state vector with orbital elements
for index = 1:length(a_trans)
    orbit_matrix(index,1) = a_trans(index); %semimajor axis
    orbit_matrix(index,2) = e_solutions(index,1); %eccentricity 
    [Omega, theta, inclination] = orbitparameters(r1_vec,r2_vec); %orbital elements
    orbit_matrix(index,3) = inclination; 
    orbit_matrix(index,4) = Omega; %this Omega is taken as RAAN
theta_star = acos(((p_solutions(index,1) / r1) - 1) / e_solutions(index,1)); %theta star is true anomaly, formula taken from PlotCycler
    orbit_matrix(index,5) = theta - theta_star;
    %longitude of periapsis to true anomaly
    orbit_matrix(index,6) = theta_star;
end

for index = 1+length(a_trans):2*length(a_trans)
    orbit_matrix(index,1) = a_trans(index); %semimajor axis
    orbit_matrix(index,2) = e_solutions(index,2); %eccentricity 
    [Omega, theta, inclination] = orbitparameters(r1_vec,r2_vec); %orbital elements
    orbit_matrix(index,3) = inclination; 
    orbit_matrix(index,4) = Omega; %this Omega is taken as RAAN
    theta_star = acos(((p_solutions(index,1) / r1) - 1) / e_solutions(index,1)); %theta star is true anomaly, formula taken from PlotCycler
    orbit_matrix(index,5) = theta - theta_star;
    %longitude of periapsis to true anomaly
    orbit_matrix(index,6) = theta_star;
end