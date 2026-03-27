function [delta_v, orbit_matrix] = lambert_solver(r1_vec,r2_vec)
mu_S = 1.989e30 * 6.674e-11; %in SI units
n = 100; %number of revolutions for lambertInitial Solver, hardcoded 

%coding chain of funciton to generate state vectors and delta v
[a_trans, c, s, r1, r2] = lambertInitial(r1_vec,r2_vec,mu_S,n); 

p_solution = lambertSolverP(a_trans, c, s, r1, r2);

TOF_solutions = lambertSolverTOF(a_trans, c, s, mu_S);

%TOF vs SMA plots
R = 6378.14; % Defined as radius of Earth, a_range will be up to 18 earth radii
lambertSolverTOFvsSMA(R,mu,c,s)

%solving for delta V
r1 = norm(r1_vec); %not sure why r1 is defined as such in all the other functions
= deltaV(r1_vec, r2_vec, r1, e_solutions, p_solutions, VE)

