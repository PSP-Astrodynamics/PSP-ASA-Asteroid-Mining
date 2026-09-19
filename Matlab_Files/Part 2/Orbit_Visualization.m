%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Orbit Information (+)
% Earth orbit
% EM Bary   1.00000261      0.01671123     -0.00001531      100.46457166    102.93768193      0.0
a_E = 1.00000261; % [AU]
e_E = 0.01671123;
i_E = deg2rad(-0.00001531); % [rad]
L_E = deg2rad(100.46457166); % [rad] M = L - longitude_perihelion
longitude_perihelion_E = deg2rad(102.93768193); % [rad] omega = longitude_perihelion - Omega
Omega_E = deg2rad(0); % [rad] Right ascension of the ascending node
omega_E = longitude_perihelion_E - Omega_E; % [rad] argument of periapsis
M_E = L_E - longitude_perihelion_E; % [rad] Mean anomaly
xepoch_keplerian_E = [a_E; e_E; i_E; Omega_E; omega_E; M_E];
r_E = 6378.1; % Earth's radius in km

% Asteroid Belt orbit
% Asteroid  (16 Psyche)    1.52371034      0.09339410      1.84969142       -4.55343205    -23.94362959     49.55953891
a_A = 2.92414501; % [AU]
e_A = 0.13411742;
i_A = deg2rad(3.09684542); % [rad]
L_A = deg2rad(622.51905687); % [rad] M = L - longitude_perihelion
longitude_perihelion_A = deg2rad(379.35905687); % [rad] omega = longitude_perihelion - Omega
Omega_A = deg2rad(150.03); % [rad]
omega_A = longitude_perihelion_A - Omega_A; % [rad] argument of periapsis
M_A = L_A - longitude_perihelion_A; % [rad] Mean anomaly
xepoch_keplerian_A = [a_A; e_A; i_A; Omega_A; omega_A; M_A];

%% Initialization

% Inputs
t0_yr = 8; % [yr] time until start of transfer
ToF_yr = 6; % [yr] transfer time of flight

% Physical constants
AU = 149597898; % [km]
mu_sun = 132712440017.99; % [km3 / s2] (1.32712440017 E20 m3/s2)
mu_E = 398600.4418; % G * mass of Earth

% Set up dimensionalization constants
char_star.mu = mu_sun; % [km3 / s2]
char_star.l = AU; % [km]
char_star.t = sqrt(char_star.l ^ 3 / char_star.mu); %[s]
char_star.v = char_star.l / char_star.t; %[km / s]
%% Search Ranges
t0_yr_vec=linspace(9,11,200);
ToF_yr_vec=linspace(1.0,2.5,200);

dV_min = Inf; 
best_solution = struct();

%% Nondimensionalize and Convert Elements
%Lambert

load_lambert();
for t0_yr = t0_yr_vec
    for ToF_yr = ToF_yr_vec

        % Nondimensionalize
        t0 =(t0_yr) * 365.25 * 24 * 3600 / char_star.t; %converting time into the solver's units
        ToF = (ToF_yr) * 365.25 * 24 * 3600 / char_star.t;

        % Get states at start of transfer (after t0)
        x0_keplerian_E = xepoch_keplerian_E;
        x0_keplerian_E(6) = x0_keplerian_E(6) + sqrt(1 / x0_keplerian_E(1) ^ 3) * t0; %mean anomaly (where Earth is along its orbit)
        xf_keplerian_A = xepoch_keplerian_A;
        xf_keplerian_A(6) = xf_keplerian_A(6) + sqrt(1 / xf_keplerian_A(1) ^ 3) * (t0 + ToF);
    
        % Convert keplerian state to cartesian
        x0_cartesian_E = keplerian_to_cartesian(x0_keplerian_E, [], 1);
        xf_cartesian_A = keplerian_to_cartesian(xf_keplerian_A, [], 1);

    %Lambert Solver (calculates trajectory)
    % Solve - look at lowest dV direction
    try
        [vel1_pos, vel2_pos, dV_pos] = best_lambert_zeroN(x0_cartesian_E, xf_cartesian_A, ToF, 0, 0, direction = 1);    
        [vel1_neg, vel2_neg, dV_neg] = best_lambert_zeroN(x0_cartesian_E, xf_cartesian_A, ToF, 0, 0, direction = -1);
    catch 
        continue
    end
    
    % Pick best direction (whichever has the smaller deltaV)
    if dV_pos > dV_neg % Choose negative direction
        vel1 = vel1_neg;
        vel2 = vel2_neg;
        dV = dV_neg;
    else % Choose positive direction
        vel1 = vel1_pos;
        vel2 = vel2_pos;
        dV = dV_pos;
    end
    %Keep best trajectory (min dV) across all runs
    if dV < dV_min
        dV_min = dV;
        best_solution.t0_yr  = t0_yr;
        best_solution.ToF_yr = ToF_yr;
        best_solution.vel1 = vel1; % transfer velocity at Earth (velocity required to meet the velocity needed)
        best_solution.vel2 = vel2; % transfer velocity at Psyche
        best_solution.x0_E = x0_cartesian_E;
        best_solution.xf_A = xf_cartesian_A;
    end
    
end

end
% Unload solver
unload_lambert();

%print best solution + use it
fprintf("Best t0 = %.3f yr\n", best_solution.t0_yr);
fprintf("Best ToF = %.3f yr\n", best_solution.ToF_yr);
fprintf("Minimum Lambert dV = %.3f km/s\n", dV_min*char_star.v);

%% Leaving Earth Orbit (Burn 1, Escape); mu_E = G * mass of Earth
r_escape_E = r_E + 200; % Earth's radius (km) + 200 km (reasonable assumption)
v_orbit_E = sqrt(mu_E/r_escape_E); % the speed of Earth's circular orbit at 200 km
v_depart_E = norm(best_solution.vel1 - best_solution.x0_E(4:6))* char_star.v; % the speed you need to escape Earth's orbit at 200 km (4:6 grabs velocity vector)
v_escape_E = sqrt(v_depart_E^2 + 2*mu_E/r_escape_E);

dV_escape = v_escape_E - v_orbit_E; % change in velocity (escape velocity - orbital velocity)
%% Entering Psych Orbit (Burn 2, Enter)

r_Psyche = 111; % (mean radius)
mu_Psyche = 1.482; % (gravitational constant of Psyche)

% orbit radii (radius of Psyche + _ km)
r_A = r_Psyche + 700;
r_B = r_Psyche + 290;
r_C = r_Psyche + 170;
r_D = r_Psyche + 85;
% circular orbital velocities
v_A = sqrt(mu_Psyche/r_A);
v_B = sqrt(mu_Psyche/r_B);
v_C = sqrt(mu_Psyche/r_C);
v_D = sqrt(mu_Psyche/r_D);

% semi major axes for transfer orbits
a_AB = (r_A + r_B) / 2;
a_BC = (r_B + r_C) / 2;
a_CD = (r_C + r_D) / 2;


v_arrive = norm(best_solution.vel2 - best_solution.xf_A(4:6)) * char_star.v;
v_hyperbolic_A = sqrt(v_arrive^2 + 2*mu_Psyche/r_A);
dV_arrive = v_hyperbolic_A - v_A;


%% Hohmann Transfer

% A to B
dv_A_AB = abs(sqrt(mu_Psyche * (2/r_A - 1/a_AB)) - v_A); % transfer ellipse - circular orbit
dv_B_AB = abs(v_B - sqrt(mu_Psyche * (2/r_B - 1/a_AB)));
 % go from v_A to elipse speed to V_B

% B to C
dv_B_BC = abs(sqrt(mu_Psyche * (2/r_B - 1/a_BC)) - v_B); % speed at start of transfer ellipse
dv_C_BC = abs(v_C - sqrt(mu_Psyche * (2/r_C - 1/a_BC)));

% C to D
dv_C_CD = abs(sqrt(mu_Psyche * (2/r_C - 1/a_CD)) - v_C);
dv_D_CD = abs(v_D - sqrt(mu_Psyche * (2/r_D - 1/a_CD)));

dv_decending = dv_A_AB + dv_B_AB + dv_B_BC + dv_C_BC + dv_C_CD + dv_D_CD;
%% Total Mission Delta-V Summary
fprintf('\n=== MISSION DELTA-V BUDGET ===\n');
fprintf('Burn 1 - Earth escape:        %.4f km/s\n', dV_escape);
fprintf('Burn 2 - Psyche capture (A):  %.4f km/s\n', dV_arrive);  
fprintf('Burn 3 - A to B decending:     %.4f km/s\n', dv_A_AB + dv_B_AB);
fprintf('Burn 4 - B to C decending:     %.4f km/s\n', dv_B_BC + dv_C_BC);
fprintf('Burn 5 - C to D decending:     %.4f km/s\n', dv_C_CD + dv_D_CD);
fprintf('Lambert interplanetary dV:    %.4f km/s\n', dV_min * char_star.v);
fprintf('==============================\n');
dV_total = dV_escape + dV_arrive + dv_decending + (dV_min * char_star.v);
fprintf('TOTAL MISSION dV: %.4f km/s\n', dV_total);
%% Calculate Transfer Orbit Properties
% Calculate Keplerian 
x0_cartesian_transfer = [best_solution.x0_E(1:3); best_solution.vel1];
xf_cartesian_transfer = [best_solution.xf_A(1:3); best_solution.vel2];

[x0_keplerian_transfer, thetastar_0_transfer] = cartesian_to_keplerian(x0_cartesian_transfer, [0; 0; 1], [1; 0; 0], 1);
[xf_keplerian_transfer, thetastar_f_transfer] = cartesian_to_keplerian(xf_cartesian_transfer, [0; 0; 1], [1; 0; 0], 1);
% All of the keplerian elements should be the same except the mean anomaly (last one)

% Fix the true anomalies' quadrants for plotting
if thetastar_0_transfer > thetastar_f_transfer
    thetastar_f_transfer = thetastar_f_transfer + 2 * pi;
end

%% Plot Results
p_E = a_E * (1 - e_E ^ 2);
p_A = a_A * (1 - e_A ^ 2);
p_transfer = x0_keplerian_transfer(1) * (1 - x0_keplerian_transfer(2) ^ 2);

figure
plotOrbit3(Omega_E, i_E, omega_E, p_E, e_E, linspace(0, 2 * pi, 200), "b", 0.5, 1, [0, 0, 0], 1, 1);
plotOrbit3(Omega_A, i_A, omega_A, p_A, e_A, linspace(0, 2 * pi, 200), "r", 0.5, 1, [0, 0, 0], 1, 1);
plotOrbit3(x0_keplerian_transfer(4), x0_keplerian_transfer(3), x0_keplerian_transfer(5), p_transfer, x0_keplerian_transfer(2), linspace(thetastar_0_transfer, thetastar_f_transfer, 200), "g", 0.5, 1, [0, 0, 0], 1, 1);
axis equal
grid on
legend("Earth", "", "Asteroid Belt", "", "Transfer", "")
xlabel("X [AU]")
ylabel("Y [AU]")
zlabel("Z [AU]")
title(sprintf("Earth to Asteroid Belt Lambert Transfer in: t0= %.1f Years, ToF = %.1f Years", best_solution.t0_yr, best_solution.ToF_yr))
subtitle(sprintf("Delta V: %.2f km / s", dV_min * char_star.v))
fprintf('dV_min = %f\n', dV_min);
disp(best_solution);