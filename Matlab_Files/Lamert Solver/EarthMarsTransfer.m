%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA
% Single Lambert Transfer Demo
% Author: Travis Hastreiter 
% Created On: 24 January, 2026
% Description: Earth to Mars two burn impulsive transfer using ivLam2
% Lambert solver
% NOTE: dV mentioned in this script is actually the relative velocity with 
% respect to the Earth at departure + the relative velocity with respect to
% Mars at arrival. To calculate actual change in velocity you must specify
% the Earth departure orbit and Mars arrival orbit and use patched conics.
% Most Recent Change: 25 January, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Orbit Information (from ssd.jpl.nasa.gov/planets/approx_pos.html)
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

% Mars orbit
% Mars      1.52371034      0.09339410      1.84969142       -4.55343205    -23.94362959     49.55953891
a_M = 1.52371034; % [AU]
e_M = 0.09339410;
i_M = deg2rad(1.84969142); % [rad]
L_M = deg2rad(-4.55343205); % [rad] M = L - longitude_perihelion
longitude_perihelion_M = deg2rad(-23.94362959); % [rad] omega = longitude_perihelion - Omega
Omega_M = deg2rad(49.55953891); % [rad]
omega_M = longitude_perihelion_M - Omega_M; % [rad] argument of periapsis
M_M = L_M - longitude_perihelion_M; % [rad] Mean anomaly
xepoch_keplerian_M = [a_M; e_M; i_M; Omega_M; omega_M; M_M];

%% Initialization
% Inputs
t0_yr = 9.8; % [yr] time until start of transfer
ToF_yr = 1; % [yr] transfer time of flight

% Physical constants
AU = 149597898; % [km]
mu_sun = 132712440017.99; % [km3 / s2]

% Set up dimensionalization constants
char_star.mu = mu_sun; % [km3 / s2]
char_star.l = AU; % [km]
char_star.t = sqrt(char_star.l ^ 3 / char_star.mu); %[s]
char_star.v = char_star.l / char_star.t; %[km / s]

%% Nondimensionalize and Convert Elements
% Nondimensionalize
t0 = year_to_sec(t0_yr) / char_star.t;
ToF = year_to_sec(ToF_yr) / char_star.t;

% Get states at start of transfer (after t0)
x0_keplerian_E = xepoch_keplerian_E;
x0_keplerian_E(6) = x0_keplerian_E(6) + sqrt(1 / x0_keplerian_E(1) ^ 3) * t0;
xf_keplerian_M = xepoch_keplerian_M;
xf_keplerian_M(6) = xf_keplerian_M(6) + sqrt(1 / xf_keplerian_M(1) ^ 3) * (t0 + ToF);

% Convert keplerian state to cartesian
x0_cartesian_E = keplerian_to_cartesian(x0_keplerian_E, [], 1) %first 3 rows position coords, next 3 are velocity vector coords
xf_cartesian_M = keplerian_to_cartesian(xf_keplerian_M, [], 1) %first 3 rows position coords, next 3 are velocity vector coords

%% Lambert Solve
% Load solver
load_lambert();

% Solve - look at lowest dV direction
[vel1_pos, vel2_pos, dV_pos] = best_lambert_zeroN(x0_cartesian_E, xf_cartesian_M, ToF, 0, 0, direction = 1);
[vel1_neg, vel2_neg, dV_neg] = best_lambert_zeroN(x0_cartesian_E, xf_cartesian_M, ToF, 0, 0, direction = -1);

if dV_pos > dV_neg % Choose negative direction
    vel1 = vel1_neg;
    vel2 = vel2_neg;
    dV = dV_neg;
else % Choose positive direction
    vel1 = vel1_pos;
    vel2 = vel2_pos;
    dV = dV_pos;
end

% Unload solver
unload_lambert();

%% Calculate Transfer Orbit Properties
% Calculate Keplerian 
x0_cartesian_transfer = [x0_cartesian_E(1:3); vel1];
xf_cartesian_transfer = [xf_cartesian_M(1:3); vel2];

[x0_keplerian_transfer, thetastar_0_transfer] = cartesian_to_keplerian(x0_cartesian_transfer, [0; 0; 1], [1; 0; 0], 1);
[xf_keplerian_transfer, thetastar_f_transfer] = cartesian_to_keplerian(xf_cartesian_transfer, [0; 0; 1], [1; 0; 0], 1);
% All of the keplerian elements should be the same except the mean anomaly (last one)

% Fix the true anomalies' quadrants for plotting
if thetastar_0_transfer > thetastar_f_transfer
    thetastar_f_transfer = thetastar_f_transfer + 2 * pi;
end

%% Plot Results
p_E = a_E * (1 - e_E ^ 2);
p_M = a_M * (1 - e_M ^ 2);
p_transfer = x0_keplerian_transfer(1) * (1 - x0_keplerian_transfer(2) ^ 2)

figure
plotOrbit3(Omega_E, i_E, omega_E, p_E, e_E, linspace(0, 2 * pi, 200), "b", 0.5, 1, [0, 0, 0], 1, 1);
plotOrbit3(Omega_M, i_M, omega_M, p_M, e_M, linspace(0, 2 * pi, 200), "r", 0.5, 1, [0, 0, 0], 1, 1);
plotOrbit3(x0_keplerian_transfer(4), x0_keplerian_transfer(3), x0_keplerian_transfer(5), p_transfer, x0_keplerian_transfer(2), linspace(thetastar_0_transfer, thetastar_f_transfer, 200), "g", 0.5, 1, [0, 0, 0], 1, 1);
axis equal
grid on
legend("Earth", "", "Mars", "", "Transfer", "")
xlabel("X [AU]")
ylabel("Y [AU]")
zlabel("Z [AU]")
title(sprintf("Earth to Mars Lambert Transfer in %.1f Years", ToF_yr))
subtitle(sprintf("Delta V: %.2f km / s", dV * char_star.v))