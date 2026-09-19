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
t0_yr = linspace(0, 20, 500); % [yr] time until start of transfer
ToF_yr = linspace(0.1, 8, 500); % [yr] transfer time of flight

% Physical constants
AU = 149597898; % [km]
mu_sun = 132712440017.99; % [km3 / s2]

% Set up dimensionalization constants
char_star.mu = mu_sun; % [km3 / s2]
char_star.l = AU; % [km]
char_star.t = sqrt(char_star.l ^ 3 / char_star.mu); %[s]
char_star.v = char_star.l / char_star.t; %[km / s]

%% Prepare Initial and Final Position and Velocity Arrays
% Nondimensionalize
t0 = year_to_sec(t0_yr) / char_star.t;
ToF = year_to_sec(ToF_yr) / char_star.t;

% Get all combinations of initial time and time of flight that must be solved for
t0_ToF_combs = combinations(t0, ToF);
t0_combs = t0_ToF_combs.t0';
ToF_combs = t0_ToF_combs.ToF';

Q = height(t0_ToF_combs); % Number of Lambert solves

%%
% Get states at start of transfer (after t0)
x0_keplerian_E = xepoch_keplerian_E .* ones([6, Q]);
x0_keplerian_E(6, :) = x0_keplerian_E(6, :) + sqrt(1 ./ x0_keplerian_E(1, :) .^ 3) .* t0_combs;
xf_keplerian_A = xepoch_keplerian_A .* ones([6, Q]);
xf_keplerian_A(6, :) = xf_keplerian_A(6, :) + sqrt(1 ./ xf_keplerian_A(1, :) .^ 3) .* (t0_combs + ToF_combs);

% Convert keplerian state to cartesian
x0_cartesian_E = keplerian_to_cartesian_array(x0_keplerian_E, [], 1);
xf_cartesian_A = keplerian_to_cartesian_array(xf_keplerian_A, [], 1);

%% Lambert Solve
% Load solver
load_lambert();

% Solve - look at lowest dV direction
vel1_posneg = zeros([3, Q, 2]);
vel2_posneg = zeros([3, Q, 2]);
dV_posneg = zeros([Q, 2]);
[vel1_posneg(:, :, 1), vel2_posneg(:, :, 1), dV_posneg(:, 1)] = best_lambert_zeroN(x0_cartesian_E, xf_cartesian_A, ToF_combs, 0, 0, direction = ones([1, Q]));
[vel1_posneg(:, :, 2), vel2_posneg(:, :, 2), dV_posneg(:, 2)] = best_lambert_zeroN(x0_cartesian_E, xf_cartesian_A, ToF_combs, 0, 0, direction = -ones([1, Q]));

[dV, dV_dir_i] = min(dV_posneg, [], 2);
vel1 = zeros([3, Q]);
vel2 = zeros([3, Q]);
for i = 1 : Q
    vel1(:, i) = vel1_posneg(:, i, dV_dir_i(i));
    vel2(:, i) = vel2_posneg(:, i, dV_dir_i(i));
end

% Unload solver
unload_lambert();

%% Calculate Transfer Orbit Properties
% Calculate Keplerian 
x0_cartesian_transfer = [x0_cartesian_E(1:3, :); vel1];
xf_cartesian_transfer = [xf_cartesian_A(1:3, :); vel2];

[x0_keplerian_transfer, thetastar_0_transfer] = cartesian_to_keplerian_array(x0_cartesian_transfer, [0; 0; 1], [1; 0; 0], 1);
[xf_keplerian_transfer, thetastar_f_transfer] = cartesian_to_keplerian_array(xf_cartesian_transfer, [0; 0; 1], [1; 0; 0], 1);
% All of the keplerian elements should be the same except the mean anomaly (last one)

%% Cutoff dV
dV_cutoff = 20; % [km / s]
dV_paired = dV;
dV_paired(dV > dV_cutoff / char_star.v) = nan;

dV_dir_i_paired = dV_dir_i;
dV_dir_i_paired(dV > dV_cutoff / char_star.v) = nan;

%% Plot Results
X = reshape(sec_to_year(t0_combs * char_star.t), [numel(t0), numel(ToF)]);
Y = reshape(sec_to_year(ToF_combs * char_star.t), [numel(t0), numel(ToF)]);
Z_dV = reshape(dV_paired * char_star.v, [numel(t0), numel(ToF)]);
Z_dir = (1.5 - reshape(dV_dir_i_paired, [numel(t0), numel(ToF)])) * 2;

figure
pcolor(X, Y, Z_dir, FaceColor="flat",EdgeAlpha=0,HandleVisibility="off");
colormap("jet")
colorbar()
xlabel("Initial Time [yr]")
ylabel("Time of Flight [yr]")
title("Best Direction vs t0 and ToF")
grid on

figure
pcolor(X, Y, Z_dV, FaceColor="flat",EdgeAlpha=0,HandleVisibility="off");
colormap("jet")
colorbar()
xlabel("Initial Time [yr]")
ylabel("Time of Flight [yr]")
title("Earth to Psyche Porkchop Plot")
subtitle("Total Relative Velocity vs t0 and ToF")
grid on
