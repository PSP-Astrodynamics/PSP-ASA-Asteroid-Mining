clc
clear
close all

% Plot formatting
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

% Define the filename
filename = 'hektor_results.txt';

% Read the file and extract data
targetEphemeris = importdata(filename);

% Extract columns from the data
time = targetEphemeris(:, 1); % Julian Date (not used for plotting, but kept for reference)

%time, x, y, z
r_target_vec = [targetEphemeris(:, 2),targetEphemeris(:, 3),targetEphemeris(:, 4)];
[minZ, idx_AN] = min(abs(r_target_vec(:,3)));

%% earth
r_E_vec = planetEphemeris(time,'Sun','Earth');

%% lambert
mu_e = 3.986e5; % Gravitational Parameter of Earth
mu_S = 132712440017.99; % Gravitational Parameter of Sun
r_earth = 149.9e6; % radius of earth orbit [km]

r1_vec = r_E_vec(idx_AN,:); % Initial Position of the Earth
r2_vec = r_target_vec(idx_AN,:); % Target Position of Hektor

[a_trans, c, s, r1, r2] = lambertIntial(r1_vec,r2_vec,mu_S);

% intialize matrix to zeros
e_solutions = zeros(length(a_trans), 4);

% create 100x2 matrix of eccentricity solutions for lambert solutions 1 & 2
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

    % this is not good and assumes that the orbit around earth is circular
    % and that they intersect at the end of r1 + radius above earth, most likly will break on
    % full 3D earth has a circular orbit
    % v_orbit = sqrt(mu_S*(2/r1 - 1/a_trans(index)));
    % h_orbit = sqrt(mu_S * [pAB_1, pAB_2, pAB_1other, pAB_2other]);
    % v_earth = sqrt(mu_S/r_earth);
    % h_earth = sqrt(mu_S * r_earth);
    % gamma_cycler = acos(h_orbit ./ (r1 .* v_orbit));
    % gamma_earth_scalar = acos(h_earth / (r_earth * v_earth));
    % gamma_earth_vector = [gamma_earth_scalar, gamma_earth_scalar, gamma_earth_scalar, gamma_earth_scalar];
    % gamma_inf = gamma_cycler-gamma_earth_vector;
    % v_inf(index,:) = sqrt(v_earth^2 + v_orbit^2 - 2 * v_orbit * v_earth * cos(gamma_inf));
end

RA = {};
RB = {};

for index = 1:100
    [x1,y1,z1] = plotCycler(e_solutions(index,1), r1_vec, r2_vec, r1, p_solutions(index,1), '#4DBEEE');
    [x2,y2,z2] = plotCycler(e_solutions(index,2), r1_vec, r2_vec, r1, p_solutions(index,2), '#1b10c2');
    RA{index} = [x1;y1;z1];
    RB{index} = [x2;y2;z2];
end

VA = {};
VB = {};

for index = 1:100
    R1A = RA{index};
    R1B = RB{index};
    VA{index} = (R1A(:,2) - R1A(:,1))./(365.25*24*3600*index/size(R1A,2));
    VB{index} = (R1B(:,2) - R1B(:,1))./(365.25*24*3600*index/size(R1B,2));
end

VE = (r_E_vec(idx_AN+1,:) - r_E_vec(idx_AN,:))./((time(idx_AN+1)-time(idx_AN))*86400);

delVAmag = {};
delVBmag = {};
delVA = {};
delVB = {};

for index = 1:100
    VAtrans = VA{index}';
    VBtrans = VB{index}';
    delVAtrans = VAtrans - VE;
    delVBtrans = VBtrans - VE;
    delVA{index} = delVAtrans;
    delVB{index} = delVBtrans;
    delVAmag{index} = norm(delVAtrans);
    delVBmag{index} = norm(delVBtrans);

end

% plot eccentricity solutions vs semi major axes
figure(1)
hold on
plot(a_trans, e_solutions(:, 1), "o")
plot(a_trans, e_solutions(:, 2), "^")
plot(a_trans, e_solutions(:, 3), 'pentagram')
plot(a_trans, e_solutions(:, 4), 'v')
legend("Solution 1", "Solution 2", location = "northwest", FontSize = 10)
title("Eccentricity vs. Semi-Major Axis", 'FontSize', 14, 'FontWeight', 'bold')
xlabel("Semi-Major Axis (AU)", 'FontSize', 12, 'FontWeight', 'bold')
ylabel("Eccentricity", 'FontSize', 12, 'FontWeight', 'bold')
grid on
grid minor
hold off
change_axis_from_Km_to_AU(true, false);

% plotting the orbits of the body of intrest and Earth
figure(2)
hold on;
plot3(r_target_vec(:,1),r_target_vec(:,2),r_target_vec(:,3), 'color', '#D95319','LineWidth', 1.5);
hold on;
plot3(r_E_vec(:,1),r_E_vec(:,2),r_E_vec(:,3), 'color', '#EDAA1A', 'LineWidth', 1.5)

% Plot the Sun at the origin
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', '#F2E01D');

plot3([0,r1_vec(1)],[0,r1_vec(2)],[0,r1_vec(3)],'Color', '#D93319', 'LineWidth', 1.5)
plot3([0,r2_vec(1)],[0,r2_vec(2)],[0,r2_vec(3)], 'Color', '#7E2F8E','LineWidth', 1.5)

% Create a 3D plot of the orbits
hold on
for index = 6:16
    [x1,y1,z1] = plotCycler(e_solutions(index,1), r1_vec, r2_vec, r1, p_solutions(index,1), '#4DBEEE');
    [x2,y2,z2] = plotCycler(e_solutions(index,2), r1_vec, r2_vec, r1, p_solutions(index,2), '#1b10c2');
    plot3(x1, y1, z1, 'Color', '#4DBEEE', 'LineWidth', 1.5);
    plot3(x2, y2, z2, 'Color', '#1b10c2', 'LineWidth', 1.5);
end

% Add labels and title
xlabel('X Position (AU)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Position (AU)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Z Position (AU)', 'FontSize', 12, 'FontWeight', 'bold');
title('A Selection Synodic Cycler Orbits', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
grid minor;
axis equal;
change_axis_from_Km_to_AU(true, true);

% Add a legend
legend('Hektor Orbit', 'Earth','Sun', 'Intial Position vector', 'Target Position vector', 'Location','southwest', FontSize = 10);

% plotting tof vs sma
lambertSolverTOFvsSMA(1495978707/10,mu_S,c,s)

% plotting delVmag vs. sma

