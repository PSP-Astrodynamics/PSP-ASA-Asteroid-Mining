clc
clear
close all

% Define the filename
filename = 'hektor_results.txt';

% Read the file and extract data
data = importdata(filename);

% Extract columns from the data
time = data(:, 1); % Julian Date (not used for plotting, but kept for reference)
r_H_vec = [data(:, 2),data(:, 3),data(:, 4)];

%% earth
r_E_vec = planetEphemeris(time,'Sun','Earth');

%% lamber
a_E = 149597898;
mu_S = 132712440017.99; % Gravitational Parameter of Sun
% mu_E = 398600.4415;     % Gravitational Parameter of Earth
N = 10;
% TOF_hoh = 8.6176e+07;
% period_E = 31558205; % sec

% r1_vec = [-71901356.638820,      -7888916.570738,     129962325.986539]; %Intial position vector of sataliete km
% r2_vec = [-22048919.391468,     -42955627.467596,     779979918.968054]; %Sataliet final position km
r1_vec = r_E_vec(1,:); % Initial Position of the Earth
r2_vec = r_H_vec(1,:); % Inital Position of Hektor

% unit vectors
r1 = norm(r1_vec);
r2 = norm(r2_vec);

c = norm(r2_vec - r1_vec);

s = 1/2 * (r1 + r2 + c);

IP_E = 365.25*24*3600; % s
IP_trans = N*IP_E;

a_trans = (mu_S*(IP_trans/(2*pi))^2)^(1/3);

alpha0 = 2*asin(sqrt(s/(2*a_trans)));

beta0 = 2*asin(sqrt((s-c)/(2*a_trans)));

%TOF = (1/sqrt(mu_S))* a_trans^(3/2) * (alpha0 - sin(alpha0) - (beta0 - sin(beta0)));
lambert_solutions = lambertSolver(a_trans, c, s, mu_S);

% pick 1A
%a_trans = lambert_solutions{2,1};

%p_trans = 4*a_trans*(s-r1)*(s-r2)*(sin(0.5*(alpha0+beta0))^2)/c^2;
p_trans = 4*a_trans*(s-r1)*(s-r2)*(sin(0.5*(alpha0-beta0))^2)/c^2;
e_trans = sqrt(1-p_trans/a_trans);

TA = acos(dot(r1_vec,r2_vec)/(r1*r2));

f = 1-(r2/p_trans)*(1-cos(TA));
g = r1*r2*sin(TA)/sqrt(mu_S*p_trans);

v1_vec = (r2_vec-f*r1_vec)/g;

theta_star1_plus = -acos((p_trans/r1-1)/e_trans);
omega_trans = -(theta_star1_plus-0);

h_hat = cross(r1_vec,r2_vec)/(r1*r2);
i_trans = acos(dot(h_hat,[0,0,1]));

r1_hat = r1_vec/r1;

RAAN_trans = asin(h_hat(1)/sin(i_trans));

figure
plotOrbit3(RAAN_trans, i_trans, omega_trans, p_trans, e_trans, linspace(0,2*pi,1000), 'g', 1, 1, [0,0,0],0,1.5)

%% plot
% MATLAB script to plot Hektor's orbit from ephemeris file
% stk.v.12.2
% BEGIN Ephemeris
%     InterpolationMethod     Lagrange
%     InterpolationOrder      5
%     DistanceUnit            Kilometers
%     CentralBody             Sun
%     CoordinateSystem        ICRF
%     TimeFormat JDate
%     EphemerisTimePosVel

% MATLAB script to plot Hektor's orbit from the ephemeris file

% Create a 3D plot of the orbit
hold on;
plot3(r_H_vec(:,1),r_H_vec(:,2),r_H_vec(:,3), 'LineWidth', 1.5);
hold on;
plot3(r_E_vec(:,1),r_E_vec(:,2),r_E_vec(:,3), 'LineWidth', 1.5)

% Plot the Sun at the origin
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');

% Add labels and title
xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');
title('3D Orbit of Hektor Around the Sun');
grid on;
axis equal;

plot3([0,r2_vec(1)],[0,r2_vec(2)],[0,r2_vec(3)])
plot3([0,r1_vec(1)],[0,r1_vec(2)],[0,r1_vec(3)])

% Add a legend
legend('Lambert Orbit','Hektor Orbit', 'Earth','Sun', 'Hektor intial Position vector', 'Earth intial position vector');

%% function

function TOF_solutions = lambertSolver(a, c, s, mu)
    % Define alpha0 and beta0
    alpha0 = 2 * asin(sqrt(s ./ (2 * a)));  % radians
    beta0 = 2 * asin(sqrt((s - c) ./ (2 * a)));  % radians

    % Define the four equations
    TOF1A = a^(3/2) * (alpha0 - sin(alpha0) - (beta0 - sin(beta0))) / sqrt(mu);
    TOF1B = a^(3/2) * ((2 * pi - (alpha0) - sin(alpha0)) - (beta0 - sin(beta0))) / sqrt(mu);
    TOF2A = a^(3/2) * (alpha0 - sin(alpha0) + (beta0) - sin(beta0)) / sqrt(mu);
    TOF2B = a^(3/2) * ((2 * pi - (alpha0) - sin(alpha0)) + ((beta0) - sin(beta0))) / sqrt(mu);
    
    % answer key: [1A; 1B; 2A; 2B]
    TOF_solutions = {'1A', '1B', '2A', '2B'; TOF1A, TOF1B, TOF2A, TOF2B};
    
end

function [] = plotOrbit3(RAAN, inc, omega, p, e, theta_star, color, scale, grade, c,arrow,W)
    
    r_vec_xyz = zeros(3,length(theta_star));
    for n=1:length(theta_star)
        theta = theta_star(n) + omega;

        r_vec_rth = p/(1+e*cos(theta - omega)) * [1, 0, 0];

        ICR = [cos(RAAN)*cos(theta) - sin(RAAN)*cos(inc)*sin(theta), -cos(RAAN)*sin(theta)...
        - sin(RAAN)*cos(inc)*cos(theta), sin(RAAN)*sin(inc);
           sin(RAAN)*cos(theta) + cos(RAAN)*cos(inc)*sin(theta),...
           -sin(RAAN)*sin(theta) + cos(RAAN)*cos(inc)*cos(theta), -cos(RAAN)*sin(inc);
           sin(inc)*sin(theta), sin(inc)*cos(theta), cos(inc)];
    
        r_vec_xyz(:,n) = (ICR*r_vec_rth');
    end

    x = r_vec_xyz(1,:) + c(1);
    y = r_vec_xyz(2,:) + c(2);
    z = r_vec_xyz(3,:) + c(3);

    plot3(x,y,z, color, LineWidth=W)
    hold on
    if (arrow == 1)
        plotOrbitWithArrows(x, y, z, length(x)/10, color, scale, grade)
    end
    hold on
    
end

function plotOrbitWithArrows(x, y, z, n, color, scale, grade)
    
    % Calculate velocity components (derivatives of position)
    vx = grade*gradient(x);  % Velocity in x direction (approximate derivative)
    vy = grade*gradient(y);  % Velocity in y direction (approximate derivative)
    vz = grade*gradient(z);
    
    hold on;
    
    % Add arrowheads every n points
    idx = 1:n:length(x);  % Select points every n points for arrow placement
    
    % Plot arrowheads only (no body)
    quiver3(x(idx), y(idx), z(idx), vx(idx), vy(idx), vz(idx), scale, color(1), 'MaxHeadSize', 1, 'AutoScale', 'off');  % Set 0 for arrow body size
    
end