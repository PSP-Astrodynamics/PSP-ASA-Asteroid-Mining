clc
clear
close all

mu_S = 132712440017.99; % km^3/s^2
AU = 149597898; % km
R_S = 695990; % km

r1_vec_xyz = 1.83*AU*[cosd(33.9)*cosd(43.7), cosd(33.9)*sind(43.7), sind(33.9)];
r1 = norm(r1_vec_xyz);

r2_vec_xyz = 3.15*AU*[cosd(48.2)*cosd(62.4), cosd(48.2)*sind(62.4), sind(48.2)];
r2 = norm(r2_vec_xyz);

%% part A
% type 1: TA < 180 deg
TOF = (3.26/12)*365.25*24*3600; % s

% type 1
TA = acos(dot(r1_vec_xyz,r2_vec_xyz)/(r1*r2)); % rad

c_vec_xyz = r2_vec_xyz-r1_vec_xyz;
c = norm(c_vec_xyz);

s = 0.5*(r1+r2+c);

% parabolic TOF check 
TOF1 = (1/3)*sqrt(2/mu_S)*(s^(3/2)-(s-c)^(3/2));
TOF2 = (1/3)*sqrt(2/mu_S)*(s^(3/2)+(s-c)^(3/2));

% TOF1 and TOF2 parabollic are both > TOF, so this is hyperbolic

% must be 1H

% answer key: [1A; 1B; 2A; 2B; 1H; 2H]
a_sol = lambertSolverSMA(TOF,c,s,mu_S);

a_1H = abs(a_sol{2,5});

a_min = s/2;
TOF_min_sol = lambertSolverTOF(a_min,c,s,mu_S);
% we know its hyperbolic
a_trans = a_1H;

p_sol = lambertSolverP(a_trans, c, s, r1, r2);
% for 1H it must be solution with higher eccentricity which is highest p
% because: e=sqrt(p/abs(a)+1) in a hyperbola
p_trans = p_sol{2,3};

e_trans = sqrt(p_trans/a_trans+1);

h_trans_hat_xyz = cross(r1_vec_xyz,r2_vec_xyz)/norm(cross(r1_vec_xyz,r2_vec_xyz));

i_trans = acos(dot(h_trans_hat_xyz,[0,0,1]));

% ascending
theta_starD_plus = acos((p_trans/r1-1)/e_trans);
theta_starA_minus = theta_starD_plus+TA;

r1_hat_xyz = r1_vec_xyz/r1;
theta_hatD_xyz = cross(h_trans_hat_xyz,r1_hat_xyz);

thetaD_plus = pi-asin(r1_hat_xyz(3)/sin(i_trans));
thetaD_plus = 2*pi-acos(theta_hatD_xyz(3)/sin(i_trans));
thetaD_plus = asin(r1_hat_xyz(3)/sin(i_trans));
thetaD_plus = acos(theta_hatD_xyz(3)/sin(i_trans));

omega_trans = thetaD_plus-theta_starD_plus;

RAAN_trans = pi-asin(h_trans_hat_xyz(1)/sin(i_trans));
RAAN_trans = 2*pi-acos(-h_trans_hat_xyz(2)/sin(i_trans));
RAAN_trans = asin(h_trans_hat_xyz(1)/sin(i_trans));
RAAN_trans = acos(-h_trans_hat_xyz(2)/sin(i_trans));

figure
% plotOrbit3(RAAN_trans, i_trans, domega, p_trans, e_trans, linspace(theta_starD_plus,theta_starA_minus,1000), 'r', 0.5, 1, [0,0,0],1,1)
plotOrbit3(RAAN_trans, i_trans, omega_trans, p_trans/AU, e_trans, linspace(theta_starD_plus,theta_starA_minus,1000), 'r', 0.5, 1, [0,0,0],1,1)
hold on
plotOrbit3(RAAN_trans, i_trans, omega_trans, p_trans/AU, e_trans, linspace(-theta_starA_minus,theta_starA_minus,1000), 'r--', 0.5, 1, [0,0,0],0,1)
plotOrbit3(RAAN_trans, i_trans, omega_trans, p_trans/AU, e_trans, linspace(0,0,1), 'k*', 0.5, 1, [0,0,0],0,1)
earthy(10000000/AU, 'Sun', 1,[0,0,0])
axis equal
grid on
title('Orbital Diagram')
xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
plot3([0,r1_vec_xyz(1)]/AU, [0,r1_vec_xyz(2)]/AU, [0,r1_vec_xyz(3)]/AU, 'k', LineWidth=1)
plot3([0,r2_vec_xyz(1)]/AU, [0,r2_vec_xyz(2)]/AU, [0,r2_vec_xyz(3)]/AU, 'k', LineWidth=1)
plot3([r1_vec_xyz(1),r2_vec_xyz(1)]/AU, [r1_vec_xyz(2),r2_vec_xyz(2)]/AU, [r1_vec_xyz(3),r2_vec_xyz(3)]/AU, 'y', LineWidth=1)
view(0,0)
xlim([-2.5,1.5])
zlim([-0.75, 3])
legend('Transfer Arc','', 'Transfer Hyperbola', 'Periapsis Radius', 'Location', 'northwest');

H_D = 2*atanh((((e_trans-1)/(e_trans+1))^0.5)*tan(theta_starD_plus/2));

H_A = 2*atanh((((e_trans-1)/(e_trans+1))^0.5)*tan(theta_starA_minus/2));

t_tpD = (e_trans*sinh(H_D) - H_D)/sqrt(mu_S/a_trans^3);

t_tpA = (e_trans*sinh(H_A) - H_A)/sqrt(mu_S/a_trans^3);

tA_tD = t_tpA-t_tpD;

check = tA_tD/TOF;

%IP_trans = 2*pi*sqrt(a_trans^3/mu_S);

%% part c

a_E = 149597898; % km

rp = -a_trans*(1-e_trans);
rp_aE = rp/a_E;

r_AN = p_trans/(1+e_trans*cos(-omega_trans));
r_AN_aE = r_AN/a_E;

hold off
lambertSolverTOFvsSMA(10000000*10,mu_S,c,s)

hold off
figure
% plotOrbit3(RAAN_trans, i_trans, domega, p_trans, e_trans, linspace(theta_starD_plus,theta_starA_minus,1000), 'r', 0.5, 1, [0,0,0],1,1)
plotOrbit3(RAAN_trans, i_trans, omega_trans, p_trans/AU, e_trans, linspace(theta_starD_plus,theta_starA_minus,1000), 'r', 0.5, 1, [0,0,0],1,1)
hold on
plotOrbit3(RAAN_trans, i_trans, omega_trans, p_trans/AU, e_trans, linspace(-theta_starA_minus,theta_starA_minus,1000), 'r--', 0.5, 1, [0,0,0],0,1)
plotOrbit3(RAAN_trans, i_trans, omega_trans, p_trans/AU, e_trans, linspace(0,0,1), 'k*', 0.5, 1, [0,0,0],0,1)
plotOrbit3(0, 0, 0, a_E/AU, 0, linspace(0,2*pi,1000), 'g', 0.5, 1, [0,0,0],0,1)
earthy(10000000/AU, 'Sun', 1,[0,0,0])
axis equal
grid on
title('Orbital Diagram')
xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
plot3([0,r1_vec_xyz(1)]/AU, [0,r1_vec_xyz(2)]/AU, [0,r1_vec_xyz(3)]/AU, 'k', LineWidth=1)
plot3([0,r2_vec_xyz(1)]/AU, [0,r2_vec_xyz(2)]/AU, [0,r2_vec_xyz(3)]/AU, 'k', LineWidth=1)
plot3([r1_vec_xyz(1),r2_vec_xyz(1)]/AU, [r1_vec_xyz(2),r2_vec_xyz(2)]/AU, [r1_vec_xyz(3),r2_vec_xyz(3)]/AU, 'y', LineWidth=1)
view(0,0)
%xlim([-2.5,1.5])
%zlim([-0.75, 3])
legend('Transfer Arc','', 'Transfer Hyperbola', 'Periapsis Radius','Earth Orbit', 'Location', 'northwest');





