function OrbitGenerator(r1_mat,r2_mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Orbit Generator
%   Author: Chloe Johnson
%   Date: 4/2/2026
%      
%   Inputs:
%   r1_vecs: Matrix of ephemerous positions (X,3)
%   r2_vecs: Matrix of final positions (X,3)
%
%   Generates visuals for given orbits
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization

orbit_matrix = zeros(size(r1_mat,1),6); % 6 is the length of the keplerian state vector
cost_vec = zeros(size(r1_mat,1)); % Cost Function Vector
deltav_vec = zeros(size(r1_mat,1)); % DeltaV Vector (km/s)
%% Orbit Matrix Generation
for i = 1:size(orbit_matrix,1)
    [orb,cost,deltav] = orbit_solver(r1_mat(i,:),r2_mat(i,:)); % Orbit Solver
    cost_vec(i) = cost; % Individual Cost Vector Assignment
    orbit_matrix(i,:) = orb; % Individual Orbital Matrix Assignment
    deltav_vec(i) = deltav; % Individual Deltav Assignment
end
sma_vec = orbit_matrix(1,:); % SMA Assignment (km)
ecc_vec = orbit_matrix(2,:); % Eccentricity Assignment


%% Visualization

figure(1) % Eccentricity Plot
plot(ecc_vec,sma_vec,"bx")
xlabel("Eccentricity")
ylabel("Semi Major Axis (km)")
grid on;
title("Eccentricity vs SMA")

figure(2) % Cost Plot
plot(cost_vec,sma_vec,"bx")
xlabel("Cost")
ylabel("Semi Major Axis (km)")
grid on;
title("Cost vs SMA")

figure(3) % Delta V Plot
plot(deltav_vec,sma_vec,"bx")
xlabel("Delta V (km/s)")
ylabel("Semi Major Axis (km)")
grid on;
title("Delta V vs SMA")

end