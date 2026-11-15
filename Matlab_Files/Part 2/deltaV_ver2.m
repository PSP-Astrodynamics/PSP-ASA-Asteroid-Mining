%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Description: This function calculates the impulsive delta V
% required for a maneuver between escaping Earth's escape velocity and
% being inserted into the Cycler orbit.
% 
%   INPUT VARIABLES:
%   NAME : Description. (unit) 
%   r1_vec : Vector of the Earth's position with respect to the Sun. (km) 
%   otbit_matrix : State vector containing the state of Lambert Solutions
%   from the lambert Solver
%  
%   OUTPUT VARIABLES:
%   NAME : Description. (unit) 
%   delVAmag: Magnitude of delta V for A solutions
%   delVA: implusive delta V vector in heliocentric inertial coordinates
%   delVBmag: Magnitude of delta V for B solutions
%   delVB: implusive delta V vector in heliocentric inertial coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_perifocal] = deltaV_ver2(r1_vec,orbit_matrix)
%Original function outputs: [delVAmag, delVBmag, delVA, delVB]
%Loading characteristic values
char_values = load_characteristic_values();

% Unpacking the orbit matrix for required elements
theta_star = orbit_matrix(:,6); % TrueAnomaly
p = orbit_matrix(:,7); %semilatus rectum
e = orbit_matrix(:,2);

v_p = []; %empty array for 'p' component of velocity
v_e = []; %empty array for 'e' component of velocity
v_z = []; %empty array for 'z' component of velocity

%creating a vector of velocities in the cycler orbit perifocal frame 
for solutionIndex = 1:size(orbit_matrix,1)
    v_p(solutionIndex) = -(sqrt(char_values.mu ./ p(solutionIndex))) .* sin(theta_star(solutionIndex));
    v_e(solutionIndex) = (sqrt(char_values.mu ./ p(solutionIndex))) .* (e(solutionIndex) + sin(theta_star(solutionIndex)));
    v_z(solutionIndex) = 0; %velocity component in the plane normal to the elliptical orbit
    v_perifocal = [v_p ;v_e ;v_z]'; %CHECKPOINT, SENSE-MAKE THESE VALUES WITH HAND CALCULATED VALUES
end

%converting perifocal velocities to heliocentric inertial cartesian frame





end