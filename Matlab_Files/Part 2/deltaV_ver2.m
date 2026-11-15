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
%   delVBmag: Manutude of delta V for B solutions
%   delVB: implusive delta V vector in heliocentric inertial coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delVAmag, delVBmag, delVA, delVB] = deltaV_ver2(r1_vec,orbit_matrix)

%Loading characteristic values
char_values = load_characteristic_values();

% Unpacking the orbit matrix semi major axes for required elements
a = orbit_matrix(:,1); % All semimajor axes from orbit matrix

%creating a vector of velocities in the cycler orbit perifocal frame 
v = sqrt(charValues.mu*((2./r1_vec)-(1./a))); 

%converting the velocities to heliocentric inertial frame





end