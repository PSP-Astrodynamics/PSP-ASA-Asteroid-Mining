function [a_trans, c, s, r1, r2] = lambertCandidateSetup(r1_vec,r2_vec,mu_S,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Description: This function sets up candidate semimajor axes for
% our cycler orbit, which each have periods that are some integer multiples
% of Earth's period
% 
%   INPUT VARIABLES:
%   NAME : Description. (unit) 
%   r1_vec : Vector of the Earth's position with respect to the Sun. (km) 
%   r2_vec : Vector of the Target's position with respect to the Sun. (km)
%   mu_S : GM value of the Sun
%   n:  number of multiples of Earth years over which we want to generate
%       over (e.g. n = 10 generates candidates for orbits with a period of 10
%       Earth years
%
%   OUTPUT VARIABLES:
%   NAME : Description. (unit) 
%   a_trans : Matrix of Semimajor axis lengths. (km)
%   c : length of the line segment connecting (km)
%   s : semi-perimeter of the space triangle (km)
%   r1 : Vector of the Earth's position with respect to the Sun. (km)
%   r2 : Vector of the Target's position with respect to the Sun. (km)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% magnitude of vectors
r1 = norm(r1_vec);
r2 = norm(r2_vec);

c = norm(r2_vec - r1_vec); % Is the magnitude of vector from Earth to Target

s = 1/2 * (r1 + r2 + c); 

IP_E = 365.25*24*3600;              % Period of the Earth (s)
IP_trans = (IP_E : IP_E : n * IP_E);  % Creating a vector of Periods: Integer multiples of the Earth's Period for our cycler orbit's Periods


a_trans = (mu_S*(IP_trans./(2.*pi)).^2).^(1/3); %Converting candidate Periods into candidate semi-major axes for which we will generate our orbits
end