function [Omega, theta, i] = orbitparameters(r1_vec, r2_vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Description: This function sets up candidate semimajor axes for
% our cycler orbit, which each have periods that are some integer multiples
% of Earth's period
% 
%   INPUT VARIABLES:
%   NAME : Description. (unit) 
%   r1_vec : Vector of the Earth's position with respect to the Sun. (km) 
%   r2_vec : Vector of the Target's position with respect to the Sun. (km)
%
%   OUTPUT VARIABLES:
%   NAME : Description. (unit) 
%   Omega : The angle from the positive x axis to the direction in which the orbit intersects the xy plane (radians)
%   theta : The true anomaly of r1 from the direction in which the orbit intersects the xy plane (radians)
%   i : The inclination of the plane containing both r1 and r2 with respect to the z axis (radians)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    r1mag = norm(r1_vec);
    r1hat = r1_vec/r1mag;

    htrans = cross(r1_vec, r2_vec);
    htransmag = norm(htrans);
    hhat = htrans/htransmag; % unit vector normal to both r1 and r2
    
    thetahat = cross(hhat, r1hat);

    i = acos(hhat(3));
    
    Omegatrans1 = asin(hhat(1)/sin(i));
    Omegatrans2 = acos(-hhat(2)/sin(i));
    Omega = nonuniqueAngle([Omegatrans1, pi - Omegatrans1, Omegatrans2, 2*pi - Omegatrans2]);

    theta_1 = asin(r1hat(3)/sin(i));
    theta_2 = acos(thetahat(3)/sin(i));
    theta = nonuniqueAngle([theta_1, pi-theta_1, theta_2, 2 * pi - theta_2]);
end