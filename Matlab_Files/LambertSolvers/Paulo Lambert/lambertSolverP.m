function p_solution = lambertSolverP(a, c, s, r1, r2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Description:
% This function outputs values for the Semi-latus rectum (p): the distance 
% from the foci to the edge of the ellipse following a line perpendicular 
% to the semi-major axis (a), defining the shape of an elliptical orbit.
% 
%  There are only three truly independant inputs: a, r1, and r2 as c and s
% depend on these. For these values, it is possible to form two ellipses,
% as there are exactly two possible positions for the empty foci. Therefore
% this function will output two unique p values.
%
%   INPUT VARIABLES:
% Units are all scalar distances, usually AU.
%
%   a : Semimajor axis of elliptical orbit.
%   c : The straight line distance between first position on the orbit p1 
% and second position on the obrbit p2. mag(r2 - r1). 
%   s : The semiperimeter, calculated as the sum of (mag(r1)
%   + mag(r2) + c)/2. It's half the perimeter of the triangle formed between
%   the sun, p1, and p2.
%   r1 : Position vector to the spacecraft's first position p1 with respect 
% to the sun.
%   r2 : Position vector of the spacecraft's second position p2 with respect
%  to the sun.
%
%   OUTPUT VARIABLES:   
%   pAB_1 : P for ellipse 1
%   pAB_2 : P for ellipse 2
%   pAB_1other : Repeated value for 1?
%   pAB_2other : Repeated value for 2?
%
%
%   NOTES:
%   * pAB_1other and PAB_2other appear to be redundant.
%   * Commented out lines of code refer to hyperbolic orbits.
%
% Code header written by Kai s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alpha0 = 2 * real(asin(sqrt(s / (2 * a))));  % radians
    beta0 = 2 * asin(sqrt((s - c) / (2 * a)));  % radians
    alphaOther = 2 * pi - alpha0;
    betaOther = 2 * pi - beta0;
    %alphaH = 2 * asinh(sqrt(s / (2 * a))); 
    %betaH = 2 * asinh(sqrt((s - c) / (2 * a)));

    pAB_1 = 4*a*(s-r1)*(s-r2)*(sin(0.5*(alpha0+beta0))^2)/c^2;
    pAB_2 = 4*a*(s-r1)*(s-r2)*(sin(0.5*(alpha0-beta0))^2)/c^2;
    pAB_1other = 4*a*(s-r1)*(s-r2)*(sin(0.5*(alphaOther+betaOther))^2)/c^2;
    pAB_2other = 4*a*(s-r1)*(s-r2)*(sin(0.5*(alphaOther-betaOther))^2)/c^2;

    %pH_1 = 4*a*(s-r1)*(s-r2)*(sinh(0.5*(alphaH+betaH))^2)/c^2;
    %pH_2 = 4*a*(s-r1)*(s-r2)*(sinh(0.5*(alphaH-betaH))^2)/c^2;

    p_solution = {'pAB_1', 'pAB_2', 'pAB_1other', 'pAB_2other'; pAB_1, pAB_2, pAB_1other, pAB_2other};
end

