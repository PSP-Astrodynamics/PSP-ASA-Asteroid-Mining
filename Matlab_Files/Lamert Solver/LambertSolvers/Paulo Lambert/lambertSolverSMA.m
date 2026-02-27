function a_solutions = lambertSolverSMA(TOF, c, s, mu)
    % Define alpha0 and beta0
    alpha0 = @(a) 2 * asin(sqrt(s ./ (2 * a)));  % radians
    beta0 = @(a) 2 * asin(sqrt((s - c) ./ (2 * a)));  % radians
    
    % Define the four equations
    lambert_eq1A = @(a) sqrt(mu) * TOF - a^(3/2) * (alpha0(a) - sin(alpha0(a)) - (beta0(a) - sin(beta0(a))));
    lambert_eq1B = @(a) sqrt(mu) * TOF - a^(3/2) * (2*pi - (alpha0(a) - sin(alpha0(a))) - (beta0(a) - sin(beta0(a))));
    lambert_eq2A = @(a) sqrt(mu) * TOF - a^(3/2) * (alpha0(a) - sin(alpha0(a)) + beta0(a) - sin(beta0(a)));
    lambert_eq2B = @(a) sqrt(mu) * TOF - a^(3/2) * (2*pi - (alpha0(a) - sin(alpha0(a))) + (beta0(a) - sin(beta0(a))));
    
    % Solve for 'a' using fzero for each equation
    a_guess = 10*s;
    
    % answer key: [1A; 1B; 2A; 2B; 1H; 2H]
    a_solutions = real([fsolve(lambert_eq1A, a_guess), fsolve(lambert_eq1B, a_guess),...
        fsolve(lambert_eq2A, a_guess), fsolve(lambert_eq2B, a_guess)]);

    a_solutions = {'1A', '1B', '2A', '2B';...
        a_solutions(1),a_solutions(2),a_solutions(3),a_solutions(4)};
    
end