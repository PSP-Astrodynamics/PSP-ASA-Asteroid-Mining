function p_solution = lambertSolverP(a, c, s, r1, r2)
    alpha0 = 2 * real(asin(sqrt(s / (2 * a))));  % radians
    beta0 = 2 * asin(sqrt((s - c) / (2 * a)));  % radians

    alphaH = 2 * asinh(sqrt(s / (2 * a))); 
    betaH = 2 * asinh(sqrt((s - c) / (2 * a)));

    pAB_1 = 4*a*(s-r1)*(s-r2)*(sin(0.5*(alpha0+beta0))^2)/c^2;
    pAB_2 = 4*a*(s-r1)*(s-r2)*(sin(0.5*(alpha0-beta0))^2)/c^2;

    pH_1 = 4*a*(s-r1)*(s-r2)*(sinh(0.5*(alphaH+betaH))^2)/c^2;
    pH_2 = 4*a*(s-r1)*(s-r2)*(sinh(0.5*(alphaH-betaH))^2)/c^2;

    p_solution = {'pAB_1', 'pAB_2', 'pH_1', 'pH_2'; pAB_1, pAB_2, pH_1, pH_2};
end

