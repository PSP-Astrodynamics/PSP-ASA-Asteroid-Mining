function p_solution = lambertSolverP(a, c, s, r1, r2)
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

