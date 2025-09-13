function plotCycler(eccentricity, r1_vec, r2_vec, r1, p, color)

    theta_star = acos(((p / r1) - 1) / eccentricity);
    
    [Omega, theta, i] = orbitparameters(r1_vec, r2_vec);
    aop = theta - theta_star;
    
    plotOrbit3(Omega, i, aop, p, eccentricity, linspace(0,2*pi,1000), color, 1, 1, [0,0,0],0,1.5)

end