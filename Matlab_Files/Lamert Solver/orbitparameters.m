function [Omega, theta, i] = orbitparameters(r1, r2)
    
    r1mag = norm(r1);
    r1hat = r1/r1mag;

    htrans = cross(r1, r2);
    htransmag = norm(htrans);
    hhat = htrans/htransmag;
    
    thetahat = cross(hhat, r1hat);

    i = acos(hhat(3));
    
    Omegatrans1 = asin(hhat(1)/sin(i));
    Omegatrans2 = acos(-hhat(2)/sin(i));
    Omega = nonuniqueAngle([Omegatrans1, pi - Omegatrans1, Omegatrans2, 2*pi - Omegatrans2]);

    theta_1 = asin(r1hat(3)/sin(i));
    theta_2 = acos(thetahat(3)/sin(i));
    theta = nonuniqueAngle([theta_1, pi-theta_1, theta_2, 2 * pi - theta_2]);
end