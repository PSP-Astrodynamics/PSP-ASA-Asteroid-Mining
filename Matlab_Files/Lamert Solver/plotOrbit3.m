function [x, y, z] = plotOrbit3(RAAN, inc, omega, p, e, theta_star, color, scale, grade, c,arrow,W)
    
    r_vec_xyz = zeros(3,length(theta_star));
    for n=1:length(theta_star)
        theta = theta_star(n) + omega;

        r_vec_rth = p / (1+ e * cos(theta - omega)) .* [1, 0, 0];

        ICR = [cos(RAAN)*cos(theta) - sin(RAAN)*cos(inc)*sin(theta), -cos(RAAN)*sin(theta)...
        - sin(RAAN)*cos(inc)*cos(theta), sin(RAAN)*sin(inc);
           sin(RAAN)*cos(theta) + cos(RAAN)*cos(inc)*sin(theta),...
           -sin(RAAN)*sin(theta) + cos(RAAN)*cos(inc)*cos(theta), -cos(RAAN)*sin(inc);
           sin(inc)*sin(theta), sin(inc)*cos(theta), cos(inc)];
    
        r_vec_xyz(:,n) = (ICR*r_vec_rth');
    end

    x = r_vec_xyz(1,:) + c(1);
    y = r_vec_xyz(2,:) + c(2);
    z = r_vec_xyz(3,:) + c(3);
    

    % plot3(x, y, z, 'Color', color, LineWidth=W)
    % hold on
    % if (arrow == 1)
    %     plotOrbitWithArrows(x, y, z, length(x)/10, color, scale, grade)
    % end
    % hold on
    
end