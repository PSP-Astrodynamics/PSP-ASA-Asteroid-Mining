function [] = lambertSolverTOFvsSMA(R,mu,c,s)

    %R is characteristic scale, typically radius of planet
    
    a_range = R*linspace(0,18,10000);
    
    a_min = s/2;
    TOF_min = lambertSolverTOF(a_min, c, s, mu);
    TOF_min1 = TOF_min{2,1};
    TOF_min2 = TOF_min{2,3};
    
    % '1A', '1B', '2A', '2B', '1P', '2P', '1H', '2H'
    TOF_range = zeros(length(a_range), 8);
    for n=1:length(a_range)
        TOF_sol = lambertSolverTOF(a_range(n), c, s, mu);
        TOF_range(n,:) = [TOF_sol{2,1},TOF_sol{2,2},TOF_sol{2,3},TOF_sol{2,4},...
            TOF_sol{2,5},TOF_sol{2,6},TOF_sol{2,7},TOF_sol{2,8}];
    end
    
    for i=1:length(a_range)
        for j=1:8
            TOF_element = TOF_range(i,j);
            if (abs(imag(TOF_element))>0.5)
                TOF_range(i,j) = NaN;
                %TOF_range(i,j) = abs(TOF_range(i,j));
            end
        end
    end

    figure
    
    interval = 100;
    colors = lines(8);

    % Plot all curves first
    plot(a_range, TOF_range(:,1)/3600, 'LineWidth', 1.5, 'Color', colors(1,:)); hold on;
    plot(a_range, TOF_range(:,2)/3600, 'LineWidth', 1.5, 'Color', colors(2,:));
    plot(a_range, TOF_range(:,3)/3600, 'LineWidth', 1.5, 'Color', colors(3,:));
    plot(a_range, TOF_range(:,4)/3600, 'LineWidth', 1.5, 'Color', colors(4,:));
    plot(a_range, TOF_range(:,5)/3600, 'LineWidth', 1.5, 'Color', colors(5,:));
    plot(a_range, TOF_range(:,6)/3600, 'LineWidth', 1.5, 'Color', colors(6,:));
    plot(a_range, TOF_range(:,7)/3600, 'LineWidth', 1.5, 'Color', colors(7,:));
    plot(a_range, TOF_range(:,8)/3600, 'LineWidth', 1.5, 'Color', colors(8,:));
    
    % Plot all markers next
    plot(a_range(1:interval:end), TOF_range(1:interval:end,1)/3600, 'o', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(1,:)); hold on;
    plot(a_range(1:interval:end), TOF_range(1:interval:end,2)/3600, 's', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(2,:));
    plot(a_range(1:interval:end), TOF_range(1:interval:end,3)/3600, '*', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(3,:));
    plot(a_range(1:interval:end), TOF_range(1:interval:end,4)/3600, 'd', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(4,:));
    plot(a_range(1:interval:end), TOF_range(1:interval:end,5)/3600, 'x', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(5,:));
    plot(a_range(1:interval:end), TOF_range(1:interval:end,6)/3600, '^', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(6,:));
    plot(a_range(1:interval:end), TOF_range(1:interval:end,7)/3600, '+', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(7,:));
    plot(a_range(1:interval:end), TOF_range(1:interval:end,8)/3600, '>', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors(8,:));

    % Plot the horizontal lines for TOF_min1 and TOF_min2
    yline(TOF_min1/3600, 'k--', 'LineWidth', 1.5, 'DisplayName', 'TOF_{min,1}');
    yline(TOF_min2/3600, 'k-.', 'LineWidth', 1.5, 'DisplayName', 'TOF_{min,2}');

    % Formatting the plot
    grid on
    xlabel('Semimajor Axis [km]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Time of Flight [hrs]', 'FontSize', 12, 'FontWeight', 'bold');

    % 18 hours for earth
    ylim([0, 3*TOF_min1/3600]);
    xlim([0, 10*R])
    
    % Add a legend with customized labels
    legend({'','','','','','','','','1A', '1B', '2A', '2B', '1P', '2P', '1H', '2H', 'TOF_{min,1}', 'TOF_{min,2}'}, ...
           'Location', 'northwest', 'FontSize', 10);
    
    % Add title (optional)
    title('Time of Flight vs Semimajor Axis', 'FontSize', 14, 'FontWeight', 'bold');
end
