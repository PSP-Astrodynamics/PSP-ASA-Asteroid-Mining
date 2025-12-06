function change_axis_from_Km_to_AU(x_axis, y_axis)
    % Astronomical Units in km
    AU = 149597870.7; % km    
    
    if x_axis
        % --- Get current x-axis limits in km and convert to AU ---
        x_km_limits = xlim;
        x_au_min = floor(x_km_limits(1) / AU);
        x_au_max = ceil(x_km_limits(2) / AU);
        
        % --- Generate tick values from floor(min) to ceil(max), always including 0 ---
        x_ticks_au = x_au_min : 3 : x_au_max;
        x_ticks_km = x_ticks_au * AU;
        
        % --- Apply ticks and labels ---
        set(gca, 'XTick', x_ticks_km);
        set(gca, 'XTickLabel', string(x_ticks_au));
        
        % --- Force axis limits to match AU tick range exactly ---
        xlim([x_au_min, x_au_max] * AU);
    end
    
    if y_axis
        % --- Get current y-axis limits in km and convert to AU ---
        y_km_limits = ylim;
        y_au_min = floor(y_km_limits(1) / AU);
        y_au_max = ceil(y_km_limits(2) / AU);
        
        % --- Generate tick values from floor(min) to ceil(max), always including 0 ---
        y_ticks_au = y_au_min : 2 : y_au_max;
        y_ticks_km = y_ticks_au * AU;
        
        % --- Apply ticks and labels ---
        set(gca, 'YTick', y_ticks_km);
        set(gca, 'YTickLabel', string(y_ticks_au));
        
        % --- Force axis limits to match AU tick range exactly ---
        ylim([y_au_min, y_au_max] * AU);
    end
end