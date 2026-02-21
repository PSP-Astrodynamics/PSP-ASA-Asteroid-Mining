function [] = cost_vs_sma(orbit_matrix, cost)
semi_major_axis = orbit_matrix(:,1);

figure(n)
plot(semi_major_axis, cost, 'or')
change_axis_from_Km_to_AU(true, false)
xlabel('Semi-Major Axis (AU)')
ylabel('Cost')
title('Semi-Major Axis vs. Cost')
grid on;
end