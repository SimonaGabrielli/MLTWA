function do_local_curve_plot(r, E1, E2, E3, x_min, x_centr, p)
% r, E1, E2, E3: Observed data (can be raw or smoothed)
% x_min: Minimum Absolute parameters (Black dashed line)
% x_centr: Centroid parameters (Orange solid line)
% p: sim_params

if isempty(r) || length(r) < 2
    title('No data available for this plot');
    return;
end

% Color Definitions
c1 = [0.3 0.3 1]; c2 = [1 0.2 0.2]; c3 = [0.1 0.7 0.1]; % Blue, Red, Green
c_orange = [1 0.5 0]; % Centroid color

hold on;

% 1. Plot OBSERVED DATA
plot(r, E1, '.', 'Color', c1, 'MarkerSize', 8, 'DisplayName', 'Obs E1'); 
plot(r, E2, '.', 'Color', c2, 'MarkerSize', 8, 'DisplayName', 'Obs E2');
plot(r, E3, '.', 'Color', c3, 'MarkerSize', 8, 'DisplayName', 'Obs E3');

r_m = linspace(min(r), max(r), 100)';

% 2. Plot CENTROID (Solid Orange - Base Model)
if ~isempty(x_centr) && all(~isnan(x_centr))
    Y1_c = computeEnergy3(x_centr, r_m, 1, p);
    Y2_c = computeEnergy3(x_centr, r_m, 2, p);
    Y3_c = computeEnergy3(x_centr, r_m, 3, p);
    plot(r_m, Y1_c, '-', 'Color', c_orange, 'LineWidth', 3, 'DisplayName', 'Centroid Fit');
    plot(r_m, Y2_c, '-', 'Color', c_orange, 'LineWidth', 3, 'HandleVisibility', 'off');
    plot(r_m, Y3_c, '-', 'Color', c_orange, 'LineWidth', 3, 'HandleVisibility', 'off');
end

% 3. Plot MINIMUM ABSOLUTE (Black Dashed - Comparison)
if ~isempty(x_min) && all(~isnan(x_min))
    Y1_m = computeEnergy3(x_min, r_m, 1, p);
    Y2_m = computeEnergy3(x_min, r_m, 2, p);
    Y3_m = computeEnergy3(x_min, r_m, 3, p);
    plot(r_m, Y1_m, 'k--', 'LineWidth', 2, 'DisplayName', 'Min Abs Fit');
    plot(r_m, Y2_m, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(r_m, Y3_m, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
end

% Aesthetics
grid off;
set(gca, 'Box', 'on', 'LineWidth', 1.1);
set(gca, 'FontSize', 14);
xlabel('Hypocentral Distance [km]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Log Normalized Energy', 'FontSize', 14, 'FontWeight', 'bold');
all_E = [E1; E2; E3];
ylim([min(all_E)-0.5, max(all_E)+0.5]);
legend('show', 'Location', 'best', 'FontSize', 12);
hold off;
end
