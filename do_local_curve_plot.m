function do_local_curve_plot(r, E1, E2, E3, x_raw, x_sm, p)

if isempty(r) || length(r) < 2
        title('No data available for this plot');
        return;
end
% Plot Energies vs Distance
plot(r, E1, '.', 'Color', [0.3 0.3 1]); hold on;
plot(r, E2, '.', 'Color', [1 0.2 0.2]);
plot(r, E3, '.', 'Color', [0.1 0.7 0.1]);

r_m = linspace(min(r), max(r), 100)';

% Theoretical energy curve - unsmoothed data

if ~isempty(x_raw) && all(~isnan(x_raw))
    Y1_r = computeEnergy3(x_raw, r_m, 1, p);
    Y2_r = computeEnergy3(x_raw, r_m, 2, p);
    Y3_r = computeEnergy3(x_raw, r_m, 3, p);
    plot(r_m, Y1_r, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Raw Fit');
    plot(r_m, Y2_r, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(r_m, Y3_r, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% Theoretical energy curve - smoothed data
if ~isempty(x_sm) && all(~isnan(x_sm))
    Y1_s = computeEnergy3(x_sm, r_m, 1, p);
    Y2_s = computeEnergy3(x_sm, r_m, 2, p);
    Y3_s = computeEnergy3(x_sm, r_m, 3, p);
    plot(r_m, Y1_s, 'b-', 'LineWidth', 2, 'DisplayName', 'Smooth Fit');
    plot(r_m, Y2_s, 'b-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(r_m, Y3_s, 'b-', 'LineWidth', 2, 'HandleVisibility', 'off');
end

title('Curve Fit', 'FontSize', 12);
xlabel('Distance [km]'); ylabel('Energy');
ylim([min(E3)-1, max(E1)+1]);
% xlim([20, 80]);
set(gca, 'YScale', 'log'); % MLTWA richiede quasi sempre la scala log
legend('show', 'Location', 'best');
axis square; grid off;
end