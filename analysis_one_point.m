function mse = analysis_one_point(x1, x2, data, params)

E1_raw = data.E1;
E2_raw = data.E2;
E3_raw = data.E3;
rn_raw = data.rn;

if isfield(params, 'min_dist')
    min_r = params.min_dist;
else
    min_r = 0;
end

if isfield(params, 'max_dist')
    max_r = params.max_dist;
else
    max_r = Inf;
end

valid_rows = (E1_raw > E2_raw) & ...
    (E2_raw > E3_raw) & ...
    (E1_raw > E3_raw) & ...
    (rn_raw >= min_r) & ...
    (rn_raw <= max_r);


% Filtering the rows fulfilling the condition
E1s = E1_raw(valid_rows);
E2s = E2_raw(valid_rows);
E3s = E3_raw(valid_rows);
r = rn_raw(valid_rows);

cut_percent = 0.5; % Remove the 0.5% higher and lower
   
    lim_low_1 = prctile(E1s, cut_percent);      lim_high_1 = prctile(E1s, 100-cut_percent);
    lim_low_2 = prctile(E2s, cut_percent);      lim_high_2 = prctile(E2s, 100-cut_percent);
    lim_low_3 = prctile(E3s, cut_percent);      lim_high_3 = prctile(E3s, 100-cut_percent);
    
    % (1 = to keep, 0 = to remove)
    mask1 = (E1s >= lim_low_1) & (E1s <= lim_high_1);
    mask2 = (E2s >= lim_low_2) & (E2s <= lim_high_2);
    mask3 = (E3s >= lim_low_3) & (E3s <= lim_high_3);
    
    % We keep the data with all three conditions
    keep_idx = mask1 & mask2 & mask3;
    
    E1s = E1s(keep_idx);
    E2s = E2s(keep_idx);
    E3s = E3s(keep_idx);
    r   = r(keep_idx);
    
    n_removed = sum(~keep_idx);
    if n_removed > 0
        fprintf('Number of outliers removed (top/bottom %.1f%%)\n', cut_percent);
    end
    % --------------------------------------------------

    if isempty(r)
        error('No data in the range %.0f - %.0f km', min_r, max_r);
    end

% Sorting r
[r_sorted, sortedIdx] = sort(r);

x(1) = x1; % B0
x(2) = x2; % Le

mse = main(x, r_sorted, sortedIdx, E1s, E2s, E3s, params);


    function U = main(x, r, sortid, E1, E2, E3, p)
        Y1 = computeEnergy3(x, r, 1, p);
        Y2 = computeEnergy3(x, r, 2, p);
        Y3 = computeEnergy3(x, r, 3, p);

        E1(isinf(E1)) = NaN;
        E1 = fillmissing(E1, "previous");
        U1 = mean((Y1 - E1(sortid)).^2);

        E2(isinf(E2)) = NaN;
        E2 = fillmissing(E2, "previous");
        U2 = mean((Y2 - E2(sortid)).^2);

        E3(isinf(E3)) = NaN;
        E3 = fillmissing(E3, "previous");
        U3 = mean((Y3 - E3(sortid)).^2);

        %MSE Total
        U = U1 + U2 + U3;

        % Visualization
        figure;

        h1 = plot(r, E1(sortid), 'b.','MarkerSize', 9, 'DisplayName', 'Obs E1');
        hold on;
        h2 = plot(r, E2(sortid), '.', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 9, 'DisplayName', 'Obs E2');
        h3 = plot(r, E3(sortid), 'g.','MarkerSize', 9, 'DisplayName', 'Obs E3');

        plot(r, Y1, 'k-', 'LineWidth', 2, 'DisplayName', 'Syn Model'); % Ho messo nero per contrasto, o usa 'b-'
        plot(r, Y2, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
        plot(r, Y3, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');

        if isinf(max_r)
            x_limit_up = max(r) + 5;
        else
            x_limit_up = max_r + 5;
        end

        if min_r == 0
            x_limit_d = min(r) - 5;
        else
            x_limit_d = min_r - 5;
        end

        % xlim([x_limit_d, x_limit_up]);
        % ylim([-1 10]);

        xlabel("Hypocentral Distance r [km]", 'FontSize', 12);
        ylabel("Log Avg Normalized Energy", 'FontSize', 12);
        % title(['Single Point Check (MSE: ' num2str(U, '%.4f') ')'], 'FontSize', 14);
        legend([h1 h2 h3], 'Location', 'best');
        grid on;
        hold off;

        set(gcf, 'Renderer', 'painters');
        set(gcf, 'Color', 'white');

        saveas(gcf, 'Energies_plot_SinglePoint.png');
    end

end