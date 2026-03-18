function plot_validation_results(target_f, WinData, data_real)
    % PLOT_VALIDATION_RESULTS: Visual comparison between MC and observations
    % target_f: Frequency for the title
    % WinData:  Matrix from MC simulation [Dist, W1, W2, W3]
    % data_real: The specific real data matrix for this frequency
    
    % 1. Process real data: log10(4*pi*r^2 * (E_win / E_ref))
    r_real = data_real(:,1);
    % Standard normalization using the 5th column (coda) as reference
    E_obs = log10(4*pi*r_real.^2 .* (data_real(:,2:4)./data_real(:,5)));

    % 2. Calculate vertical offset
    % We align the first MC window (WinData column 2) with real data at short distance
    offset = mean(E_obs(r_real < 30, 1)) - WinData(1,2);

    % 3. Visualization
    figure('Color', 'w', 'Name', sprintf('Validation %.1f Hz', target_f));
    hold on;
    colors = [0.2 0.2 1; 1 0.2 0.2; 0.2 0.8 0.2]; % Blue (W1), Red (W2), Green (W3)

    for w = 1:3
        % Plot observed data as points (HandleVisibility off to clean the legend)
        plot(r_real, E_obs(:,w), '.', 'Color', colors(w,:), 'MarkerSize', 4, 'HandleVisibility', 'off');
        % Plot MC simulation results as thick lines
        plot(WinData(:,1), WinData(:,w) + offset, 'Color', colors(w,:), 'LineWidth', 2.5, ...
             'DisplayName', sprintf('MC Window %d', w));
    end

    grid on; xlim([20 80]); % Focus on the standard study range [cite: 152, 286]
    xlabel('Epicentral Distance (km)'); 
    ylabel('log Normalized Energy');
    title(sprintf('Validation at %.1f Hz: Moho Model vs Obs. Data', target_f));
    legend('Location', 'northeast');
    hold off;
end