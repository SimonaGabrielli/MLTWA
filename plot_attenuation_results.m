function plot_attenuation_results(ResultsTable)
    % plot_attenuation_results: figure with B0, Le and Q-inv vs frequency
    freqs = ResultsTable.Freq_Hz;
    n = length(freqs);
    
    B0 = zeros(n, 1); B0_err = zeros(n, 1);
    Le = zeros(n, 1); Le_err = zeros(n, 1);
    Qi = zeros(n, 1); Qi_err = zeros(n, 1);
    Qs = zeros(n, 1); Qs_err = zeros(n, 1);
    Qt = zeros(n, 1); Qt_err = zeros(n, 1);
    
    % Helper to parse "value ± error" strings
    parse_val = @(str) sscanf(strrep(str, ' ± ', ' '), '%f %f');
    for i = 1:n
        res_B0 = parse_val(ResultsTable.B0(i));
        B0(i) = res_B0(1); B0_err(i) = res_B0(2);
        res_Le = parse_val(ResultsTable.Le_inv(i));
        Le(i) = res_Le(1); Le_err(i) = res_Le(2);
        res_i = parse_val(ResultsTable.Qi_inv(i));
        Qi(i) = res_i(1); Qi_err(i) = res_i(2);
        res_s = parse_val(ResultsTable.Qs_inv(i));
        Qs(i) = res_s(1); Qs_err(i) = res_s(2);
        res_t = parse_val(ResultsTable.Qt_inv(i));
        Qt(i) = res_t(1); Qt_err(i) = res_t(2);
    end
    
    % --- Q0 and n (Q = Q0*f^n) ---
    % Fit log10(Q^-1) = intercept - n*log10(f)
    % Q0 = 1 / 10^(intercept)
    fit_q = @(f, y_inv) polyfit(log10(f), log10(y_inv), 1);

    p_i = fit_q(freqs, Qi); n_i = -p_i(1); Q0_i = 1 / (10^p_i(2));
    p_s = fit_q(freqs, Qs); n_s = -p_s(1); Q0_s = 1 / (10^p_s(2));
    p_t = fit_q(freqs, Qt); n_t = -p_t(1); Q0_t = 1 / (10^p_t(2));

    % function [Q0, n, dQ0, dn] = fit_with_errors(x, y_inv)
    %     logx = log10(x);
    %     logy = log10(y_inv);
    %     [p, S] = polyfit(logx, logy, 1);
    %     % Errori sui coefficienti del fit lineare (p1 è la pendenza, p2 è l'intercetta)
    %     % Usiamo la matrice di covarianza derivata dalla struttura S
    %     stats = sqrt(diag(inv(S.R)*inv(S.R)')*S.normr^2/S.df);
    %     dn = stats(1);
    %     da = stats(2);
    % 
    %     n = -p(1);
    %     Q0 = 1 / (10^p(2));
    %     % Errore propagato per Q0 (derivata di 10^-a)
    %     dQ0 = abs(log(10) * Q0 * da);
    % end
    % 
    % [Q0_i, n_i, dQ0_i, dn_i] = fit_with_errors(freqs, Qi);
    % [Q0_s, n_s, dQ0_s, dn_s] = fit_with_errors(freqs, Qs);
    % [Q0_t, n_t, dQ0_t, dn_t] = fit_with_errors(freqs, Qt);

    fig = figure('Color', 'w', 'Name', 'Attenuation Analysis', 'Position', [100, 100, 1000, 800]);
    
    % --- B0 ---
    subplot(2, 2, 1);
    hold on; grid off; box on;
    set(gca, 'XScale', 'log', 'TickDir', 'out', 'FontSize', 12);
    errorbar(freqs, B0, B0_err, 'mo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
    ylabel('Seismic Albedo (B_0)', 'FontSize', 14, 'FontWeight', 'bold');   
    xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    ylim([0 1]); xlim([min(freqs)*0.8, max(freqs)*1.2]);
    xticks(freqs); xticklabels(string(freqs));
    title('Seismic Albedo', 'FontSize', 15);
    
    % --- Le^-1 ---
    subplot(2, 2, 2);
    hold on; grid off; box on;
    set(gca, 'XScale', 'log', 'TickDir', 'out', 'FontSize', 12);
    errorbar(freqs, Le, Le_err, 'go-', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
    ylabel('L_e^{-1} [km^{-1}]', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([min(freqs)*0.8, max(freqs)*1.2]);
    xticks(freqs); xticklabels(string(freqs));
    title('Extinction Length Inverse', 'FontSize', 15);
    
    % --- Q^-1 (LOG-LOG) ---
    subplot(2, 2, 3); 
    hold on; grid off; box on;
    set(gca, 'XScale', 'log', 'YScale', 'log', 'TickDir', 'out', 'FontSize', 12);
    % 
    lab_t = sprintf('Total: Q_t(f) = %.0f f^{%.2f}', Q0_t, n_t);
    lab_i = sprintf('Intrinsic: Q_i(f) = %.0f f^{%.2f}', Q0_i, n_i);
    lab_s = sprintf('Scattering: Q_s(f) = %.0f f^{%.2f}', Q0_s, n_s);

    % lab_t = sprintf('Total: Q_t(f) = (%.0f \\pm %.0f) f^{%.2f \\pm %.2f}', Q0_t, dQ0_t, n_t, dn_t);
    % lab_i = sprintf('Intrinsic: Q_i(f) = (%.0f \\pm %.0f) f^{%.2f \\pm %.2f}', Q0_i, dQ0_i, n_i, dn_i);
    % lab_s = sprintf('Scattering: Q_s(f) = (%.0f \\pm %.0f) f^{%.2f \\pm %.2f}', Q0_s, dQ0_s, n_s, dn_s);
    
    errorbar(freqs, Qt, Qt_err, 'ks-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', lab_t);
    errorbar(freqs, Qi, Qi_err, 'r^-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', lab_i);
    errorbar(freqs, Qs, Qs_err, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', lab_s);
    
    xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Attenuation (Q^{-1})', 'FontSize', 14, 'FontWeight', 'bold');
    title('Frequency Dependence of Attenuation Factors', 'FontSize', 15);
    legend('Location', 'southwest', 'FontSize', 12, 'FontWeight', 'normal'); % Font 12 per far stare tutto
    
    ylim([1e-4, 5e-2]); xlim([min(freqs)*0.8, max(freqs)*1.2]);
    xticks(freqs); xticklabels(string(freqs));
    
    filename = 'Attenuation_Summary.png';
    exportgraphics(fig, filename, 'Resolution', 300);
    savefig(fig, 'Attenuation_Summary.fig')
    fprintf('Summary figure saved as: %s\n', filename);
    
    hold off;
end