function Results = extract_windows(p, dists, filename_histo, tW, tc)
    % EXTRACT_PISCONTI_WINDOWS (Updated with integrated plotting)
    ET = load(filename_histo);
    Nd = length(dists);
    Nt = size(ET, 1);
    
    % Time vector based on histogram bins (0 to tmax)
    time_vec = (0:Nt-1)'; 
    v_S = p.A1; 
    Results = zeros(Nd, 4);
    Results(:, 1) = dists(:); 
    
    for nd = 1:Nd
        r = dists(nd);
        t_S = r / v_S; % S-wave arrival time
        
        % --- Define Time Windows (basate sull'input tW) ---
        mask1 = (time_vec >= t_S) & (time_vec < t_S + tW);
        mask2 = (time_vec >= t_S + tW) & (time_vec < t_S + 2*tW);
        mask3 = (time_vec >= t_S + 2*tW) & (time_vec < t_S + 3*tW);
        
        % Coda Reference Window: da tc a tc + tW (origin time) 
        mask_coda = (time_vec >= tc) & (time_vec < tc + tW);
        
        % 2. Sum energy in the windows (Integrals)
        E1 = sum(ET(mask1, nd));
        E2 = sum(ET(mask2, nd));
        E3 = sum(ET(mask3, nd));
        Ecoda = sum(ET(mask_coda, nd));
        
        geo_spread = 4 * pi * r^2;
        
        % 3. Coda Normalization (log10 scale)
        if Ecoda > 0
            Results(nd, 2) = log10(geo_spread * (E1 / Ecoda));
            Results(nd, 3) = log10(geo_spread * (E2 / Ecoda));
            Results(nd, 4) = log10(geo_spread * (E3 / Ecoda));
        else
            Results(nd, 2:4) = NaN;
        end
    end

    % --- INTEGRATED PLOTTING BLOCK ---
    figure('Color', 'w'); hold on;
    
    % Results(:,1) è la distanza, (:,2) Win1, (:,3) Win2, (:,4) Win3
    plot(Results(:,1), Results(:,2), 'bs-', 'LineWidth', 1.5, 'DisplayName', sprintf('Win 1 (0-%.0fs)', tW));
    plot(Results(:,1), Results(:,3), 'rs-', 'LineWidth', 1.5, 'DisplayName', sprintf('Win 2 (%.0f-%.0fs)', tW, 2*tW));
    plot(Results(:,1), Results(:,4), 'gs-', 'LineWidth', 1.5, 'DisplayName', sprintf('Win 3 (%.0f-%.0fs)', 2*tW, 3*tW));

    % Dynamic axis limits based on distances
    xlim([min(dists)-5, max(dists)+5]);
    
    xlabel('Source-Receiver Distance (km)');
    ylabel('log Normalized Energy Integral');
    title('Monte Carlo Energy Windows - Layered Model');
    legend('Location', 'northeast');
    grid on;
    hold off;
    
    fprintf('Extraction and Plotting completed for %d distances (tW=%ds, tc=%ds).\n', Nd, tW, tc);
end