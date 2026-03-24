function Results = extract_windows(p, dists, filename_histo, tW, tc, data_real, target_freq)
    % 1. Initialization and Defaults
    % If target_freq is not provided, set a default or dummy value for the plot title
    if nargin < 7
        target_freq = NaN; 
    end
    if nargin < 6
        data_real = []; 
    end

    % 2. Load Monte Carlo simulation (HCTM Model) 
    ET_data = load(filename_histo);
    Nd = length(dists);
    Nt = size(ET_data, 1);
    time_vec = (0:Nt-1)'; 
    v_S = p.A1; 
    
    Results = zeros(Nd, 4);
    Results(:, 1) = dists(:); 
    
    for nd = 1:Nd
        r = dists(nd);
        t_S = r / v_S; 
        
        % --- Time Windows Definition (MLTWA standard)  ---
        mask1 = (time_vec >= t_S) & (time_vec < t_S + tW);
        mask2 = (time_vec >= t_S + tW) & (time_vec < t_S + 2*tW);
        mask3 = (time_vec >= t_S + 2*tW) & (time_vec < t_S + 3*tW);
        mask_coda = (time_vec >= tc) & (time_vec < tc + tW); 
        
        % 3. Sum energy particles and apply Geometrical Spreading 
        E1 = sum(ET_data(mask1, nd));
        E2 = sum(ET_data(mask2, nd));
        E3 = sum(ET_data(mask3, nd));
        Ecoda = sum(ET_data(mask_coda, nd));
        
        geo_spread = 4 * pi * r^2; 
        
        if Ecoda > 0
            Results(nd, 2) = log10(geo_spread * (E1 / Ecoda)); 
            Results(nd, 3) = log10(geo_spread * (E2 / Ecoda));
            Results(nd, 4) = log10(geo_spread * (E3 / Ecoda));
        else
            Results(nd, 2:4) = NaN;
        end
    end

    % --- INTEGRATED PLOTTING BLOCK ---
    figName = 'Theoretical Energies (HCTM)';
    if ~isnan(target_freq)
        figName = sprintf('Theoretical vs Observed - %.1f Hz', target_freq);
    end
    
    figure('Color', 'w', 'Name', figName); 
    hold on;
    
    % Plot MC Simulation (Layered HCTM) as Lines
    plot(Results(:,1), Results(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'MC Win 1');
    plot(Results(:,1), Results(:,3), 'r-', 'LineWidth', 2, 'DisplayName', 'MC Win 2');
    plot(Results(:,1), Results(:,4), 'g-', 'LineWidth', 2, 'DisplayName', 'MC Win 3');
    
    % Optional: Plot Observed Data (UHS) if provided 
    if ~isempty(data_real) && ~isnan(target_freq)
        field_name = sprintf('f%.1fHz', target_freq);
        field_name = strrep(field_name, '.', '_'); 

        if isfield(data_real, field_name)
            target_data = data_real.(field_name);
            if isstruct(target_data) || istable(target_data)
                scatter(target_data.distances, target_data.E1, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility', 'off');
                scatter(target_data.distances, target_data.E2, 40, 'r', 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility', 'off');
                scatter(target_data.distances, target_data.E3, 40, 'g', 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility', 'off');
            else
                scatter(target_data(:,1), target_data(:,2), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility', 'off');
                scatter(target_data(:,1), target_data(:,3), 40, 'r', 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility', 'off');
                scatter(target_data(:,1), target_data(:,4), 40, 'g', 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility', 'off');
            end
        end
    end

    xlim([min(dists)-5, max(dists)+5]);
    xlabel('Distance (km)');
    ylabel('log_{10} [Norm. Energy \cdot 4\pi r^2]');
    title(figName);
    legend('Location', 'best');
    grid off;
    hold off;
end
