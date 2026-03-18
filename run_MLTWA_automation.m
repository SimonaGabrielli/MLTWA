function results_summary = run_MLTWA_automation(data_struct, sim_params, grid_settings, freqs, labels)
  
    
    results_summary = struct();
    fisher_levels = [1, 1.05, 1.10, 1.25];
    
    for i = 1:length(freqs)
        freq = freqs(i);
        label = labels{i};
        
        % Data extraction
        if ~isfield(data_struct, label), continue; end
        raw_data = data_struct.(label);
        rn=raw_data(:,1); E1=raw_data(:,2); E2=raw_data(:,3); E3=raw_data(:,4);
        
        mask = (E1>E2) & (E2>E3) & (E1>E3) & (rn >= sim_params.min_dist) & (rn <= sim_params.max_dist);
        if sum(mask)==0, continue; end
        
        [r_sort, idx] = sort(rn(mask));
        E1_s=E1(mask); E1_s=E1_s(idx);
        E2_s=E2(mask); E2_s=E2_s(idx);
        E3_s=E3(mask); E3_s=E3_s(idx);
        
        % Smoothing
        [r_sm, E1_sm, E2_sm, E3_sm] = prepare_smoothed_data(r_sort, E1_s, E2_s, E3_s, 50, 20);
        
        % Grid Search
        fprintf(' -> Working on frequency: %.1f Hz...', freq);
        res_RAW = gridSearch_Core(r_sort(1:5:end), E1_s(1:5:end), E2_s(1:5:end), E3_s(1:5:end), sim_params, grid_settings);
        res_SMOOTH = gridSearch_Core(r_sm, E1_sm, E2_sm, E3_sm, sim_params, grid_settings);
        
        results_summary.(label).RAW = res_RAW;
        results_summary.(label).SMOOTH = res_SMOOTH;
        
        % Plotting
        fig = figure('Name', ['Report ' label], 'Color', 'w', 'Position', [50, 50, 1200, 900], 'Visible', 'off');
        
        % 1. Fisher RAW (Top Left)
        subplot(2,2,1); 
        do_local_fisher_plot(res_RAW, fisher_levels, 'Fisher Map - RAW');
        
        % 2. Fisher SMOOTH (Top Right)
        subplot(2,2,2); 
        do_local_fisher_plot(res_SMOOTH, fisher_levels, 'Fisher Map - SMOOTH');
        
        % 3. Energies RAW + RAW Model (Bottom Left)
        subplot(2,2,3);
        do_local_curve_plot(r_sort, E1_s, E2_s, E3_s, res_RAW.x_opt, [], sim_params);
        title('Raw Data vs Theoretical Fit');
        
        % 4. Energies SMOOTH + SMOOTH Model (Bottom Right)
        subplot(2,2,4);
        do_local_curve_plot(r_sm, E1_sm, E2_sm, E3_sm, res_SMOOTH.x_opt, [], sim_params);
        title('Smoothed Data vs Theoretical Fit');
        
        sgtitle(['MLTWA Analysis - ' num2str(freq) ' Hz'], 'FontSize', 18, 'FontWeight', 'bold');
        
        % Salvataggio
        saveas(fig, ['Report__' label '.png']);
        savefig(fig, ['Report__' label '.fig']);
        close(fig);
        fprintf(' (PNG 2x2 saved).\n');
    end
    fprintf('--- COMPLETED ---\n');
end