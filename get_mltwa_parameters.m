function [etai, g0] = get_mltwa_parameters(resultsTable, target_f)
    % Find frequency index
    idx = find(resultsTable.Freq_Hz == target_f);
    if isempty(idx), error('Frequency %.1f not found.', target_f); end
    
    B0_str = string(resultsTable.B0{idx});
    Le_inv_str = string(resultsTable.Le_inv{idx});
    
    parts_B0 = strsplit(B0_str, ' ');
    parts_Le = strsplit(Le_inv_str, ' ');
    
    B0_val = str2double(parts_B0(1));
    Le_inv_val = str2double(parts_Le(1));
    
    if isnan(B0_val) || isnan(Le_inv_val)
        error('Error in number conversion: check the table format.');
    end
    
    % calculation of etai and etas
    g0 = Le_inv_val * B0_val;
    etai = Le_inv_val * (1 - B0_val);
    
    fprintf('Bridge Success: f=%.1f, Le^-1=%.3f, B0=%.2f -> g0=%.4f, etai=%.4f\n', ...
            target_f, Le_inv_val, B0_val, g0, etai);
end