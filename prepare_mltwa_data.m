function data_out = prepare_mltwa_data(raw_file, freqs, labels, is_data_processed)
    % prepare_mltwa_data: Loads energies and applies log10/geom spreading 
    % ONLY if the user specifies that the data is raw.
    
    data_out = struct();
    fprintf('\n--- Loading and Processing Energies ---\n');
    
    for k = 1:length(labels)
        var_out = labels{k};      
        f_val = freqs(k);         
        
        f_str = num2str(f_val); 
        f_us  = strrep(f_str, '.', '_'); 
        possible_names = {var_out, sprintf('E%sHz', f_us), sprintf('E%sHz', f_str), ...
                          sprintf('Freq_%s', f_us), sprintf('f%sHz', f_us)};
        possible_names = unique(possible_names, 'stable');
        
        mat_raw = [];
        
        % Looking for the variable
        for j = 1:length(possible_names)
            if isfield(raw_file, possible_names{j})
                mat_raw = raw_file.(possible_names{j});
                break; 
            end
        end
        
        if ~isempty(mat_raw)
        
            if size(mat_raw, 2) >= 5
                mask_zeros = all(mat_raw(:, 2:5) > 0, 2);
            else
                mask_zeros = all(mat_raw(:, 2:4) > 0, 2);
            end
            mat_raw = mat_raw(mask_zeros, :);
        
            dist_vec = mat_raw(:, 1);
            
            if is_data_processed
                % TRUE: Data already processed
                E_final = mat_raw(:, 2:4);
                fprintf('  -> %s: Loaded (Data assumed as already processed).\n', var_out);
            else
                % FALSE: Raw data, we need to process (normalization and geom. spread.)

                E_raw = mat_raw(:, 2:4);
                Ecoda = ones(size(E_raw,1), 1);
                if size(mat_raw, 2) >= 5, Ecoda = mat_raw(:, 5); end
                Ecoda(Ecoda == 0) = eps; 
                
                R = E_raw ./ Ecoda;
                geom_corr = 4 * pi * (dist_vec .^ 2);
                E_final = log10( max(R .* geom_corr, eps) );
                
                fprintf('  -> %s: Processed (Coda norm + Geom spread + Log10 applied).\n', var_out);
            end
            
            data_out.(var_out) = [dist_vec, E_final];
        else
            warning('Data for frequency %g Hz (%s) not found in the file.', f_val, var_out);
        end
    end
end
