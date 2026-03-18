function ResultsTable = extract_mltwa_results_table(final_results, freqs, labels, v, use_centroid, use_smooth)
% extract_mltwa_results_table: Formats results into a clean table using inverse Q (Q^-1).

Freq_Col   = [];
B0_Col_Str = strings(0, 1);
Le_Col_Str = strings(0, 1);
Qt_Col_Str = strings(0, 1);
Qs_Col_Str = strings(0, 1);
Qi_Col_Str = strings(0, 1);

for i = 1:length(labels)
    lbl = labels{i};
    f   = freqs(i); 
    
    if isfield(final_results, lbl)
        res_freq = final_results.(lbl);
        
        if use_smooth && isfield(res_freq, 'SMOOTH')
            res_inner = res_freq.SMOOTH;
        elseif isfield(res_freq, 'RAW')
            res_inner = res_freq.RAW;
        else
            res_inner = res_freq;
        end
        
        if use_centroid && isfield(res_inner, 'x_opt')
            vals = res_inner.x_opt;
            B0_val = vals(1); Le_val = vals(2);
        elseif ~use_centroid && isfield(res_inner, 'x_min_absolute')
            vals = res_inner.x_min_absolute;
            B0_val = vals(1); Le_val = vals(2);
        else
            B0_val = NaN; Le_val = NaN;
        end
        
        if isfield(res_inner, 'err_B0'), err_B = res_inner.err_B0; else, err_B = 0; end
        if isfield(res_inner, 'err_Le'), err_L = res_inner.err_Le; else, err_L = 0; end
        
        if ~isnan(B0_val) && Le_val > 0
            % --- MODIFICA: Calcolo dei valori inversi (Q^-1) ---
            % Qt_inv = (Le * v) / (2 * pi * f)
            Qt_inv = (Le_val * v) / (2 * pi * f);
            Qs_inv = Qt_inv * B0_val;
            Qi_inv = Qt_inv * (1 - B0_val);
            
            if isfield(res_inner, 'good_B') && length(res_inner.good_B) > 1
                g_B = res_inner.good_B;
                g_L = res_inner.good_L;
                
                % Calcolo errori sui valori inversi
                all_Qt_inv = (g_L * v) ./ (2 * pi * f);
                all_Qs_inv = all_Qt_inv .* g_B;
                all_Qi_inv = all_Qt_inv .* (1 - g_B);
                
                err_Qt_inv = (prctile(all_Qt_inv, 84) - prctile(all_Qt_inv, 16)) / 2;
                err_Qs_inv = (prctile(all_Qs_inv, 84) - prctile(all_Qs_inv, 16)) / 2;
                err_Qi_inv = (prctile(all_Qi_inv, 84) - prctile(all_Qi_inv, 16)) / 2;
            else
                err_Qt_inv = 0; err_Qs_inv = 0; err_Qi_inv = 0;
            end
            
            Freq_Col = [Freq_Col; f];
            str_B0 = sprintf('%.2f ± %.2f', B0_val, err_B);
            str_Le = sprintf('%.3f ± %.3f', Le_val, err_L);
            
            % --- MODIFICA: Formattazione scientifica per Qt, Qi, Qs (inversi) ---
            str_Qt = sprintf('%.2e ± %.2e', Qt_inv, err_Qt_inv);
            str_Qs = sprintf('%.2e ± %.2e', Qs_inv, err_Qs_inv);
            str_Qi = sprintf('%.2e ± %.2e', Qi_inv, err_Qi_inv);
            
            B0_Col_Str = [B0_Col_Str; string(str_B0)];
            Le_Col_Str = [Le_Col_Str; string(str_Le)];
            Qt_Col_Str = [Qt_Col_Str; string(str_Qt)];
            Qs_Col_Str = [Qs_Col_Str; string(str_Qs)];
            Qi_Col_Str = [Qi_Col_Str; string(str_Qi)];
        end
    end
end

ResultsTable = table(Freq_Col, Le_Col_Str, B0_Col_Str, ...
    Qt_Col_Str, Qi_Col_Str, Qs_Col_Str, ...
    'VariableNames', {'Freq_Hz', 'Le_inv', 'B0', 'Qt_inv', 'Qi_inv', 'Qs_inv'});
end