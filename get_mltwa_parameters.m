function [etai, g0] = get_mltwa_parameters(resultsTable, target_f)
    % Trova l'indice della frequenza
    idx = find(resultsTable.Freq_Hz == target_f);
    if isempty(idx), error('Frequenza %.1f non trovata nella tabella.', target_f); end
    
    % Estrazione stringhe dalle celle
    B0_str = string(resultsTable.B0{idx});
    Le_inv_str = string(resultsTable.Le_inv{idx});
    
    % Usiamo 'strsplit' per prendere solo la prima parte prima dello spazio
    parts_B0 = strsplit(B0_str, ' ');
    parts_Le = strsplit(Le_inv_str, ' ');
    
    B0_val = str2double(parts_B0(1));
    Le_inv_val = str2double(parts_Le(1));
    
    % Verifichiamo che la conversione sia riuscita
    if isnan(B0_val) || isnan(Le_inv_val)
        error('Errore nella conversione numerica: controlla il formato della tabella.');
    end
    
    % ricalcolo etai e g0 (eta_s)
    g0 = Le_inv_val * B0_val;
    etai = Le_inv_val * (1 - B0_val);
    
    fprintf('Bridge Success: f=%.1f, Le^-1=%.3f, B0=%.2f -> g0=%.4f, etai=%.4f\n', ...
            target_f, Le_inv_val, B0_val, g0, etai);
end