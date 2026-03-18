function histogram_strip(p, dists_input, filename_in, filename_out)
    % HISTOGRAM_STRIP_CREATION
    % Crea l'istogramma usando la normalizzazione "Strip" (Sferica)
    % replicando la logica del codice Fortran/Snippet fornito.
    % Mantiene la sintassi compatibile con le funzioni precedenti.
    
    % --- 1. ESTRAZIONE PARAMETRI DA 'p' ---
    dt_sim = p.dt;      % Passo di simulazione (es. 0.05, 0.1, 0.2)
    tmax   = p.tmax;
    Npart  = p.Npart;
    
    % Recuperiamo la velocità superficiale dalla struttura, se esiste.
    % Altrimenti usiamo il valore del tuo snippet (2.45) o del paper (3.5).
    if isfield(p, 'A0')
        vel = p.A0;
    else
        vel = 2.45; % Valore di default dal tuo snippet (ATTENZIONE: Verifica!)
        warning('Velocità A0 non trovata in p. Uso valore default: %.2f', vel);
    end
    
    % --- 2. PARAMETRI SPECIFICI STRIP ---
    Dint   = 1.0;          % Semilarghezza del guscio (Delta r)
    TLint  = 4.0 * dt_sim; % Larghezza Bin Temporale (Time Bin Width) - Dal tuo snippet
    
    % --- 3. DEFINIZIONE BIN ---
    % Vettore dei tempi (Bordi dei bin)
    Edges = 0 : TLint : tmax;
    Nt = length(Edges) - 1;
    
    % Vettore delle distanze
    dists = dists_input;
    % Rimuoviamo lo zero se presente per evitare divisioni o errori logici
    if any(dists == 0), dists(dists == 0) = []; end
    Nd = numel(dists);
    
    % --- 4. CALCOLO NORMALIZZAZIONE (STRIP NORM) ---
    % Formula dal tuo snippet: 
    % StripNorm = 4 * Pi * TLint * Npart * Dist * (Dint + vel*dt)^2
    
    thickness_term = (Dint + (vel * dt_sim))^2;
    NormScalar = 4 * pi * TLint * Npart * thickness_term;
    
    % Vettore normalizzazione per ogni distanza (scalare * distanza)
    StripNormVec = NormScalar * dists; 
    
    % --- 5. CARICAMENTO DATI ---
    try
        M = load(filename_in);
    catch
        warning('File %s non trovato o vuoto. Creo output vuoto.', filename_in);
        dlmwrite(filename_out, [], 'delimiter', ' ');
        return;
    end
    
    if isempty(M)
        dlmwrite(filename_out, [], 'delimiter', ' ');
        return;
    end
    
    % Coordinate e Energia
    % r = distanza radiale (sferica)
    r = sqrt(M(:,1).^2 + M(:,2).^2 + M(:,3).^2);
    t = M(:,4);
    e = M(:,5);
    
    % Preallocazione matrice output
    EnergyTime = zeros(Nd, Nt);
    
    % --- 6. CICLO DI CALCOLO VELOCE ---
    for nd = 1:Nd
        rmin = dists(nd) - Dint;
        rmax = dists(nd) + Dint;
        
        % Filtro Spaziale (Particelle nel guscio sferico)
        idx_space = (r >= rmin) & (r < rmax);
        
        if any(idx_space)
            t_subset = t(idx_space);
            e_subset = e(idx_space);
            
            % Binning Temporale Veloce (sostituisce il ciclo 'for nt')
            % histc restituisce l'indice del bin per ogni particella
            [~, bin_idx] = histc(t_subset, Edges);
            
            % Filtriamo indici validi (1 a Nt)
            valid = (bin_idx > 0) & (bin_idx <= Nt);
            
            if any(valid)
                % accumarray somma l'energia (e_subset) raggruppandola per bin
                % Questo è equivalente al ciclo for temporale ma istantaneo
                counts = accumarray(bin_idx(valid), e_subset(valid), [Nt, 1]);
                
                % Normalizzazione Strip
                EnergyTime(nd, :) = counts / StripNormVec(nd);
            end
        end
    end
    
    % --- 7. SALVATAGGIO ---
    % Trasposta: Righe=Tempo, Colonne=Distanze
    dlmwrite(filename_out, EnergyTime', 'delimiter', ' ', 'precision', 6);
    fprintf('Istogramma Strip salvato in: %s (vel=%.2f, TLint=%.3f)\n', ...
            filename_out, vel, TLint);
end