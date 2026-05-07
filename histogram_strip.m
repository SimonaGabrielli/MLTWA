function histogram_strip(p, dists_input, filename_in, filename_out)

    dt_sim = p.dt;    
    tmax   = p.tmax;
    Npart  = p.Npart;
    

    if isfield(p, 'A0')
        vel = p.A0;
    else
        vel = 2.45; 
        warning('Velocità A0 non trovata in p. Uso valore default: %.2f', vel);
    end
    
    % --- Specific parameters for the histogram ---
    Dint   = 1.0;          % Semi-width (delta r)
    TLint  = 4.0 * dt_sim; % Time Bin Width
    
    % --- Defining the Bin ---
    % Time vector (bin edges)
    Edges = 0 : TLint : tmax;
    Nt = length(Edges) - 1;
    
    % Distance vector
    dists = dists_input;
    if any(dists == 0), dists(dists == 0) = []; end
    Nd = numel(dists);
    
    % --- Normalization ---
    % StripNorm = 4 * Pi * TLint * Npart * Dist * (Dint + vel*dt)^2
    
    thickness_term = (Dint + (vel * dt_sim))^2;
    NormScalar = 4 * pi * TLint * Npart * thickness_term;
    
    StripNormVec = NormScalar * dists; 
    
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
    
    % r = radial distance (spherical)
    r = sqrt(M(:,1).^2 + M(:,2).^2 + M(:,3).^2);
    t = M(:,4);
    e = M(:,5);
    
    EnergyTime = zeros(Nd, Nt);
    
    for nd = 1:Nd
        rmin = dists(nd) - Dint;
        rmax = dists(nd) + Dint;
        
        % Particles in the sphere
        idx_space = (r >= rmin) & (r < rmax);
        
        if any(idx_space)
            t_subset = t(idx_space);
            e_subset = e(idx_space);
            
            %histc gives the index of the bin for each particle
            [~, bin_idx] = histc(t_subset, Edges);
            
            valid = (bin_idx > 0) & (bin_idx <= Nt);
            
            if any(valid)
               
                counts = accumarray(bin_idx(valid), e_subset(valid), [Nt, 1]);
                
                % Strip normalization
                EnergyTime(nd, :) = counts / StripNormVec(nd);
            end
        end
    end
    
    dlmwrite(filename_out, EnergyTime', 'delimiter', ' ', 'precision', 6);
    fprintf('Istogramma Strip salvato in: %s (vel=%.2f, TLint=%.3f)\n', ...
            filename_out, vel, TLint);
end