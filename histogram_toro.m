function histogram_toro(p, dists_input, filename_in, filename_out)
    % HISTOGRAM TOROIDAL
    % Reads raw data and creates the space-time histogram using toroidal geometry.
    
    % --- PARAMETERS ---
    % Extract necessary parameters from struct 'p'
    tmax  = p.tmax;
    Npart = p.Npart;       % Integer
    
    % Fixed internal parameters
    Dint   = 1.0; 
    Pi     = pi;
    
    dists = dists_input; 
    Nd = numel(dists);
    
    % Time bins
    Times = 0:tmax; 
    Nt = fix(tmax); 
    
    % Calculate normalization factor (Toroidal)
    ConstantFactor = 2 * Npart * Dint^2 * Pi^2;
    ToroNorm = ConstantFactor * dists;
    
    EnergyTime = zeros(Nd, Nt, 'single');
    
    try
        M = load(filename_in); 
        % M = M;        
    catch
        error("Unable to load data file '%s'.", filename_in);
    end
    
    if isempty(M)
        warning('Data file is empty. Generating an empty histogram file.');
        dlmwrite(filename_out, [], 'delimiter', ' ');
        return;
    end

    % Column extraction (vectorized)
    % M(:,1)=x, M(:,2)=y, M(:,3)=z -> r = sqrt(x^2+y^2+z^2)
    r = sqrt(M(:,1).^2 + M(:,2).^2 + M(:,3).^2); 
    t = M(:,4); 
    e = M(:,5); 
    
    for nd = 1:Nd
        Rmin = dists(nd) - Dint;
        Rmax = dists(nd) + Dint;
        
        for nt = 1:Nt
            Tmin = Times(nt);     
            Tmax = Times(nt+1);   
            
            idx = (r >= Rmin) & (r < Rmax) & (t >= Tmin) & (t < Tmax);
            
            % Sum energies in the bin
            count = sum(e(idx)); 
            
            % Normalization
            EnergyTime(nd, nt) = count / ToroNorm(nd);
        end
    end
    
    % Save TRANSPOSED (Times on rows, Distances on columns)

    dlmwrite(filename_out, EnergyTime', 'delimiter', ' ', 'precision', 6);
    disp(['Histogram saved to: ', filename_out]);
end