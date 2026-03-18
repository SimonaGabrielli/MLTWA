function n_reflect = scattering_layer(p, filename_in, filename_out)
    % SCATTERING_FUNCTION_V2 (Multiple Layers)
    % Version with correct physics: v and g change with depth z.
    
    % --- PARAMETER EXTRACTION ---
    etai = p.etai; dt = p.dt; tmax = p.tmax;
    xs = p.xs; ys = p.ys; zs = p.zs;
    deltar = p.deltar; Npart = p.Npart;
    
    A0 = p.A0; A1 = p.A1; A2 = p.A2;
    g0 = p.g0; g1 = p.g1; g2 = p.g2;
    aa1 = p.aa1; aa2 = p.aa2;
    b1 = p.b1; b2 = p.b2;
    
    fprintf('V2 (Multiple Layers): Generating theoretical model plot...\n');
    generate_theoretical_plot(p); 

    % --- DATA LOADING ---
    try
        raw_data = load(filename_in); 
        casual = raw_data; 
    catch
        error(['Unable to read input file: ', filename_in]);
    end
    Ncasual = length(casual);
    
    pi_single = pi;
    one_single = 1;
    two_single = 2;
    
    results_cell = cell(Npart, 1);
    
    parfor j = 1:Npart
        totaltime = 0; totalspace = 0;
        xx = xs; yy = ys; zz = zs; 
        particle_history = []; 
        
        idx_rand = randi(Ncasual); 
        U1 = casual(idx_rand); U2 = casual(randi(Ncasual));
        theta = acos(one_single - two_single * U1);
        phi = two_single * pi_single * U2;
        
       while totaltime <= tmax
            
            v_val = vel(A0, A1, A2, aa1, aa2, b1, b2, zz);
            g_val = g(g0, g1, g2, aa1, aa2, b1, b2, zz);
            dvdz_val = dvdz(A0, A1, A2, aa1, aa2, b1, b2, zz);

            U3 = casual(randi(Ncasual));
            if (g_val * v_val * dt > U3)
                U1 = casual(randi(Ncasual)); 
                U2 = casual(randi(Ncasual));
                theta = acos(one_single - two_single * U1);
                phi = two_single * pi_single * U2;
            end
            
            delta_theta = dvdz_val * dt * sin(theta); 
            delta_x_prime = v_val * dt * sin(theta + delta_theta);
            delta_z = v_val * dt * cos(delta_theta + theta);
            
            dx = delta_x_prime * cos(phi);
            dy = delta_x_prime * sin(phi);
            dz = delta_z;
            
            xx_new = xx + dx;
            yy_new = yy + dy;
            zz_new = zz + dz; 
            
            totaltime = totaltime + dt;
            totalspace = totalspace + v_val * dt; 
            
            % --- DETECTION AND REFLECTION ---
            % If the particle attempted to cross the surface (z < 0)
            % Or if it ended up inside the collection volume (z < deltar)
            
            if zz_new <= 0 
                % REFLECTION OCCURRED!
                % This counts as a receiver arrival.
                
                energy = exp(-etai * totalspace);
                
                new_row = [xx_new, yy_new, 0, totaltime, energy];
                particle_history = [particle_history; new_row]; 
                
                zz_new = -zz_new;  
                theta = pi_single - theta; 
                
                xx = xx_new; yy = yy_new; zz = zz_new;

            elseif zz_new <= deltar
                % Inside the surface volume but did not touch zero (rare case with large dt)
                % We record this as well for safety.
                energy = exp(-etai * totalspace);
                new_row = [xx_new, yy_new, zz_new, totaltime, energy];
                particle_history = [particle_history; new_row];
                
                xx = xx_new; yy = yy_new; zz = zz_new;
                
            else
                xx = xx_new; yy = yy_new; zz = zz_new;
            end
        
        end
        if ~isempty(particle_history), results_cell{j} = particle_history; end
    end

    n_reflect = write_output_file(filename_out, results_cell);
end


function val = vel(A0, A1, A2, aa1, aa2, b1, b2, z)
    term1 = (A1 - A0) / (1 + exp(-b1 * (z - aa1)));
    term2 = (A2 - A1) / (1 + exp(-b2 * (z - aa2)));
    val = A0 + term1 + term2;
end

function val = dvdz(A0, A1, A2, aa1, aa2, b1, b2, z)
    sig1 = 1 / (1 + exp(-b1 * (z - aa1)));
    sig2 = 1 / (1 + exp(-b2 * (z - aa2)));
    dv1 = (A1 - A0) * sig1 * (1 - sig1) * b1;
    dv2 = (A2 - A1) * sig2 * (1 - sig2) * b2;
    val = dv1 + dv2;
end

function val = g(g0, g1, g2, aa1, aa2, b1, b2, z)
    term1 = (g1 - g0) / (1 + exp(-b1 * (z - aa1)));
    term2 = (g2 - g1) / (1 + exp(-b2 * (z - aa2)));
    val = g0 + term1 + term2;
end

function generate_theoretical_plot(p)
    z_plot = linspace(0, 60, 1000); 
    v_curve = arrayfun(@(z) vel(p.A0, p.A1, p.A2, p.aa1, p.aa2, p.b1, p.b2, z), z_plot);
    g_curve = arrayfun(@(z) g(p.g0, p.g1, p.g2, p.aa1, p.aa2, p.b1, p.b2, z), z_plot);
    % dv_curve = arrayfun(@(z) dvdz(p.A0, p.A1, p.A2, p.aa1, p.aa2, p.b1, p.b2, z), z_plot);

    %fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 600]);

    yyaxis left; 
    h1 =plot(z_plot, g_curve, 'b-', 'LineWidth', 2);
    ylabel('Scattering g'); ylim([0, max(g_curve)*1.1]); ax=gca; ax.YColor='b';
    yyaxis right; 
    h2 = plot(z_plot, v_curve, 'r--', 'LineWidth', 2); 
    hold on;

    % plot(z_plot, dv_curve, 'k-', 'LineWidth', 1.5);
    ylabel('Velocity v'); ylim([0, max(max(v_curve))*1.2]); ax=gca; ax.YColor='r';
        % ylabel('Velocity v & dv/dz'); ylim([0, max(max(v_curve), max(dv_curve))*1.2]); ax=gca; ax.YColor='k';

    grid on; 
    
    title('Theoretical Model (Dynamic Version)'); 
    % legend({'g','v','dv/dz'},'Location','best');
    legend([h1, h2],{'g','v'},'Location','best');

    % exportgraphics(fig, 'Theoretical_Model.png', 'Resolution', 300);
    % close(fig);

    % fprintf('Model Visualization:\n'); imshow('Theoretical_Model_V2.png'); axis image; drawnow; fprintf('\n');
end

function n = write_output_file(fname, results)
    results = results(~cellfun('isempty', results));
    fid = fopen(fname, 'w');
    if isempty(results)
        n = 0;
    else
        final_mat = vertcat(results{:});
        n = size(final_mat, 1);
        fprintf(fid, '%.6f %.6f %.6f %.6f %.6e\n', final_mat');
    end
    fclose(fid);
end