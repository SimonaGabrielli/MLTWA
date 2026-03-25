function do_local_fisher_plot(res, levels, title_text)
    % Plot of errors

    B_GRID = res.b0_grid; 
    L_GRID = res.le_grid; 
    [BB, LL] = meshgrid(B_GRID, L_GRID);
    norm_err = res.normalized_errors;
    
    % Colormap custom (Yellow, Orange, White)
    colors = [1 1 0; 1 0.5 0; 1 1 1];
    x_c = linspace(0, 10, 3);
    cmap = interp1(x_c, colors, linspace(0,10,256));
    
    contourf(BB, LL, norm_err, levels, 'LineColor', 'none');
    colormap(gca, cmap); hold on;
    contour(BB, LL, norm_err, levels, 'LineColor', 'k');
    
    % Minimum and Centroid
    scatter(res.x_min_absolute(1), res.x_min_absolute(2), 150, 'r', 'p',  'filled', 'MarkerEdgeColor', 'k'); 
    scatter(res.x_opt(1), res.x_opt(2), 150, 'b', 'p', 'filled', 'MarkerEdgeColor', 'k'); 
    
    title({title_text, ...
       sprintf('Min: B0=%.3f Le=%.2f | Centr: B0=%.3f Le=%.2f', ...
       res.x_min_absolute(1), res.x_min_absolute(2), res.x_opt(1), res.x_opt(2))}, ...
       'FontSize', 12);
    
    set(gca, 'FontSize', 14); 
    xlabel('B_0', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('L_e^{-1} [km^{-1}]', 'FontSize', 14, 'FontWeight', 'bold');
    axis square; grid on;
end
