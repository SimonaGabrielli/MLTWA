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
    scatter(res.x_min_absolute(1), res.x_min_absolute(2), 50, 'r', 'x', 'LineWidth', 2); 
    scatter(res.x_opt(1), res.x_opt(2), 120, 'b', 'p', 'filled', 'MarkerEdgeColor', 'k'); 
    
    title({title_text, sprintf('Centr: B0=%.3f Le=%.3f', res.x_opt(1), res.x_opt(2))}, 'FontSize', 12);
    xlabel('Scattering B_0'); ylabel('Inv. Extinction Le^{-1}');
    axis square; grid on;
end