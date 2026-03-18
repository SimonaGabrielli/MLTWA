function plot_envelopes(filename_in, dists_vec)

% Reads the histogram file and plots the energy curves.

try
    EnergyOut = load(filename_in);
catch
    error(['Unable to load histogram file: ', filename_in]);
end

if isempty(EnergyOut)
    warning('Histogram file is empty. No plots generated.');
    return;
end

num_dists_file = size(EnergyOut, 2);
num_dists_input = length(dists_vec);

if num_dists_file ~= num_dists_input
    warning('Warning: The number of columns in the file (%d) does not match the distances vector (%d).', ...
        num_dists_file, num_dists_input);

    dists_vec = 1:num_dists_file;
end

% Assumes each row is 1 second, centered at -0.5 (e.g., row 1 = 0.5s)
Times = (1:size(EnergyOut, 1))' - 0.5;

figure('Name','Envelope Compact', 'Color','w', 'Position',[100, 100, 1000, 550]);

nplots = size(EnergyOut, 2);
ncols = 3;
nrows = ceil(nplots / ncols);

t = tiledlayout(nrows, ncols, 'Padding','compact', 'TileSpacing','compact');

for nsta = 1:nplots
    nexttile;

    y_data = EnergyOut(:, nsta);

    % --- VERTICAL LINE ---
    % This will add a vertical line from the x to the maximum energy
    % value. Not mandatory, it could be removed

    y_data(y_data <= 0) = 1e-30;

    semilogy(Times, y_data, 'r-', 'LineWidth', 1.2);

    pbaspect([2 1 1])


    xlim([0, max(Times)]);

    ylim([1e-12, 0.1]);

    grid off;

    title(sprintf('r = %.1f km', dists_vec(nsta)), 'FontSize', 10, 'FontWeight', 'bold');

    if mod(nsta-1, ncols) == 0
        ylabel('Energy');
    end
    if nsta > (nrows-1)*ncols
        xlabel('Time [s]');
    end

    % Custom Ticks - this could vary
    ax = gca;
    ax.YScale = 'log';
    ytickvals = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2];
    yticklabels_vals = ["10^{-9}", "10^{-8}", "10^{-7}", "10^{-6}", "0.00001", "0.0001", "0.001", "0.01"];

    yticks(ytickvals);
    yticklabels(yticklabels_vals);
end

title(t, 'Scattering Envelopes', 'FontSize', 14, 'FontWeight', 'bold');
end