function results = gridSearch_Core(r_vec, E1_vec, E2_vec, E3_vec, params, grid_cfg)
% gridSearch_Core: Calculates the error and returns the centroid.

b0_vec = grid_cfg.b0_vec;
le_vec = grid_cfg.le_vec;

% Grid Generation
[B_grid, L_grid] = meshgrid(b0_vec, le_vec);
errors = zeros(numel(B_grid), 1);

p_loc = params;

if isfield(params, 'num_workers')
    if strcmpi(params.num_workers, 'auto')
        num_w = Inf; % MATLAB default (uses all available cores)
    else
        num_w = round(params.num_workers); % User defined (0 means sequential)
    end
else
    num_w = Inf; % Fallback
end

% Safety: If user wants parallel but doesn't have the Parallel Toolbox installed
if num_w > 0 && exist('parpool', 'file') == 0
    fprintf('   [!] Parallel Computing Toolbox not found. Forcing sequential run (0 workers).\n');
    num_w = 0;
end

% The loop: if num_w is 0, it acts exactly like a standard 'for' loop
parfor (k = 1:numel(B_grid), num_w)
    
    x_try = [B_grid(k), L_grid(k)];

    % Energy Calculation
    Y1 = computeEnergy3(x_try, r_vec, 1, p_loc);
    Y2 = computeEnergy3(x_try, r_vec, 2, p_loc);
    Y3 = computeEnergy3(x_try, r_vec, 3, p_loc);

    % NaN removal
    if any(isnan(Y1)), Y1(isnan(Y1)) = -100; end
    if any(isnan(Y2)), Y2(isnan(Y2)) = -100; end
    if any(isnan(Y3)), Y3(isnan(Y3)) = -100; end

    % Sum of Squared Errors (SSE)
    errors(k) = sum((Y1 - E1_vec).^2) + sum((Y2 - E2_vec).^2) + sum((Y3 - E3_vec).^2);
end

results.min_error = min(errors(:));
results.error_matrix = reshape(errors, size(B_grid));
results.normalized_errors = results.error_matrix / results.min_error;

% Find the Absolute Minimum
[r_idx, c_idx] = find(results.error_matrix == results.min_error, 1);
x_absolute_min = [b0_vec(c_idx), le_vec(r_idx)];

% --- Fisher Calculation ---
n_data = length(r_vec);
p_params = 2; % B0 and Le

results.F_critic1 = finv(0.68, p_params, n_data - p_params); % Visual "1 sigma"
results.F_critic2 = finv(0.90, p_params, n_data - p_params); % Visual "90%"

% --- Local Minimum and best Solution ---
threshold_norm = results.F_critic1;

% 1. Find local minima in the matrix (neighbor scan)
[rows, cols] = size(results.error_matrix);
local_min_mask = false(rows, cols);

for r = 1:rows
    for c = 1:cols

        r_min = max(1, r-1); r_max = min(rows, r+1);
        c_min = max(1, c-1); c_max = min(cols, c+1);

        neighborhood = results.error_matrix(r_min:r_max, c_min:c_max);

        % If the point is the lowest (or equal to the lowest) in its neighborhood,
        % we classify it as a local minimum.
        if results.error_matrix(r,c) == min(neighborhood(:))
            local_min_mask(r,c) = true;
        end
    end
end


% 2. We keep only the "good" ones (below the F_critic threshold)
significant_minima = local_min_mask & (results.normalized_errors <= threshold_norm);

num_minima = sum(significant_minima(:));

if num_minima == 1
    % Case A: There is only one clear minimum. No centroid!
    results.x_opt = x_absolute_min;

elseif num_minima > 1
    % Case B: There are multiple significant local minima
    avg_B = mean(B_grid(significant_minima));
    avg_L = mean(L_grid(significant_minima));
    results.x_opt = [avg_B, avg_L];

else
    % Case C: Safety fallback (if all goes wrong, we use the absolute minimum)
    results.x_opt = x_absolute_min;
end

results.x_min_absolute = x_absolute_min;
results.b0_grid = b0_vec;
results.le_grid = le_vec;

% --- Uncertainty Calculation ---
confidence_area = results.normalized_errors <= threshold_norm;
if sum(confidence_area(:)) > 1
    results.err_B0 = (max(B_grid(confidence_area)) - min(B_grid(confidence_area))) / 2;
    results.err_Le = (max(L_grid(confidence_area)) - min(L_grid(confidence_area))) / 2;

    results.good_B = B_grid(confidence_area);
    results.good_L = L_grid(confidence_area);
else
    results.err_B0 = 0;
    results.err_Le = 0;
    results.good_B = results.x_opt(1);
    results.good_L = results.x_opt(2);
end

results.b0_grid = b0_vec;
results.le_grid = le_vec;
end