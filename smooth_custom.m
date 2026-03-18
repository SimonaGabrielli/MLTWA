function smoothed_list = smooth_custom(list, n, s)
% SMOOTH_CUSTOM: Moving Average
    smoothed_list = [];
    if isempty(list) || length(list) < n, return; end
    idx = 1;
    while (idx + n - 1) <= length(list)
        smoothed_list(end+1) = mean(list(idx : idx + n - 1));
        idx = idx + s;
    end
end