function [r_sm, E1_sm, E2_sm, E3_sm] = prepare_smoothed_data(r, E1, E2, E3, n_win, s_win)
    
    % Application of smoothing
    r_sm  = smooth_custom(r, n_win, s_win)';
    E1_sm = smooth_custom(E1, n_win, s_win)';
    E2_sm = smooth_custom(E2, n_win, s_win)';
    E3_sm = smooth_custom(E3, n_win, s_win)';
    
    r_sm = r_sm(:);
    E1_sm = E1_sm(:);
    E2_sm = E2_sm(:);
    E3_sm = E3_sm(:);
end