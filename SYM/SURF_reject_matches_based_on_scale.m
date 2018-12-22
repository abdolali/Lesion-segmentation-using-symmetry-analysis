function [match_ind,match_ind_m] = SURF_reject_matches_based_on_scale(...
    tol, keyVector, keyVector_m, match_ind, match_ind_m);

scale = keyVector(:,3);
scale_m = keyVector_m(:,3);
scale_sim = (abs(scale(match_ind) - scale_m(match_ind_m))) ./ ...
    max([scale(match_ind), scale_m(match_ind_m)]')' < tol;
match_ind = match_ind(scale_sim);
match_ind_m = match_ind_m(scale_sim);

end