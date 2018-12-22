function [match_ind,match_ind_m] = meisam_reject_matches_based_on_scale(...
    tol, keyVector, keyVector_m, match_ind, match_ind_m);
% Reject matches at different scales
%
% INPUTS:
%   tol                     - percentage by which scale can vary within a matched pair
%   keyVector               - matrix of feature vectors (n x 4)
%   keyVector_m             - matrix of mirrored feature vectors (n x 4)
%   match_ind, match_ind_m :-
%       indices of matched pairs, i.e. keyDescriptor(match_ind) matches
%       keyDescriptor_m(match_ind_m)
%
% OUTPUTS:
%   match_ind, match_ind_m :-
%       indices of matched pairs, with matches at across scales removed.
%
scale = keyVector(:,3);
scale_m = keyVector_m(:,3);
scale_sim = (abs(scale(match_ind) - scale_m(match_ind_m))) ./ ...
    max([scale(match_ind), scale_m(match_ind_m)]')' < tol;
match_ind = match_ind(scale_sim);
match_ind_m = match_ind_m(scale_sim);

end