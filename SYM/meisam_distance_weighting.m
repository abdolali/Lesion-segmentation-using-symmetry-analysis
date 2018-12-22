
function dist_weight = meisam_distance_weighting(keyVector, keyVector_m, ...
    match_ind, match_ind_m)
% Compute optional distance weighting.  This function is not used at
% present.
%
% INPUTS:
%   keyVector               - matrix of feature vectors (n x 4)
%   keyVector_m             - matrix of mirrored feature vectors (n x 4)
%   match_ind, match_ind_m :-
%       indices of matched pairs, i.e. keyDescriptor(match_ind) matches
%       keyDescriptor_m(match_ind_m)
%
% OUTPUT:
%  dist_weight = distance weighting
%
dist_sq = (keyVector(match_ind,1) - keyVector_m(match_ind_m,1)).^2 + ...
    (keyVector(match_ind,2) - keyVector_m(match_ind_m,2)).^2;
sigma = sqrt(max(dist_sq))/6;
dist_weight = 1/(sigma*sqrt(2*pi)) * exp(-dist_sq/(2*sigma^2));

end