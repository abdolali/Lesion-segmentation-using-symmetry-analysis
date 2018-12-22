function [match_ind,match_ind_m,ang,phase_weight] = meisam_angular_constraint(...
    keyVector, keyVector_m, match_ind, match_ind_m)
% Reject matches that are not symmetrically oriented
%
% INPUTS:
%   keyVector               - matrix of feature vectors (n x 4)
%   keyVector_m             - matrix of mirrored feature vectors (n x 4)
%   match_ind, match_ind_m :- indices of matched pairs, i.e.
%                             keyDescriptor(match_ind) matches
%                             keyDescriptor_m(match_ind_m)
%
% OUTPUTS:
%   match_ind, match_ind_m :- indices of matched pairs, after removing
%                             non-symmetric matches
%   ang                     - ang of lines joining feature pairs
%   phase_weight            - phase weight contribution to symmetry
%                             magnitude (n x 1) vector.

% determine orientation of line joining the pair
ang = atan2( (keyVector(match_ind,1) - keyVector_m(match_ind_m,1)), ...
    (keyVector(match_ind,2) - keyVector_m(match_ind_m,2)) );

% flip angle of mirrored points about y-axis.
% necessary because not a right hand coordinate frame.
r = ones(size(keyVector_m(:,4)));
[uxm,uym] = pol2cart(keyVector_m(:,4), r);
uxm = -uxm;
[keyVector(:,4),r] = cart2pol(uxm,uym);
[keyVector_m(:,4),r] = cart2pol(uxm,uym);

% compute angular component of first part of Reisfeld's phase weight function [Reisfeld94]
ang_i = keyVector(match_ind,4);
ang_j = keyVector_m(match_ind_m,4);
phase = ang_i + ang_j - 2*ang;
phase_weight = - cos(phase);
phase_weight(phase_weight < 0) = 0;

% discard pairs that are assymetrically oriented
sym_orient = phase_weight > 0;
match_ind = match_ind(sym_orient);
match_ind_m = match_ind_m(sym_orient);
ang = ang(sym_orient);
phase_weight = phase_weight(sym_orient);
end