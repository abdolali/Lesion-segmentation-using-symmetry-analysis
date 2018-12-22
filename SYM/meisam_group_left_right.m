function [left_ind,right_ind,pts,pts_m] = meisam_group_left_right(...
    keyVector,keyVector_m,match_ind,match_ind_m,r,ang,max_r,max_ang,...
    grouping_thresh,ind,i)
% Group left and right points of symmetric feature pairs into two distinct
% sets.
%
% INPUTS:
%   keyVector,keyVector_m,match_ind,match_ind_m,r,ang,max_r,max_ang,
%   grouping_thresh
%
% OUTPUTS:
%   ind,left_ind,right_ind,pts,pts_m
%

pts = [keyVector(match_ind(ind),2), keyVector(match_ind(ind),1)];
pts_m = [keyVector_m(match_ind_m(ind),2), keyVector_m(match_ind_m(ind),1)];

% divide matched points into two sets about symmetry axis
norm_vec = [cos(max_ang) sin(max_ang)];
left_ind = find(((pts - pts_m) * norm_vec')>0);
right_ind = find(((pts - pts_m) * norm_vec')<0);
end