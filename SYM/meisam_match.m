function [match_ind,match_ind_m] = meisam_match(keyDescriptor, num_matches_per_feature, keyVector, keyDescriptor_m, keyVector_m)
% Match keypoint pairs.
%
% INPUTS:
%   keyDescriptor           - matrix of descriptor vectors (n x 128)
%   num_matches_per_feature - number of matches per feature
%   keyVector               - matrix of feature vectors (n x 4)
%   keyDescriptor_m         - matrix of mirrored descriptor vectors (n x 128)
%   keyVector_m             - matrix of mirrored feature vectors (n x 4)
%
% OUTPUTS:
%   match_ind, match_ind_m :-
%       indices of matched pairs, i.e. keyDescriptor(match_ind) matches
%       keyDescriptor_m(match_ind_m)
%
% Only need to compute upper triangular portion of distance matrix between
% the point sets since points cannot match with themselves and we don't
% want double matches. This assumes the kps in keyDescriptor and keyDescriptor_m are aligned,
% i.e., keyDescriptor_m(i,:) is the mirrored version of keyDescriptor(i,:)
%
if nargin < 5
    keyDescriptor_m = keyDescriptor;
    keyVector_m = keyVector;
end

% Build distance matrix describing distances between mirrored and
% unmirrored keypoints.
sz = size(keyDescriptor,1);
dist_mx = zeros(sz);
for k1 = 1:size(keyDescriptor,1)

    % compare a single keypoint with all other keypoints
    diff = ones(size(keyDescriptor_m,1)-k1,1) * keyDescriptor(k1,:) - keyDescriptor_m(k1+1:end,:);
    dist_mx(k1,k1+1:end) = sum(diff.^2,2)';
end

% make dist_mx symmetric
dist_mx = dist_mx + dist_mx';

dist_mx(dist_mx==0) = max(dist_mx(:));

match_ind = [];
match_ind_m = [];
for i=1:num_matches_per_feature
    match_ind = [match_ind, 1:size(keyDescriptor,1)];
    [temp, match_ind_m_new] = min(dist_mx,[],2);
%     dist_mx(:,match_ind_m_new) = 1e4;               original (meisam edited)
    dist_mx(:,match_ind_m_new) = 1e1;
    match_ind_m = [match_ind_m; match_ind_m_new];
end
end
