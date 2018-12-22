function ind = meisam_assign_to_axis(r,ang,max_r,max_ang,grouping_thresh)
% Assign pairs to the dominant axis of symmetry. Tag all pairs that agree
% with the dominant axis, and determine extent of axis in image.
%
% INPUTS:
%   r, ang          - Hough space polar coordinates of axes of symmetry
%   max_r, max_ang  - Hough space polar coordinates of dominant axis of symmetry
%   grouping_thresh - angular and radial tolerances for sym particles
%                     associated with the same symmetry axis.
%
% OUTPUT:
%   ind - indices of pairs assigned to dominant axis of symmetry

% compute angular (d_ang1) radial (d_r1) difference between each pair and
% the dominant axis.
temp = unwrap( [ang, repmat(max_ang,length(ang),1)], [], 2);
d_ang1 = abs(temp(:,1) - temp(:,2));
d_r1 = abs(r - max_r);

% Do the same again, but with the axis oriented 180 degrees in the
% opposite directioncompute angular (d_ang2) radial (d_r2) difference
% between each pair and the dominant axis.
temp = unwrap( [ang+pi, repmat(max_ang,length(ang),1)], [], 2);
d_ang2 = abs(temp(:,1) - temp(:,2));
d_r2 = abs(-r - max_r);

% choose the result that is closest to the dominant axis
[d_ang iind] = min([d_ang1, d_ang2],[],2);
d_r_comb = [d_r1, d_r2];
d_r(iind==1) = d_r1(iind==1);
d_r(iind==2) = d_r2(iind==2);

% Assign pairs that are sufficiently close to the dominant axis to this
% axis.
ind = find((d_ang<grouping_thresh(1)) & (d_r'<grouping_thresh(2)));

end