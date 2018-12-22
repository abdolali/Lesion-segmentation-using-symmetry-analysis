function  [surfingout,segmout,sym_measure, symstrength , symangle] = SURF_symmetry(im, match_th)

addpath('siftDemoV4/');

if size(im,3)>1
    im = sum(im,3)/3;
end
im1 = uint8(double(im));

grouping_thresh = [3/180*pi, 3];
tol = 0.2;

[xsurfi, ysurfi] = size(im1);
if round(ysurfi/2) ~=ysurfi/2
    im1 = im1(:,1:end-1);
    ysurfi=ysurfi-1;
end

im2 = zeros(xsurfi, ysurfi);
for j=1:ysurfi
    im2(:,ysurfi - j+1)=im1(:,j);
end
im2 = uint8(im2);

pointsF_1 = detectSURFFeatures(im1);
pointsF_2 = detectSURFFeatures(im2);

pFp(:,1) = double(pointsF_1.Scale);
pFp(:,2) = double(pointsF_1.SignOfLaplacian);
pFp(:,3) = double(pointsF_1.Orientation);
pFp(:,4:5) = double(pointsF_1.Location);

keyVector(:,1) = pFp(:,5) ; 
keyVector(:,2) = pFp(:,4) ; 
keyVector(:,3) = double(pointsF_1.Scale);
keyVector(:,4) = double(pointsF_1.Orientation);

pFp2(:,1) = double(pointsF_2.Scale);
pFp2(:,2) = double(pointsF_2.SignOfLaplacian);
pFp2(:,3) = double(pointsF_2.Orientation);
pFp2(:,4:5) = double(pointsF_2.Location);

keyVector_m(:,1) = keyVector(:,1);
keyVector_m(:,2) = keyVector(:,2);
keyVector_m(:,3) = keyVector(:,3);
keyVector_m(:,4) = keyVector(:,4); 


% [im, keyDescriptor, keyVector] = siftV4_mod(im);

if (~pointsF_1.Count || ~pointsF_2.Count)
%     output=0;
        surfingout=[0;0]; segmout=0; sym_measure=0; symstrength =0; symangle=0;
    return 
end

%%
% correct keypoint locatoins to be 1 pixel to the right and lower
% keyVector(:,[1,2]) = keyVector(:,[1,2]) + 1;

% directly modify SIFT feature vector to represent reflected regions
[fF1,vpF1] = extractFeatures(im1,pointsF_1);
[fF2,vpF2] = extractFeatures(im2,pointsF_2);

keyDescriptor = fF1 ;
keyDescriptor_m = fF2 ;

% keyDescriptor_m = meisam_mirror_SIFT_descriptor(keyDescriptor);
% keyVector_m = keyVector;



% REFLECTIVE symmetry
% num_matches_per_feature = 1;

ipF = matchFeatures(fF1,fF2,'MatchThreshold', match_th);
match_ind = ipF(:,1) ;
match_ind_m = ipF(:,2) ;

% [match_ind,match_ind_m] = meisam_match(keyDescriptor, num_matches_per_feature, keyVector, keyDescriptor_m, keyVector_m);
[match_ind,match_ind_m] = SURF_reject_matches_based_on_scale(tol, keyVector,keyVector_m,match_ind,match_ind_m);
[match_ind,match_ind_m,ang,phase_weight] = SURF_angular_constraint(keyVector, keyVector_m, match_ind, match_ind_m);
 % Compute optional distance weighting:
 %   dist_weight = meisam_distance_weighting(keyVector, keyVector_m, match_ind, match_ind_m);
if isempty(match_ind)
     surfingout=[0;0]; segmout=0; sym_measure=0; symstrength =0; symangle=0;
     return 
end
           
        % identify sym axes particles
        sym_x = mean([keyVector(match_ind,2), keyVector_m(match_ind_m,2)],2);
        sym_y = mean([keyVector(match_ind,1), keyVector_m(match_ind_m,1)],2);
        [Hough_im, max_r, max_ang, r, sym_strength,mhmax] = linear_hough(im,sym_x,sym_y,ang,phase_weight);

        for i=1:length(max_r) %min(3,length(max_r))
            ind{i} = meisam_assign_to_axis(r,ang,max_r(i),max_ang(i),grouping_thresh);
            [left_ind{i},right_ind{i},pts{i},pts_m{i}] = meisam_group_left_right( keyVector,keyVector_m,match_ind, match_ind_m,r,ang,max_r,max_ang(i),grouping_thresh,ind{i});
        end

        [surfingout,segmout]=meisam_display_output(im,keyVector,ang,sym_x,sym_y, ind,match_ind,match_ind_m,Hough_im,keyVector_m, max_ang,max_r,left_ind,right_ind,pts, pts_m,phase_weight,sym_strength);
%         output = [max_r', max_ang', sym_strength'];
          symangle = max_ang(1);
% % %          output (2)= OuT;
         symstrength=sym_strength(1);
%%%%mh 304 symmetry measure:
         sym_measure = mhmax/size(match_ind,2);
%%%%%
%             output=OuT;           %%%% for surfacing
end