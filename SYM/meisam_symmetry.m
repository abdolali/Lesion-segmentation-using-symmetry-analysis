function  [surfingout,segmout,sym_measure, symstrength , symangle] = meisam_symmetry(im);

addpath('siftDemoV4/');

if size(im,3)>1
    im = sum(im,3)/3;
end
im = uint8(double(im));

% angular and radial tolerances for sym particles associated with the same symmetry axis
grouping_thresh = [3/180*pi, 3];
% percentage by which scale can vary within a matched pair
tol = 0.2;

% ----------------------------------
% COMPUTE KEYPOINTS
% ----------------------------------
[im, keyDescriptor, keyVector] = siftV4_mod(im);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~, keyDescriptor_m, keyVector_m] = siftV4_mod(im_m);
% keyVector_m = keyVector_m (1:end-2,:);
% keyDescriptor_m = keyDescriptor_m (1:end-2,:);
% keyVector_m(:,[1,2]) = keyVector_m(:,[1,2]) + 1;
% %%%%  mh:
% im =imread('E:\2.jpg');
% [im_x, im_y] = size(im);
% im_m = zeros(size(im));
% for j=1:im_y
%     im_m(:,im_y-j+1)=im(:,j);
% end
% im_m = uint8(im_m);
% points = detectBRISKFeatures(im);
% points = detectSURFFeatures(im);
% [features,validPoints] = extractFeatures(im,points);
% points_m = detectSURFFeatures(im_m);
% points_m = detectBRISKFeatures(im_m);
% 
% [features_m,validPoints_m] = extractFeatures(im_m,points_m);
% [indexPairs,matchmetric] = matchFeatures(features,features_m,'MatchThreshold',100);
% matched_pts1 = validPoints(indexPairs(:, 1));
% matched_pts2 = validPoints_m(indexPairs(:, 2));
% figure; showMatchedFeatures(im, im_m, matched_pts1, matched_pts2);
% title('Putative point matches');
% legend('matchedPts1','matchedPts2');
% 
% 
% figure;imshow(im); hold on;plot(points);
% figure;imshow(im_m); hold on;plot(points_m);
% keyVector(:,1:2) = validPoints.Location;
% keyVector(:,3)  = validPoints.Scale;
% keyVector(:,4) = validPoints.Orientation;
% keyDescriptor(:,1:64) = features;
% keyDescriptor(:,65:128) = features;

%%

if isempty(keyVector)
%     output=0;
        surfingout=[0;0]; segmout=0; sym_measure=0; symstrength =0; symangle=0;
    return 
end

%%
% correct keypoint locatoins to be 1 pixel to the right and lower
keyVector(:,[1,2]) = keyVector(:,[1,2]) + 1;

% directly modify SIFT feature vector to represent reflected regions
keyDescriptor_m = meisam_mirror_SIFT_descriptor(keyDescriptor);
keyVector_m = keyVector;

% REFLECTIVE symmetry
num_matches_per_feature = 1;
[match_ind,match_ind_m] = meisam_match(keyDescriptor, num_matches_per_feature, keyVector, keyDescriptor_m, keyVector_m);
[match_ind,match_ind_m] = meisam_reject_matches_based_on_scale(tol, keyVector,keyVector_m,match_ind,match_ind_m);
[match_ind,match_ind_m,ang,phase_weight] = meisam_angular_constraint(keyVector, keyVector_m, match_ind, match_ind_m);
 % Compute optional distance weighting:
 %   dist_weight = meisam_distance_weighting(keyVector, keyVector_m, match_ind, match_ind_m);
if isempty(match_ind)
     surfingout=[0;0]; segmout=0; sym_measure=0; symstrength =0; symangle=0;
     return 
end
           
        % identify sym axes particles
        sym_x = mean([keyVector(match_ind,2), keyVector_m(match_ind_m,2)],2);
        sym_y = mean([keyVector(match_ind,1), keyVector_m(match_ind_m,1)],2);
        [Hough_im, max_r, max_ang, r, sym_strength,mhmax, mhStrength] = linear_hough(im,sym_x,sym_y,ang,phase_weight);

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
            sym_measure = mhStrength ;
%          sym_measure = mhmax/size(match_ind,2);
%%%%%
%             output=OuT;           %%%% for surfacing
end

% set(gcf,'Color',[1 1 1])
% format_figure(im);

% draw scale of features
% sc = 2;
% for i = 1:length(ind)
%     h = plot(keyVector(match_ind(ind(i)),2), keyVector(match_ind(ind(i)),1),...
%         'ro','MarkerSize',sc*keyVector(match_ind(ind(i)),3));
%     h = plot(keyVector_m(match_ind_m(ind(i)),2), keyVector_m(match_ind_m(ind(i)),1),...
%         'ro','MarkerSize',sc*keyVector_m(match_ind_m(ind(i)),3));
% end;

% disp([gcf, sym_strength]');
% disp(sprintf('Figure %g: ', gcf));
% disp(sprintf('Symmetry Magnitude = %g  ', sym_strength));
% hold off
% drawnow

