
function [surfingout,segmout]=meisam_display_output_6(im_dark,max_r,max_ang,ind,sym_x,sym_y,left_ind,right_ind,pts,pts_m,sym_strength);
% display matched pairs of mirrored keypoints associated with dominant axis

% figure(6);
% imshow(im_dark);
% hold on;
% color =  colormap(gray); % colormap(hsv);
centre = size(im_dark)/2;
tol = 1500; %10;
% if mean(size(im_dark)) < 400
%     linewidth = mean(size(im_dark))*0.00005; %*0.01; %20;%10;
%     markersize = mean(size(im_dark))*0.000035; %*0.075; %20;%10;
% else
%     linewidth = mean(size(im_dark))*0.000075; %*0.015; %20;%10;
%     markersize = mean(size(im_dark))*0.0005; %*0.1; %20;%10;
% end

for i = 1:1%length(max_r)

    if isempty(ind{i}),surfingout=[0;0];  segmout=mean(0);  continue, end

%     col = color(floor((i-1)/length(ind)*size(color,1))+1, :);
%     shade = sym_strength(i)/max(sym_strength); %max([sym_strength(i)/max(sym_strength), min([sym_strength(i)/0.3, 1])]);
%     col = [1 1 1] * shade;
%     c = mod(color_cycle_rate*i,size(color,1))+1;

    % construct set of points describing current symmetry axis and draw
    % axis.
    if abs(diff(unwrap([max_ang(i),pi]))) < (pi/4)
        Y = [min(sym_y(ind{i})) - tol : ...
            (max(sym_y(ind{i})) - min(sym_y(ind{i})))/100 : ...
            max(sym_y(ind{i})) + tol];
        X = ((-Y + centre(1))*sin(max_ang(i)) + max_r(i)) / cos(max_ang(i)) + centre(2);
%         display('F                  I                            RSTS')
    else
        X = [min(sym_x(ind{i})) - tol : (max(sym_x(ind{i})) - min(sym_x(ind{i})))/100: max(sym_x(ind{i})) + tol];
        Y = ((-X + centre(2))*cos(max_ang(i)) + max_r(i)) / sin(max_ang(i)) + centre(1);
%         display('S   E    C    O   ND')
    end
    
    p_ind = find((X>=(min(sym_x(ind{i}))-tol)) & (X<=(max(sym_x(ind{i}))+tol)) & ...
        (Y>=(min(sym_y(ind{i}))-tol)) & (Y<=(max(sym_y(ind{i}))+tol)));
    p_ind = find((X>=(min(sym_x(ind{i}))-tol)) & (X<=(max(sym_x(ind{i}))+tol)) & ...
        (Y>=(min(sym_y(ind{i}))-tol)) & (Y<=(max(sym_y(ind{i}))+tol)));
%      plot(X(p_ind),Y(p_ind),'y-','LineWidth', linewidth,'Color',col);
%     save('X.mat')
    meiX=X(p_ind);
    meiY=Y(p_ind);
    surfingout=[meiX;meiY];     %%%%%%%%  for surfacing
     segmout=mean(X);        %%%%% for segmenting
     if (isempty(X))
           segmout= max(sym_x(ind{i}));
           meiX=max(sym_x(ind{i}));
           meiY=Y(p_ind);
           surfingout=[meiX;meiY];
     end
%      hold off
    
    %%% meisam _ 2:
%     ss=zeros(512,512);
%     sss=ss;
%     figure(333);imshow(ss)
% 
%     hold on
%     plot(X(p_ind),Y(p_ind),'y-','LineWidth', linewidth,'Color',col); 
%     hold off
%     
%     display(X(p_ind))
%     display(Y(p_ind))
%     display(round((Y(p_ind)-min(min(Y(p_ind))))))
%     sss(round(X(p_ind)),round((1+Y(p_ind)-min(min(Y(p_ind))))))=1;
%     figure(13444);imshow(sss)
%     save('line.mat')
    %%% meisam _ 2
    
    
% % %     draw left and right point sets (currently uses the same colour for left and right sets)
%     plot(pts{i}(right_ind{i},1), pts{i}(right_ind{i},2), 'r.','Color',col,'MarkerSize',markersize);
%     plot(pts_m{i}(left_ind{i},1), pts_m{i}(left_ind{i},2), 'r.','Color',col,'MarkerSize',markersize);
%     plot(pts{i}(left_ind{i},1), pts{i}(left_ind{i},2), 'g.','Color', col,'MarkerSize',markersize);
%     plot(pts_m{i}(right_ind{i},1), pts_m{i}(right_ind{i},2), 'g.','Color', col,'MarkerSize',markersize);
end
end