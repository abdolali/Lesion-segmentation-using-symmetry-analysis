function [surfingout,segmout]= meisam_display_output(im,keyVector,ang,sym_x,sym_y, ind,match_ind,match_ind_m,Hough_im,keyVector_m, max_ang,max_r,left_ind,right_ind,pts, pts_m,phase_weight,sym_strength)
% Disply output for bilateral symmetry detection

color = hsv(256);
% color_cycle_rate = 50;

% DISPLAY OUTPUT
im_dark = uint8(round(double(im))); % im_dark = uint8(round(double(im)*0.55));

%         figure(1);
%         imshow(im); 
%         format_figure(im);
%         print('-djpeg','results/butterfly_eg1.jpg');
% 
        
        %%%%%%%%%%%    mh304 : try this to show :
% 
        figure; imshow(im_dark); 
        hold on;         draw_keypoints_fast(im, keyVector, [0 1 0]); 
        plot([keyVector(:,2), keyVector_m(:,2)]', [keyVector(:,1), keyVector_m(:,1)]', 'r.'),
        title('keypoints');
        hold off;

        figure;
        r = 5*ones(size(ang));
        [u,v] = pol2cart(ang+pi/2,r);
        imshow(im_dark);
        hold on,
        axis_width = 1;
        
        % draw sym axes
%         quiver(sym_x,sym_y,u,v,0,'g.'), title('symmetry axes');
%         Disable the above line and enable the lines below to make output
%         with colour of axis proportional to symmetry magnitude.  This
%         takes much longer to run!!
                for i = 1:length(sym_x)
                    col = [1 0 0];
                    h1 = quiver(sym_x(i),sym_y(i),u(i),v(i),0,'r.');
                    set(h1,'LineWidth',axis_width,'Color',phase_weight(i)*col);
                    h2 = quiver(sym_x(i),sym_y(i),-u(i),-v(i),0,'r.');
                    set(h2,'LineWidth',axis_width,'Color',phase_weight(i)*col);
                end

                title([num2str(length(sym_x)),' local symmetry axes']);
%         xlabel(sprintf('%g reflective matches',length(sym_x)));
        hold off;
        
        figure(4);
        imagesc(Hough_im), title('2D Hough space'),
%         hold on;
        
        % draw sym axis particles
        figure(5);
        set(gca,'FontSize',48);%, FontName, 'Times');
        r = 5*ones(size(ang));
        [u,v] = pol2cart(ang+pi/2,r);
        imshow(im_dark); 
        hold on,
        axis_width = 1.5;%3;
        for i = 1:length(ind)
            col1 = [1 0 0];  % col = color(floor((i-1)/length(ind)*size(color,1))+1, :);
            col2 = [0 1 0];
            col3 = [1 1 0];
            for j = 1:length(ind{i})
                h1 = quiver(sym_x(ind{i}(j)),sym_y(ind{i}(j)),u(ind{i}(j)),v(ind{i}(j)),0,'g.');
                set(h1,'LineWidth',axis_width,'Color',phase_weight(ind{i}(j))*col1);
                h2 = quiver(sym_x(ind{i}(j)),sym_y(ind{i}(j)),-u(ind{i}(j)),-v(ind{i}(j)),0,'g.');
                set(h2,'LineWidth',axis_width,'Color',phase_weight(ind{i}(j))*col1);

                %             % draw lines of correspondance
                            plot([keyVector(match_ind(ind{i}(j)),2), keyVector_m(match_ind_m(ind{i}(j)),2)]', ...
                                [keyVector(match_ind(ind{i}(j)),1), keyVector_m(match_ind_m(ind{i}(j)),1)]', ...
                                'g:','Color',phase_weight(ind{i}(j))*col3);
                
                            % draw feature points
                            plot([keyVector(match_ind(ind{i}(j)),2), keyVector_m(match_ind_m(ind{i}(j)),2)]', ...
                                [keyVector(match_ind(ind{i}(j)),1), keyVector_m(match_ind_m(ind{i}(j)),1)]', ...
                                'g+','Color',phase_weight(ind{i}(j))*col1);
            end

            % plot orientation of features
            draw_keypoints_fast(im,keyVector(match_ind(ind{i}),:));
            draw_keypoints_fast(im,keyVector_m(match_ind_m(ind{i}),:));

            %disp(sprintf('%g symmetric features',length(ind{i})));
        end
        title('sym axis particles(red lines) and lines of correspondances(yellow lines)');
%         hold off;

        % display 2D Hough space
        
      
        %%%%%%%%%%%%%%%%%   <----mh304 
        
        rad_upper_bound = round(sqrt((size(im,2)/2).^2 + (size(im,1)/2).^2));

        % convert sym_x and sym_y to polar coordinates with origin in centre of
        % image
        ang_h = mod(ang,pi);
        xx = sym_x - size(im,2)/2;
        yy = sym_y - size(im,1)/2;
        r = xx .* cos(ang_h) + yy .* sin(ang_h);

        % convert angles and radii to points in Hough vote image
        ang_hs = floor(180/pi*ang_h)+1;
        r_hs = round(r + rad_upper_bound);
        
        %%%%%%%%%%%%%%%%%%%%%    mh 304  . imshow :
%         plot(r_hs,ang_hs,'w.','MarkerSize',1);
%         for i = 1:length(ind)
%             col = [1 1 1];  %col = color(floor((i-1)/length(ind)*size(color,1))+1, :);
%             % c = mod(color_cycle_rate*i,size(color,1))+1;
%             plot(max_r(i) + rad_upper_bound, max_ang(i)*180/pi, 'go', 'MarkerSize',10, ...
%                 'Color',col);
%             plot(r_hs(ind{i}),ang_hs(ind{i}),'g.','Color',col,'MarkerSize',1);
%         end
%         colormap(gray);
%         hold off;
%         set(gcf,'Color',[1 1 1])
%         axis square; axis on;
%         set(gca,'FontSize',18);%, FontName, 'Times');
%         xlabel('r');
%         ylabel('\theta');
%         tick_labels = get(gca,'XTickLabel');
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   <-----mh304

r_extrema = 100;%size(Hough_im,2)/2,
 r_middle = size(Hough_im,2)/2;
 
 %%%%%%%%%%%%%%%%%%%%%%%  mh304. imshow :
 
%  
%         set(gca,'XTick',[r_middle-r_extrema, r_middle, r_middle+r_extrema],...
%             'XTickLabel',[-r_extrema,0,r_extrema]);
%         set(gcf,'InvertHardcopy','off');
%         print('-djpeg','results/butterfly_eg4.jpg');
% 
% 
%         figure(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%5     mh304  <-----

        % meisam St:
%        display(max_r);
%         display(max_ang);
%         display(ind);
%         save('sym_x.mat');
%         display(sym_y);
%         display(sym_strength);
%         display(left_ind);
%         display(right_ind);
%         display(r);
%         display(pts);
%         display(pts_m);
        % meisam End
        [surfingout,segmout]=meisam_display_output_6(im_dark,max_r,max_ang,ind,sym_x,sym_y,left_ind,right_ind,pts,pts_m,sym_strength);
        
         y=surfingout(1,:);
         x=surfingout(2,:);
         p = polyfit(x,y,1);
         x=1:512 ;
         y = polyval(p,x);
         plot(y,x,'Color',[0 0 1],'LineWidth',2);
         title('Global symmetry line');
         hold off;
%         format_figure(im);
%         print('-djpeg','results/butterfly_eg6.jpg');

%         clf;
%         [surfingout,segmout]=meisam_display_output_6(im_dark,max_r,max_ang,ind,sym_x,sym_y,left_ind,right_ind,pts,pts_m,sym_strength);
end
