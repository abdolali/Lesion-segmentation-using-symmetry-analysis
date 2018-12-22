clc
% clear all
% % % % load('C:\Users\meisam\Documents\Master_Project\DATA\lesion\romina keshavarz\romania.mat');
% load('C:\Users\meisam\Documents\Master_Project\DATA\impaction\amir hosein piriyai\amir.mat');
 load('C:\Users\meisam\Documents\Master_Project\DATA\impaction\sajad akrami\sajad.mat');
% % load('C:\Users\meisam\Documents\Master_Project\DATA\gholami\gholami.mat');
% J3=J;
% J=J2;
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\DATA\impaction\amir hosein piriyai\piriyai.nii');
% I_disorder = zeros(size(IMAGE.img));
% 
% Main settings
% main.similarity='SSD';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI
main.similarity='SSD';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI
main.MIbins=64;         % number of bins for the Mutual Information similarity measure
main.subdivide=3;       % use 3 hierarchical levels
main.lambda = 0.2;    % transformation regularization weight, 0 for none
optim.maxsteps = 70;   % maximum number of iterations at each hierarchical level
optim.fundif = 4e-10;    % tolerance (stopping criterion)
optim.gamma = 0.4;       % initial optimization step size 
optim.anneal = 0.05;       % annealing rate on the optimization step  
main.okno = 10 ;            % mesh window size, the smaller it is the more complex deformations are possible
hg = fspecial('gaussian', 5 ,2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% bwAC=zeros(x,y,z);
middle = y/2;
whole = y;
I_num  = 243;   unt = 0 ;
for i = I_num : I_num + unt
    close all ;

    if ~isnan(xx(i))

    Im=zeros(x,y);
    Im_f=Im;
    Im_fr=Im;
    clear left right lenght L Zp q ss qq c CC CH im im1 im2 num_Pixels idx_s ss0 BW ssi
    clear img img1 s Ch extra_cols imdi Im_f ACmask 
    
    a = imfilter(J(:,:,i),hg);
    a =imadjust(a);

    L=round(xx(i));
    if L>middle
       lenght=whole-L;
        left=a(:,L-lenght+1:L);
        right=a(:,L+1:end);
        extra_cols=L-middle;
%     extra_cols=(L-middle)*2;
    else if L==middle
            lenght=L;
            left=a(:,1:L);
            right=a(:,L+1:end);
            extra_cols=0;
        else
            lenght=L;
            left=a(:,1:L);
            right=a(:,L+1:L+lenght);
            extra_cols = L-middle;
%             extra_cols=(L-middle)*2;
        end
    end
    Zp=(whole-2*lenght)/2;
    Im(:,Zp+1:Zp+lenght)=left;
    Im(:,Zp+1+lenght:whole-Zp)=right(:,1:(whole-2*Zp-lenght));
    Im_f(:,1:middle)=Im(:,1:middle);
    for j=middle+1:whole
        Im_f(:,whole-j+middle+1)=Im(:,j);
        Im_fr(:,whole-j+1)=Im(:,j-middle);
    end
    
Im_f=double(Im_f)/255;
Im_fr=double(Im_fr)/255;

Im_r=zeros(size(Im_f));
Im_r(:,1:middle)=Im_f(:,1:middle);

Im_rr=zeros(size(Im_fr));
Im_rr(:,1:middle)=Im_fr(:,1:middle);

im1=Im_f(:,1:middle);
im2=Im_f(:,middle+1:whole);

[res, im]=mirt2D_register(im1, im2,main,optim);
   
figure;imshow(im),title('reged');
figure;imshow(im1),title('ref im');
figure;imshow(im2),title('float im');

Im_r(:,middle+1:end)=im;
imdi=(im-im2);

%%    Difference

img1=(im1);%>0.12;
img=(im);%>0.12;

c=(im2-im1) ;
d=(im-im1);
%  d=(im1-im) ;
% d=c;

 %%%%   for 512:
% d(380:end,:)=0;
% d(1:70,:)=0;
% d(:,1:90)=0;
 
%%%   for 462:
% d(250:end,:)=0;
% d(1:50,:)=0;
d(:,1:70)=0;

%  figure;imshow(img1),title('img1')
%  figure;imshow(im),title('img');

 figure;imshow(c),title('diff')
 figure;imshow(d),title('regdiff');

s=c>0.09;
ss0=d>0.2;

figure; imshow(ss0),title('thresholded');

ss = imerode(ss0,strel('disk',1));
ss = imdilate(ss,strel('disk',5));
ss = imerode(ss,strel('disk',5));
ss = imdilate(ss,strel('disk',3));

ss=imclose(ss, strel('disk',4));

ss = imreconstruct(ss, ss0);

figure;imshow(ss),title('regdifferodeddilat');
q=ss;
qr=imfilter(q,hg);
q=imfilter(q,hg);
CC = bwconncomp(q,18);
qq=q;
numPixels = cellfun(@numel,CC.PixelIdxList);
[num_Pixels,idx_s]=sort(numPixels,'descend');
clear Ch BW qr
CH=zeros(size(q));
jj_jj=1;
for jj=1:size(num_Pixels,2)
    if num_Pixels(jj)>0 && jj_jj<=20
        qr=q;
        qr(CC.PixelIdxList{idx_s(jj)}) = 0;
        Ch(:,:,jj)=qq-qr;
        Ch(:,:,jj)=bwconvhull(Ch(:,:,jj),'objects',8);
        CH=CH+Ch(:,:,jj);
        if jj_jj==1
            break
        end
        jj_jj=jj_jj+1;
    end
end

BW = edge(CH);
for ii=1:middle
    F1(:,ii+middle,i)=BW(:,middle+1-ii);
end

BW=(double(BW));
ACmask=zeros(x,y);
nnn=size(CH,2);
nnsh = 1;
CH2 = zeros(size(CH));
CH2(:,1:nnn-nnsh+1) = CH(:,nnsh:end);
for ii=1:middle
    F2(:,ii+middle,i)=BW(:,middle+1-ii);
    F2(:,ii,i)=BW(:,ii);
%     ACmask(:,ii+middle)=CH2(:,middle+1-ii);
%         ACmask(:,ii)=CH2(:, ii);
        ACmask(:,ii+middle)=CH(:,middle+1-ii);
%         ACmask(:,ii)=CH(:,ii);
end

figure, imshow(F2(:,:,i)),title('F2');

FF=zeros(x,whole);
F(:,:,i)=FF;
% ACmask(160:end, :) = 0;
% ACmask(:,320:end) = 0;
ACMASK = zeros(size(F));
if extra_cols>0
%     F(:,extra_cols+middle+1:whole,i)=255*F2(:,middle+1:middle+lenght,i);
%     F(:,extra_cols:extra_cols+middle-1,i)=255*F2(:,1:middle,i);
        F(:,extra_cols+middle+1:whole,i)=255*F2(:,middle+1:middle+lenght,i);
    F(:,extra_cols:extra_cols+middle-1,i)=255*F2(:,1:middle,i);
%     ACMASK(:,extra_cols+middle+1:whole,i)=255*ACmask(:,middle+1:middle+lenght,i);
 else if extra_cols<0
         extra_cols
         extra_cols=-extra_cols;
         F(:,L+1:2*L,i)=255*F2(:,middle+1:L+middle,i);
         F(:,1:L,i)=255*F2(:,1:L,i);
%          ACMASK(:,L+1:2*L,i)=255*ACmask(:,middle+1:L+middle,i);
    end
end

F2(:,1:middle,i)=F2(:,1:middle,i)+im1;

FFF(:,:,i) = imrotate(F(:,:,i),-angl_J(i),'bicubic','crop');
FFF(:,:,i) = imrotate(FFF(:,:,i),-angl_I(i),'bicubic','crop');

% close all
hgf = fspecial('gaussian', 15 ,1);
hgf2 = fspecial('log', 70 ,0.51);
af = I(:,:,i);
% af = 255 -I(:,:,i);
af= imfilter(af,hgf);
 af= imfilter(af,hgf2);
 af=af+imadjust(af);

ACtemp = zeros(size(ACmask));
% ACtemp(60:200,230:360) = ACmask(60:200,230:360);
ACtemp = ACmask;
% ACtemp = FFF(:,:,i);
ACmask = ACtemp;

% ACm=activecontour(af,ACmask,40);
% bwAC(:,:,i) = ACm;
bwAC(:,:,i)=ACmask;
%      figure, imshow((af)),title('filtered version of the original image');
%      figure, imshow(ACMASK),title('whole masks');
%       figure, imshow(uint8(I(:,:,i))),title('original');
%   figure, imshow(uint8(FFF(:,:,i))),title('End of End');
% %   figure, imshowpair(uint8(FFF(:,:,i)),uint8(I(:,:,i))),title('pair of original and the contour');
%   figure, imshowpair(uint8(bwAC(:,:,i)),uint8(I(:,:,i))),title('pair of original and the final contour');
% imgs1=uint8(FFF(:,:,i))+uint8(I(:,:,i));

% Cf = imfuse(255*edge(bwAC(:,:,i)),I(:,:,i),'ColorChannels',[2 1 0]);
%%%%%%%%%%%%%%%%%%%%55
% for iii = 1:253
%       bbb(:,iii) = aaa(:, 254-iii);
% end      

contou = zeros(x,y,3);
contou(:,:,2) = 255*edge(bwAC(:,:,i)) ;
Cf = imfuse(contou,I(:,:,i),'ColorChannels',[2 1 0]);
figure, imshowpair(contou, I(:,:,i), 'blend'),title('pair of original and the final contour');


disord = 255*(bwAC(:,:,i));
disord_2 = uint8(zeros(size(disord')));

    for iix=1 : x
        for jjx = 1 : y
            I_disorder(jjx,(x-iix)+1,i) = disord(iix,jjx);
        end
    end
    
% figure;imshow(I_disorder(:,:,i)),title('for nii');

figure;imshow(Cf),title('pair');
imgs1=(Cf);
    imwrite(imgs1,['C:\Users\meisam\Documents\Master_Project\DATA\impaction\sajad_R\',num2str(i),'.bmp']);
% imwrite(imgs1,['C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\hamideh_R\',num2str(i),'.bmp']);
% imwrite(imgs1,['C:\Users\meisam\Documents\Master_Project\DATA\impaction\piriyai_R\R2\',num2str(i),'.bmp']);
end
end
% I_disorder(:,:,i)=zeros(size(I_disorder(:,:,i)));
IMAGE.img=I_disorder;
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\hamideh_R\hami_auto.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\DATA\impaction\piriyai_R\R2\amir_auto.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\deylamidata\Shafayi jokandan_Maryam Y19(Dr Adham)_1993_12_26\R\joka_auto.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\DATA\impaction\sajad_R\sajad_dis.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\DATA\lesion\roman2.nii');

% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\DATA\gholami\R\gholami_auto.nii');