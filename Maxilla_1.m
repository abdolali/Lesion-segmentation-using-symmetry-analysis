% % % % % % % % % % % % %           Maxilla     1
clc
clear all
close all
% imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\Janzadeh_Mahtab  28Y (Dr.Moshtaghi)_1983_8_23\30if17a23000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Sadeghzadeh_Elyas 21Y  (Dr Sigaroudi)_1991_8_9\2j7a08be3000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Derogar_Mohamad 55y (Dr Motevaseli)_1957_5_29\2j7a07ig3000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Eprib_Zahra 40Y ( Dr.Adham)_1971_1_12\30if18073000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\Khademi_Sareieh 47 Y(Dr Nezhati )_1963_6_29\30if158d3000.dcm';
% imgpath='E:\Master_Project\DATA\impaction\sajad akrami\';    imgname ='000.dcm';

% imgpath='C:\Users\meisam\Documents\Master_Project\DATA\cleft\masoome babaeii\000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\DATA\cleft\mohammad mahdi arab\000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\DATA\cleft\vahid farokh\000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\DATA\impaction\amir hosein piriyai\000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\DATA\lesion\romina keshavarz\000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\DATA\lesion\romina keshavarz\000.dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\DATA\impaction\sajad akrami\000.dcm';
% [V,info]=ReadData3D([imgpath,imgname]);

% imgpath='C:\Users\meisam\Documents\Master_Project\DATA\gholami\1 (1).dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\Shafayi jokandan_Maryam Y19(Dr Adham)_1993_12_26\1 (1).dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Saeidi_S.Hamideh 31Y (Dr.Adham)_1981_9_20\1 (1).dcm';
% imgpath ='C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Gheravi_Mohammad  39Y (Dr.Sigaroodi)_1972_12_1  -   odd number of pixel\1 (1).dcm';
%  imgpath ='C:\Users\meisam\Documents\Master_Project\latest_data\reza ghasemi\000.dcm';
% imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\behrooz nazari\000.dcm';
% imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\kamal khanizadeh\000.dcm';
%  imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\ghazal soltani\000.dcm';
imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\zahra zadehnasir\000.dcm';

[V,info]=ReadData3D([imgpath]);

% info = dicom_read_header('E:\Softwares\DATA\cleft\masoome babaeii\000.dcm');
% V = dicom_read_volume(info);
% IMAGE=load_nii('E:\Softwares\DATA\cleft\mohammad mahdi arab\arab.nii\arab.nii');
% V=IMAGE.img;

m=min(V(:));
M=max(V(:));
I=(V-m)/M;
I=im2uint8(I);
% I=V;
[x,y,z]=size(I);
J=zeros(x,y,z);
% J=im2uint16(J);
angl_I=zeros(z,1);
angl_J=zeros(z,1);
sym_str1=zeros(z,3);
sym_str2=zeros(z,3);
numKeys=zeros(z,2);
close all;
I_num=1;
unt=z-1;
for i=I_num:I_num+unt
    i
        [~,~,sym_measure,~ , symangle, numkeys]=symmetry(I(:,:,i),'mirror',1);
        sym_str1(i,:) = sym_measure';
        numKeys(i,1) = numkeys;
        angl_I(i)=symangle;
        angl_I(i)=180*angl_I(i)/pi;
        if angl_I(i)>90
            angl_I(i)=angl_I(i)-180;
        end
        J(:,:,i) = imrotate(I(:,:,i),angl_I(i),'bicubic','crop');
end
J=uint8(J);
% J=uint16(J);
J2=J;
for i=I_num:I_num+unt
            [~, ~, sym_measure, ~ , symangle, numkeys]=symmetry(J(:,:,i),'mirror',1);
            sym_str2(i,:) = sym_measure';
            numKeys(i,2) = numkeys;
            angl_J(i)=symangle;
            angl_J(i) =180 * angl_J(i) / pi;
            if angl_J(i)>90
                angl_J(i)=angl_J(i)-180;
            end
            J(:,:,i) = imrotate(J(:,:,i),angl_J(i),'bicubic','crop');
end

xx=zeros(z,1);

for i=I_num:I_num+unt
   
    [~,segmout,~,~,~]=symmetry(J(:,:,i),'mirror',1);
    xx(i)=segmout;
    
end

 xxi = 1;
for i=1:z
    if (~isnan(xx(i)))
        xxx(xxi) = xx(i);
        xxi = xxi+1;
    end
end
for i =1:z
    if (abs(xx(i) - median(xxx))>5) || (isnan(xx(i)))
        xx(i)=median(xxx);
    end
end


%  save('C:\Users\meisam\Documents\Master_Project\DATA\lesion\romina keshavarz\romania.mat');
%  save('C:\Users\meisam\Documents\Master_Project\DATA\impaction\amir hosein piriyai\amir.mat');
%  save('C:\Users\meisam\Documents\Master_Project\DATA\impaction\sajad akrami\sajad.mat');
%  save('C:\Users\meisam\Documents\Master_Project\DATA\lesion\romania_full.mat');
%  save('C:\Users\meisam\Documents\Master_Project\DATA\gholami\gholami.mat');

%  save('C:\Users\meisam\Documents\Master_Project\deylamidata\Shafayi jokandan_Maryam Y19(Dr Adham)_1993_12_26\jokandan.mat');
%  save('C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\hamideh_R\hamideh.mat');

% save('C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Gheravi_R\Gheravi.mat');
% save('C:\Users\meisam\Documents\Master_Project\latest_data\reza ghasemi\reza_ghasemi.mat');
% save('C:\Users\meisam\Documents\Master_Project\cyst\behrooz nazari\behrooz_nazari.mat');
% save('C:\Users\meisam\Documents\Master_Project\cyst\kamal khanizadeh\kamal_khanizadeh.mat');
% save('C:\Users\meisam\Documents\Master_Project\cyst\ghazal soltani\ghazal_soltani.mat');
save('C:\Users\meisam\Documents\Master_Project\cyst\zahra zadehnasir\zahra_zadehnasir.mat');
%%%      71  - 422
   
%%%%%%%%%%%%%%%%%%            Use the m.file "Maxilla_1_loop_after_flip"
%%%%%%%%%%%%%%%%%%            instead of the below block

% %%  Flipping :      
%  % Main settings
% main.similarity='RC';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI
% % main.MIbins=64;         % number of bins for the Mutual Information similarity measure
% main.subdivide=3;       % use 3 hierarchical levels
% main.okno=10;            % mesh window size, the smaller it is the more complex deformations are possible
% main.lambda = 0.001;    % transformation regularization weight, 0 for none
% % main.single=1;          % show mesh transformation at every iteration
% 
% % Optimization settings
% optim.maxsteps = 125;   % maximum number of iterations at each hierarchical level
% optim.fundif = 4e-4;    % tolerance (stopping criterion)
% optim.gamma = 0.15;       % initial optimization step size 
% optim.anneal=0.6;       % annealing rate on the optimization step  
% 
% 
%  main.okno=20;            % mesh window size, the smaller it is the more complex deformations are possible
%  F1=zeros(512,512,512);
%  F2=zeros(512,512,512);
%  
% for i=250:250 %I_num:I_num+unt
% 
% %  Im=uint16(Im);
%  Im=zeros(x,y);
%  Im_f=Im;
%     clear left right lenght L Zp q ss qq c CC CH im im1 im2 num_Pixels idx_s
%     clear img img1 s Ch
%     a=J(:,:,i);
%     L=round(xx(i));
%     if L>256
%        lenght=512-L;
%         left=a(:,L-lenght+1:L);
%         right=a(:,L+1:end);
%         extra_cols=L-lenght;
%     else if L==256
%             lenght=L;
%             left=a(:,1:L);
%             right=a(:,L+1:end);
%             extra_cols=0;
%         else
%             lenght=L;
%             left=a(:,1:L);
%             right=a(:,L+1:L+lenght);
%             extra_cols=L+lenght-512;
%         end
%     end
%     Zp=(512-2*lenght)/2;
%     Im(:,Zp+1:Zp+lenght)=left;
%     Im(:,Zp+1+lenght:512-Zp)=right(:,1:(512-2*Zp-lenght));
%     Im_f(:,1:256)=Im(:,1:256);
%     for j=257:512
%         Im_f(:,512-j+257)=Im(:,j);
%     end
% %     imshow(left),title('left')
% %     figure;imshow(right),title('R')
% %     figure;imshow(uint8(Im(:,:,i))),title('Im')
% %     figure;imshow(uint8(Im_f(:,:,i))),title('Im_f')
% 
% 
% 
% 
%  Im_f=double(Im_f)/255;
% Im_r=zeros(size(Im_f));
% Im_r(:,1:256)=Im_f(:,1:256);
% display('1');
%  % for i=105:105
%     im1=Im_f(:,1:256);
%     im2=Im_f(:,257:512);
%     [res, im]=mirt2D_register(im1, im2,main,optim);
%     Im_r(:,257:end)=im;
% % end
% 
% % res is a structure of resulting transformation parameters
% % newim is a deformed image 
% %
% % you can also apply the resulting transformation directly as
% % newim=mirt3D_transform(im, res);
% % 
% % figure,imshow((im1)); title('Reference (fixed) image slice');
% % figure,imshow((im2));title('Source (float) image slice');
% % figure,imshow((im)); title('Registered (deformed) image slice');
% 
% %%    Difference
% img1=(im1);%>0.2;
% img=(im);%>0.2;
% c=(im1-im2);
% d=(img1-img);
% % figure;imshow(c),title('diff')
% % figure;imshow(d),title('regdiff');
% 
% % level = graythresh(c);  %%% otsu
% % s = im2bw(c,level);    %%% otsu
% % 
% % level = graythresh(d);  %%% otsu
% % ss = im2bw(d,level);    %%% otsu
% s=c>0.1;
% % imshow(s)
% ss0=d>0.1;
% se = strel('disk',2);        
% ss = imerode(ss0,se);
% 
%  figure;imshow(s),title('diff')
%  figure;imshow(ss0),title('regdiff')
% figure;imshow(ss),title('erodedregdiff')
% 
% q=ss;
% q(:,end-10:end)=0 ; q(370:end,:)=0;
% CC = bwconncomp(q,18);
% qq=q;
% numPixels = cellfun(@numel,CC.PixelIdxList);
% % [biggest,idx] = max(numPixels);
% % q(CC.PixelIdxList{idx}) = 0;
% [num_Pixels,idx_s]=sort(numPixels,'descend');
% clear Ch BW qr
% CH=zeros(size(q));
% for jj=1:size(num_Pixels,2)
%     if num_Pixels(jj)>60
%         qr=q;
%         qr(CC.PixelIdxList{idx_s(jj)}) = 0;
%         Ch(:,:,jj)=qq-qr;
%         Ch(:,:,jj)=bwconvhull(Ch(:,:,jj));
%         CH=CH+Ch(:,:,jj);
%     end
% end
% %  qq=qq-q;
% % figure, imshow(qq);
% % figure, imshow(CH);
% BW = edge(CH);
% for ii=1:256
%     F1(:,ii+256,i)=BW(:,257-ii);
% end
% % figure, imshow(BW);
% BW=(double(BW))+im2double(im2);
% % figure, imshow(BW),title('BW');
% 
% for ii=1:256
%     F2(:,ii+256,i)=BW(:,257-ii);
% end
% if extra_cols>0
%     F(:,1:extra_cols,i)=a(:,1:extra_cols);
%     F(:,extra_cols+1:extra_cols+256,i)=left;
%     F(:,extra_cols+257:512,i)=right;
% end
% 
% F2(:,1:256,i)=im1;
% figure, imshow(F2(:,:,i)),title('End');
% end
% 
% % IMAGE.img=F2;
% % save_nii(IMAGE,'ssss.nii');
% 
% %   Otsu Thresh
% 
% %% Eval
% % 
% % T1=imread('E:\Softwares\DATA\Projects_Test\R260.png');
% % T1=logical(T1);
% % T1=T1(:,:,1);
% % T1=T1(:);
% % F=F(:);
% % common=sum(F&T1);
% % union=(F|T1);
% % d_f=sum(F);
% % d_t=sum(T1);
% % Dice=2*common/(d_f+d_t)
% % Jaccard=sum(common)/sum(union)
% % lmf=1;
% % [hd ,D]= HausdorffDist(T1,F, lmf);
% % hd
% 
% 
% 
% 
% %%          Evaluation
% % BW=(double(BW))+im2double(b);
% % for i=1:n
% %     a2(:,n+1-i)=a(:,i);
% % end
% % fi(:,1:size(a,2))=im2double(a2);
% % fi(:,1+size(a,2):2*size(a,2))=BW;
% % figure, imshow(fi)
% 
% % auto=imread('cleft_auto.bmp');
% % auto=auto(60:290,61:450);
% % aa=auto;
% % aa_d=aa(:);
% % bb = imfill(bb,'holes');
% % bb_d=bb;
% % bb_d=bb_d(:);
% % common=sum(aa_d & bb_d);
% % union=sum(aa_d | bb_d);
% % cm=sum(aa_d);
% % co=sum(bb_d);
% % Jaccard=common/union
% % Dice=(2*common)/(cm+co)
% % % rfp false pasitive ratio
% % % rfn false negative ratio
% % Y=aa;
% % X=bb;
% % %If X is your binary image and Y is the ideal binary image
% % TPR = common/cm
% % FPR = sum(and(X(:),not(Y(:)))/sum(not(Y(:))))
% % % For a perfect detection, TPR = 1, and FPR = 0.
% % 
% % rfp=(abs(co-common))/cm
% % rfn=(abs(cm-common))/cm
% % 
% % %Hausdorff distance
% % lmf=1;
% % [hd D]= HausdorffDist(aa, bb, lmf);
% % hd

