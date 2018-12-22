clc
clear all
close all
% % % dir('C:\Users\meisam\Documents\Master_Project\DATA');
% % % comm = 'reg_f3d -ref  L_romania -flo R_romania';
% % % system(comm);
% % % 
% % % imgpath1='L_romania.nii';
% % % imgpath2='R_romania.nii';
% % % imgpath3='outputResult.nii';
% % % 
% % VV=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Shafayi jokandan_Maryam Y19(Dr Adham)_1993_12_26\R\joka_auto.nii');
% % V1=VV.img;
% % load('C:\Users\meisam\Documents\Master_Project\deylamidata\Shafayi jokandan_Maryam Y19(Dr Adham)_1993_12_26\jokandan.mat');
% % VV=load_nii('C:\Users\meisam\Documents\Master_Project\DATA\lesion\romina keshavarz\rom\rom2.nii');
% % VV=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Shafayi jokandan_Maryam Y19(Dr Adham)_1993_12_26\R\joka_manu_2.nii');
% VV=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\hamideh_R\hamideh_manu.nii');
VV=load_nii('C:\Users\meisam\Documents\Master_Project\DATA\amir_manu.nii');
V2=VV.img;
V3=zeros(size(V2));
[y,x,z] = size(V2);
for k=1:z
    for j=1:y
        V3(y-j+1,:,k)=V2(j,:,k);
    end
%     k
end

%%%%%%%%%%%%%%%%%%%  for making new segmentation, uncomment this part zzz ;
imgpath1='L_amir.nii';
imgpath2='R_amir.nii';
imgpath3='out_amir.nii';

% [L,info]=ReadData3D([imgpath1]);
% [R,info]=ReadData3D([imgpath2]);
% [R_Led,info]=ReadData3D([imgpath3]);

Ls=load_nii(imgpath1);
Rs=load_nii(imgpath2);
R_Leds=load_nii(imgpath3);
L=(double(Ls.img));
R=(double(Rs.img));
R_Led=255*mat2gray(double(R_Leds.img));
regdiff = (L-R_Led);
diff = (R-L);
regdiff_n = 255* mat2gray(regdiff);
bw_regdiff = regdiff_n > 145;

nii_temp =  load_nii(imgpath1);
nii_temp.img = regdiff_n;
% save_nii(nii_temp,'C:\Users\meisam\Documents\Master_Project\DATA\regdiff_hamideh.nii');
nii_temp.img = 255*bw_regdiff;
% save_nii(nii_temp,'C:\Users\meisam\Documents\Master_Project\DATA\regdiff_hamideh_th.nii');

d5 = 5; e11 = 11; d11 = 11;
% dl_1 = ones([d5,d5,d5]);
% imdl_1 = imerode(bw_regdiff,dl_1);
% imdl_1 = imdilate(imdl_1,dl_1);
% er_1 = ones([e11,e11,e11]);
% imer_1 = imerode(imdl_1,er_1);
% dl_2 = ones([d11,d11,d11]);
% imdl = imdilate(imer_1,dl_2);
[ly,lx,lz] = size(bw_regdiff);
d5 =7; e11 = 9;
dl_1 = ones([d5,d5,d5]);
imdl_0= bw_regdiff;
imdl_0(ly-100:ly,:,:) = 0;      imdl_0(:,:,lz-2:lz) = 0;
imdl_1 = imerode(imdl_0,dl_1);
imdl_1 = imreconstruct(imdl_1,bw_regdiff,26);
 d11 = 3; dl_2 = ones([d11,d11,d11]);
imdl_1 = imdilate(imdl_1,dl_2);
d11 = 7; dl_2 = ones([d11,d11,d11]);
imdl = imerode(imdl_1,dl_2);

V1=255*imdl;
nii_temp.img = 255*imdl;
save_nii(nii_temp,'C:\Users\meisam\Documents\Master_Project\DATA\imdl.nii');

% % V1= I_disorder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of zzz.

nii_temp = load_nii('C:\Users\meisam\Documents\Master_Project\DATA\amir_auto_3d.nii');
V1 = nii_temp.img;


sigma= 10;
[~,Faces_2,Vertices_2,Faces2,Vertices2]=GetImageSurface(V2,sigma);
FV2.vertices = Vertices2;       FV2.faces = Faces2; 
close all
[~,Faces_1,Vertices_1,Facesl,Verticesl]=GetImageSurface(V1,sigma);
FV1.vertices = Verticesl;       FV1.faces = Facesl; 

figure,patch(FV1,'facecolor',[0 1 0],'edgecolor','none');
hold on; patch(FV2,'facecolor',[0 0 1],'edgecolor','none');
lighting phong
camlight
hold off;

figure; patch(FV2,'facecolor',[0 1 0],'edgecolor','none'),title('manual');
lighting phong
camlight

figure,patch(FV1,'facecolor',[0 0 1],'edgecolor','none'),title('algorithm');
lighting phong
camlight

% VV1= V1(:,:,1:48);

common=nnz(V2 & V1);
union=nnz(V1 | V2);
cm=nnz(V1);
co=nnz(V2);
Jaccard=common/union
Dice=(2*common)./(cm+co)
TPR = common/cm
FPR = sum(and(V1(:),not(V2(:)))/sum(not(V2(:))))

% VV1= V1(:,:,1:48);
VV1=V1;
common=nnz(V3(1:ly,:,:) & VV1);
union=nnz(VV1 | V3(1:ly,:,:));
cm=nnz(VV1);
co=nnz(V3);
Jaccard=common/union
Dice=(2*common)./(cm+co)
