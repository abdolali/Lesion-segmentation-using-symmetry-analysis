imgpath1='L_amir.nii';
imgpath2='R_amir.nii';
imgpath3='outputResult.nii';

VV=load_nii('C:\Users\meisam\Documents\Master_Project\DATA\amir_auto_3d.nii');
V2=VV.img(:,:,221:305);
V4=zeros(size(V2));
[y,x,z] = size(V2);
for k=1:z
    for j=1:y
        V4(y-j+1,:,k)=V2(j,:,k);
    end
%     k
end
V3=V4(ly+1:end,:,:);
Ls=load_nii(imgpath1);
Rs=load_nii(imgpath2);
R_Leds=load_nii(imgpath3);
L=(double(Ls.img));
R=(double(Rs.img));
R_Led=255*mat2gray(double(R_Leds.img));

RL=(R_Led);
regdiff = (RL-L);
diff = (R-L);
regdiff_n = 255* mat2gray(regdiff);
bw_regdiff = regdiff_n > 187;

nii_temp =  load_nii(imgpath1);
nii_temp.img = regdiff_n;
save_nii(nii_temp,'C:\Users\meisam\Documents\Master_Project\DATA\regdiff_amir.nii');
nii_temp.img = 255*bw_regdiff;
save_nii(nii_temp,'C:\Users\meisam\Documents\Master_Project\DATA\regdiff_amir_th.nii');

[ly,lx,lz] = size(bw_regdiff);
d5 =5; e11 = 11;
dl_1 = ones([d5,d5,d5]);
imdl_0= bw_regdiff;
imdl_0(1:20,:,:) = 0;   
imdl_0(:,1:140,:) = 0;  
% imdl_0(:,:,lz-2:lz) = 0;
imdl_1 = imerode(imdl_0,dl_1);
imdl_1 = imreconstruct(imdl_1,bw_regdiff,26);
 d11 = 11; dl_2 = ones([d11,d11,d11]);
imdl_1 = imdilate(imdl_1,dl_2);
% d11 = 7; dl_2 = ones([d11,d11,d11]);
imdl = imerode(imdl_1,dl_2);

V1=255*imdl;
    for j = 1:ly
        V0(ly-j+1,:,:) = V1(j,:,:);
    end
nii_temp.img = 255*imdl;
save_nii(nii_temp,'C:\Users\meisam\Documents\Master_Project\DATA\imdl_amir.nii');

sigma= 5;
[~,Faces_2,Vertices_2,Faces2,Vertices2]=GetImageSurface(V3,sigma);
FV2.vertices = Vertices2;       FV2.faces = Faces2; 
close all
[~,Faces_1,Vertices_1,Facesl,Verticesl]=GetImageSurface(V0,5);
FV1.vertices = Verticesl;       FV1.faces = Facesl; 

figure,patch(FV1,'facecolor',[0 0 1],'edgecolor','none');
hold on; patch(FV2,'facecolor',[0 1 0],'edgecolor','none');
lighting phong
camlight
hold off;

figure; patch(FV2,'facecolor',[0 1 0],'edgecolor','none'),title('manual');
lighting phong
camlight

figure,patch(FV1,'facecolor',[0 0 1],'edgecolor','none'),title('algorithm');
lighting phong
camlight

VV1= V0(:,:,1:z);
common=nnz(V3 & VV1);
union=nnz(VV1 | V3);
cm=nnz(VV1);
co=nnz(V3);
Jaccard=common/union
Dice=(2*common)./(cm+co)
% 
% %%%%%%%%%%%%%%%
% right_nii = load_nii('C:\Users\meisam\Documents\Master_Project\DATA\L_romania.nii');
% left_nii = load_nii('C:\Users\meisam\Documents\Master_Project\DATA\R_romania.nii');
% %%%%%%%%%%%%%%%%%%
% % Main settings
% main.similarity='SSD';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI
% % main.MIbins=64;         % number of bins for the Mutual Information similarity measure
% main.subdivide=3;       % use 3 hierarchical levels
% main.lambda = 0.001;    % transformation regularization weight, 0 for none
% % main.single=1;          % show mesh transformation at every iteration
% % Optimization settings
% optim.maxsteps = 70;   % maximum number of iterations at each hierarchical level
% optim.fundif = 4e-10;    % tolerance (stopping criterion)
% optim.gamma = 0.4;       % initial optimization step size 
% optim.anneal = 0.05;       % annealing rate on the optimization step  
% main.okno = 150;            % mesh window size, the smaller it is the more complex deformations are possible
% 



%  [res, im_3d]=mirt3D_register(right_nii.img, left_nii.img, main, optim);