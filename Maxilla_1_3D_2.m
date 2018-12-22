% % % % % % % % % % % %           Maxilla     1
clc
clear all
close all
% % imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\Janzadeh_Mahtab  28Y (Dr.Moshtaghi)_1983_8_23\30if17a23000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Sadeghzadeh_Elyas 21Y  (Dr Sigaroudi)_1991_8_9\2j7a08be3000.dcm';
% % imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Derogar_Mohamad 55y (Dr Motevaseli)_1957_5_29\2j7a07ig3000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Eprib_Zahra 40Y ( Dr.Adham)_1971_1_12\30if18073000.dcm';
% imgpath='C:\Users\meisam\Documents\Master_Project\deylamidata\Khademi_Sareieh 47 Y(Dr Nezhati )_1963_6_29\30if158d3000.dcm';
% % imgpath='E:\Master_Project\DATA\impaction\sajad akrami\';    imgname ='000.dcm';
% 
% % imgpath='C:\Users\meisam\Documents\Master_Project\DATA\cleft\masoome babaeii\000.dcm';
% % imgpath='C:\Users\meisam\Documents\Master_Project\DATA\cleft\mohammad mahdi arab\000.dcm';
% % imgpath='C:\Users\meisam\Documents\Master_Project\DATA\cleft\vahid farokh\000.dcm';
imgpath='C:\Users\meisam\Documents\Master_Project\DATA\impaction\amir hosein piriyai\000.dcm';
% 
% % [V,info]=ReadData3D([imgpath,imgname]);
% 
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Saeidi_S.Hamideh 31Y (Dr.Adham)_1981_9_20\1 (1).dcm';
[V,info]=ReadData3D(imgpath);
% 
% % info = dicom_read_header('E:\Softwares\DATA\cleft\masoome babaeii\000.dcm');
% % V = dicom_read_volume(info);
% % IMAGE=load_nii('E:\Softwares\DATA\cleft\mohammad mahdi arab\arab.nii\arab.nii');
% % V=IMAGE.img;
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
close all;
I_num=1;
unt=z-1;
for i=I_num:I_num+unt
    i
        [~,~,~,~ , symangle]=symmetry(I(:,:,i),'mirror',1);
        
        angl_I(i)=symangle;
        angl_I(i)=180*angl_I(i)/pi;
        if angl_I(i)>90
            angl_I(i)=angl_I(i)-180;
        end
        J(:,:,i) = imrotate(I(:,:,i),angl_I(i),'bicubic','crop');
end
J=uint8(J);
% J=uint16(J);

for i=I_num:I_num+unt
            [~,~,~,~ , symangle]=symmetry(J(:,:,i),'mirror',1);
            
            angl_J(i)=symangle;
            angl_J(i)=180*angl_J(i)/pi;
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

% save('C:\Users\meisam\Documents\R5\amir_hosein_piriyai.mat');

%%%%%%%%%%%%%%%%%%%%%%%

% load('C:\Users\meisam\Documents\R5\amir_hosein_piriyai.mat');
  load('C:\Users\meisam\Documents\Master_Project\latest_data\reza_ghasemi.mat');
 
Sym_Strength=zeros(z,1);
Sym_Center_line=zeros(z,1);
close all;
X=zeros(x,y);
Y=zeros(x,z);
line_coe = zeros(z,2);

for i=1:z
    i
        clear surfingout
        clear xs
        clear ys
        clear ps
        [surfingout,~,~,~,~]=symmetry(J(:,:,i),'mirror',1);
        ys=surfingout(1,:);
        xs=surfingout(2,:);
        ps = polyfit(xs,ys,1);
        line_coe(i,:)=ps;
        for j=1:x
%             if polyval(ps,j) <200  || polyval(ps,j) > 270
%                 Y(j,i) = 0;
%             else
%                 Y(j,i) = polyval(ps,j);
%             end
            Y(j,i) = polyval(ps,j);
        end
end
for j=1:x
    mm=median(Y(j,:));
    for i=1:z
        if Y(j,i)<mm-5 || Y(j,i) > mm+5
            Y(j,i) = mm;
        end
    end
end
            
xx=1:x;       
xx=repmat(xx,[1,z]);
% zz=zeros(size(xx));
clear zz yy
for i=1:z
    i
    zz(x*(i-1)+1:x*(i-1)+x)=i;
    yy(x*(i-1)+1:x*(i-1)+x)=Y(:,i)';
end

sf = fit([xx',zz'],yy','poly11');
plot(sf,[xx',zz'],yy')

% 
% ii = 1;
% for i = 1 : size(yy,2)
%     if yy(i)~=0
%         xxx(ii) = xx(i);
%         yyy(ii) = yy(i) ;
%         zzz(ii) = zz(i) ;
%         ii = ii +1;
%     end
% end
% 
% sf = fit([yyy',zzz'],xxx','poly11');
% plot(sf,[yy',zz'],xx')

%%%%    in DICOM "x" is vertical and "y" is horizontal but in NIFTII they're reverse.

PI=zeros(x,y,z);
for i=1:x
    for j=1:z
        ap=round(feval(sf,[i,j]));
        PI(ap,i,j)=255;
    end
end

IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\DATA\impaction\amir hosein piriyai\piriyai.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Saeidi_S.Hamideh 31Y (Dr.Adham)_1981_9_20\1 (1).nii');

% for i=1:x
%     for j=1:z
%         ap=round(feval(sf,[i,j]));
%         for k=1:y
%             if k==ap
%                 IMAGE.img(k,i,j)=3070;
%             end
%         end
%     end
% end
% % 
% % %  IMAGE.img = PI;
% 
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\DATA\impaction\amir hosein piriyai\piriyai_sf0.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a3 = zeros(z,1);
for j = 1: z
    a1 = feval(sf,[1,j]) ;
    a2 = feval(sf,[x,j]) ;
    a3(j) = (a1+a2)/2;
end

hg = fspecial('gaussian', 5 ,2);
middle = round(y/2);
whole = y;

I3D=zeros(y,x,z);
J3D=zeros(y,x,z);
J2 = zeros(y,x,z);

for zcont = 1:z
    a = imfilter(J(:,:,zcont),hg);
    b = uint8(zeros(size(a')));
    for i=1 : x
        for j = 1 : y
            b(j,(x-i)+1) = a(i,j);
        end
    end
    J2(:,:,zcont) = b ;
end

for i=1:z
    i
       L_diff = middle - round(a3(i));
%         a = imfilter(J(:,:,i),hg);
        if L_diff > 0
            I3D(L_diff+1:whole,:,i) =J2(1:whole-L_diff,:,i);
        else if  L_diff<0
                I3D(1:whole+L_diff,:,i) =J2(1-L_diff:whole,:,i);
            else
                I3D(:,:,i) =J2(:,:,i);
            end
        end
end

IMAGE.img = I3D;
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\DATA\impaction\amir hosein piriyai\IMAGE.nii');

for i=1:z
    J3D(1:middle,:,i) = I3D(1:middle,:,i);
    for j = 1:middle
        J3D(whole-j+1,:,i) = I3D(middle+j,: ,i);
    end
end

J3D = uint16(J3D);
% test=load_nii('E:\2j.nii');

% left_nii = IMAGE;
% left_nii.hdr.dime.dim(2) = middle ;
% left_nii.hdr.hist.originator(1) = middle ;
% left_nii.original.hdr.dime.dim(2) = middle;
% left_nii.img = J3D(1:middle,:,:);

left_nii = IMAGE;
left_nii.hdr.dime.dim(2) = middle ;
left_nii.hdr.dime.dim(4) = 150;
left_nii.hdr.hist.originator(1) = middle ;
left_nii.original.hdr.dime.dim(2) = middle;
left_nii.original.hdr.dime.dim(4) = 150;
left_nii.img = J3D(1:middle,:,1:150);

% % % save_nii(left_nii,'C:\Users\meisam\Documents\Master_Project\DATA\L_romania.nii');
% save_nii(left_nii,'C:\Users\meisam\Documents\Master_Project\DATA\L_hamideh.nii');
% save_nii(left_nii,'C:\Users\meisam\Documents\Master_Project\DATA\L_amir_full.nii');
% save_nii(left_nii,'C:\Users\meisam\Documents\Master_Project\DATA\L_amir.nii');
save_nii(left_nii,'C:\Users\meisam\Documents\Master_Project\latest_data\L_rezghasem.nii');
% right_nii = left_nii;
% right_nii.img = J3D(middle+1:whole,:,:);
right_nii = left_nii;
right_nii.img = J3D(middle+1:whole,:, 1:150);
% % % save_nii(right_nii,'C:\Users\meisam\Documents\Master_Project\DATA\R_romania.nii');

% save_nii(right_nii,'C:\Users\meisam\Documents\Master_Project\DATA\R_hamideh.nii');
save_nii(right_nii,'C:\Users\meisam\Documents\Master_Project\latest_data\R_rezghasem.nii');
% save_nii(right_nii,'C:\Users\meisam\Documents\Master_Project\DATA\R_amir_full.nii');
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
 [res, im_3d]=mirt3D_register(right_nii.img, left_nii.img, main, optim);
% 
% i3d_result = left_nii;
% i3d_result.img = im_3d;
% 
% i3d_result.img = uint16(i3d_result.img);
% 
% save_nii(i3d_result,'C:\Users\meisam\Documents\Master_Project\DATA\registered_romania.nii');
%     
% i3d_result.img = right_nii.img - left_nii.img;
% 
%             