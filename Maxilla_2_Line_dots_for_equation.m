% clc
% clear all
% close all
% imgpath='E:\Master_Project\DATA\impaction\yasaman moadab(mesiodens)\';
% [V,info]=ReadData3D([imgpath,'000.dcm']);
% info = dicom_read_header('E:\Softwares\DATA\impaction\yasaman moadab(mesiodens)\000.dcm');
% V = dicom_read_volume(info);
% IMAGE=load_nii('E:\Softwares\DATA\Projects_Test\vahid100a000.nii.nii');
% V=IMAGE.img;
m=min(V(:));
M=max(V(:));
I=(V-m)/M;
I=im2uint8(I);
[x,y,z]=size(I);
J=zeros(x,y,z);
angl_I=zeros(z,1);
angl_J=zeros(z,1);
Sym_Strength=zeros(z,1);
Sym_Center_line=zeros(z,1);
close all;
I_num=100;
unt=10;
X=zeros(512,512);
Y=zeros(512,512);
for i=I_num:I_num+unt
    clear surfingout
    clear x
    clear y
    clear p
     [surfingout,~,~,~,~]=symmetry(I(:,:,i),'mirror',1);
        x=surfingout(1,:);
        y=surfingout(2,:);
        p = polyfit(y,x,1);
        for j=1:512
            X(i,j) = polyval(p,j);
        end
end
        
yy=1:512;       yy=repmat(yy,[1,483]);
zz=zeros(size(yy));
for i=1:483
    zz(512*(i-1)+1:512*(i-1)+512)=i+15;
    xx(512*(i-1)+1:512*(i-1)+512)=X(i+15,:);
end

sf = fit([zz',yy'],xx','poly11');
plot(sf,[zz',yy'],xx')

PI=zeros(512,512,512);
for i=1:512
    for j=1:512
        ap=round(feval(sf,[i,j]));
        PI(i,ap,j)=255;
    end
end


IMAGE=load_nii('E:\Master_Project\DATA\Projects_Test\yasaman.nii');

for i=1:512
    for j=1:512
        ap=round(feval(sf,[i,j]));
        for k=1:512
            if k==ap
                IMAGE.img(i,k,j)=3070;
            else
%                 IMAGE.img(i,k,j)=-940;
            end
        end
    end
end

save_nii(IMAGE,'yasaman moadab.nii');
   