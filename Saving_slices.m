clc
clear all 
close all
imgpath ='C:\Users\meisam\Documents\Master_Project\DATA\lesion\romina keshavarz\000.dcm';
[V,info]=ReadData3D([imgpath]);
m=min(V(:));
M=max(V(:));
I=(V-m)/M;
I=im2uint8(I);
[x,y,z]=size(I);

for i = 1:z
      imwrite(I(:,:,i),['C:\Users\meisam\Documents\Master_Project\PICS\Hami_Zahra\',num2str(i),'.bmp']);
end
