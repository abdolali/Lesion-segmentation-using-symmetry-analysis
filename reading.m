clear 
[V,info]=ReadData3D('E:\Softwares\DATA\lesion\romina keshavarz\000.dcm');
VVVV=(min(min(min(V))));
VVV=V-VVVV;
VV=VVV/max(max(max(V)));
imshow(VV(:,:,255))
imwrite(VV(:,:,255),'V.jpg');