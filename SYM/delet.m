matitk('?')
matitk('f')
matitk('s')
matitk('r')
%------
load mri;
D=squeeze(D);
b=matitk('FCA',[5 0.0625 3],double(D));
c=matitk('SCC',[1.4 10 255],double(b),[],[102 82 25]);subplot(131);
imagesc(squeeze(D(:,:,15)));
axis image;
colormap gray
subplot(132);
imagesc(squeeze(b(:,:,15)));
axis image;
colormap gray
subplot(133);

imagesc(squeeze(c(:,:,15)));
axis image;
colormap gray
%------
