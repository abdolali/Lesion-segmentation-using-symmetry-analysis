%%%%%%%%%%%%%%%%%%%%%   Dimension explaining :
%   in Dicom :
%   z == number of slices. It starts from neck and goes up to the cranial
%   y == represents number of columns in axial view. It starts from right side of the head to the left side
%   x == represents number of rows in axial view. It starts from the nose and ends at occipital
%
%   in niifti that is loaded with load_nii:
%   z == the same
%   y == is Dicom's "x"
%   x == is Dicom's "y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
beep off
% imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\afsaneh bakhshandeh(post)\000.dcm';
% imgpath ='C:\Users\meisam\Documents\Master_Project\latest_data\reza ghasemi\000.dcm';
% imgpath ='C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Derogar_Mohamad 55y (Dr Motevaseli)_1957_5_29\1 (1).dcm';
% imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\ali reza mahdavi rad\000.dcm';
% imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\behnam yousefi\000.dcm';
%  imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\zahra zadehnasir\000.dcm';
%  imgpath ='C:\Users\meisam\Documents\Master_Project\cyst\kamal khanizadeh\000.dcm';
% imgpath ='C:\Users\meisam\Documents\Master_Project\deylamidata\Jaereh_Ozra 58 Y(Dr Moshtaghi )_1954_1_1\1 (1).dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\Janzadeh_Mahtab  28Y (Dr.Moshtaghi)_1983_8_23\1 (1).dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\Hami_Zahra 35Y (dr.Ahmadian)_1976_5_9\1 (1).dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\Khademi_Sareieh 47 Y(Dr Nezhati )_1963_6_29\1 (1).dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\Kiaeifar_Arezoo 38Y (dr.Sigaroodi)_1973_9_1\1 (1).dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\Mirtalebi_Aghil 21Y  (Dr.Sigaroodi)_1991_6_17\1 (1).dcm';
% imgpath = 'C:\Users\meisam\Documents\Master_Project\deylamidata\Nemani_Shirin  9Y  (Dr.Sigaroodi)_2003_8_12\1 (1).dcm';
imgpath = 'C:\Users\meisam\Documents\Master_Project\DATA\lesion\romina keshavarz\000.dcm';
[V,info]=ReadData3D([imgpath]);

m=min(V(:));
M=max(V(:));
I=(V-m)/M;
I=im2uint8(I);
[x_di,y_di,z]=size(I);
J=zeros(x_di,y_di,z);
angl_I=zeros(z,1);
angl_J=zeros(z,1);
Sym_Strength=zeros(z,1);
Sym_Center_line=zeros(z,1);
close all;
I_num=1;
unt=z-1;
X=zeros(z,x_di);
sym_measure = zeros(3,z);
symstrength = zeros(1,z) ; symangle = zeros(1,z);

for i=I_num:I_num+unt
    clear surfingout
    clear x
    clear y
    clear p
     [surfingout, ~,symeasur, symstrength(i) , symangle(i)]=symmetry(I(:,:,i),'mirror',1);
     sym_measure(:,i) = symeasur;
        y=surfingout(1,:);
        x=surfingout(2,:);
        p = polyfit(x,y,1);
        for j=1:x_di
            X(i,j) = polyval(p,j);
        end
        i
end

j = 1; jj=1; jjj= 1;
for i = 1: z
    if (sym_measure(1,i))
        sym_measure_2(:,j) = sym_measure(:,i) ;
        j = j +1 ;
    end
    if (symstrength(i))
        symstrength_2(jj) = symstrength(i) ;
        jj = jj +1 ;
    end
%     if (symangle(i))
        if (symangle(i) < pi/2)
            symangle_2(jjj) = 180*symangle(i) / pi ;
        else if (symangle(i) > pi/2)
                symangle_2(jjj) = 180*(symangle(i) - pi) / pi  ;
            end
        end
        jjj = jjj +1 ;
%     end
end

sym_mean = mean(sym_measure_2,2)
sym_std(1) = std(sym_measure_2(1,:))    ;
sym_std(2) = std(sym_measure_2(2,:))     ;

symstrength_mean = mean(symstrength_2)      
symstrength_std = std(symstrength_2)    ;

symangle_mean = mean(symangle_2)    ;
symangle_std = std(symangle_2)            ;

all_median = mean(X(:));
all_std = std(X(:));
Y = X;
for i=1:z
    for j = 1:x_di
        if (X(i ,j) < all_median - 20 ||  X(i ,j) > all_median + 20)       
              Y(i ,j) = all_median;
        end
        if (X(i ,j) == 0 )
            Y(i ,j) = all_median;
        end
    end
end
xx=1:x_di;       xx=repmat(xx,[1,z]);
zz=zeros(size(xx));
for i=1:z
    zz(x_di*(i-1)+1:x_di*(i-1)+x_di)=i;
    yy(x_di*(i-1)+1:x_di*(i-1)+x_di)=Y(i ,:);
end

sf = fit([xx',zz'],yy','poly11');
plot(sf,[xx',zz'],yy')

PI=zeros(y_di,x_di,z);

qq = X;
for i=1:z
    for j = 1:x_di
            qq(i ,j) = round(X(i ,j));
        if (qq(i ,j) == 0 )
            qq(i ,j) = 1;
        end
    end
end

for i=1:z
    for j=1:x_di
        ap=round(feval(sf,[i,j]));
        PI(ap,j,i)=255;
%             PI(qq(i,j),j,i)=255;
    end
end

D = sf.p00 ; A = sf.p10 ; B = sf.p01; C = -1;
denum = sqrt(A*A + B*B + C*C) ; 
dd = zeros(z,x_di);
for i=1:z
    for j=1:x_di
          if (X(i,j))
                dd(i,j) = abs(D + i*A + j*B + C*X(i,j)) / denum;
          end
    end
end

SymPower =  100000000*mean(sym_measure_2(2,:)) / sum(dd(:)) ;
SymPower3D = sum(dd(:)) / (x_di * y_di * z )
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\second serei\Derogar_Mohamad 55y (Dr Motevaseli)_1957_5_29\Derogar_Mohamad.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\latest_data\reza ghasemi\rezaghasemi.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\cyst\afsaneh bakhshandeh(post)\afsaneh_bakhshandeh(post).nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\cyst\ali reza mahdavi rad\ali reza mahdavi rad.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\cyst\behnam yousefi\behnam_yousefi.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\cyst\zahra zadehnasir\zahra_zad.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\cyst\kamal khanizadeh\kamal_khanizadeh.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Jaereh_Ozra 58 Y(Dr Moshtaghi )_1954_1_1\Ozra.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Janzadeh_Mahtab  28Y (Dr.Moshtaghi)_1983_8_23\Janzadeh_Mahtab.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Hami_Zahra 35Y (dr.Ahmadian)_1976_5_9\hami_zahra.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Khademi_Sareieh 47 Y(Dr Nezhati )_1963_6_29\sareieh.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Kiaeifar_Arezoo 38Y (dr.Sigaroodi)_1973_9_1\arezoo.nii');
%  IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Mirtalebi_Aghil 21Y  (Dr.Sigaroodi)_1991_6_17\aghil.nii');
% IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\deylamidata\Nemani_Shirin  9Y  (Dr.Sigaroodi)_2003_8_12\shirin.nii');
IMAGE=load_nii('C:\Users\meisam\Documents\Master_Project\DATA\lesion\romina keshavarz\romania.nii');

IMAGE.img = IMAGE.img + (1.5*max(max(max(IMAGE.img)))) * int16((PI / 255));

% VVVV = double(IMAGE.img(200:250,200:250,200:250));
% viewer3d(VVVV);
% VVVVV= render(VVVV);
% save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\VVVV.mat','VVVV');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\Derogar_Mohamad_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\rezaghasemi_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\afsaneh_bakhshandeh(post)_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\mahdavi_rad_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\behnam_yousefi_Sp.nii');

save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\romania_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\shirin_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\aghil_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\arezoo_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\Janzadeh_Mahtab_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\zahra_zad_Sp.nii');
% save_nii(IMAGE,'C:\Users\meisam\Documents\Master_Project\Symmetry Planes\Ozra_Sp.nii');
% save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\afsaneh_bakhshandeh(post)_Sp.mat');
%    save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\mahdavi_rad_Sp.mat');
%    save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\behnam_yousefi_Sp.mat');

% save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\zahra_zad_Sp_Sp.mat');
% save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\Ozra_Sp.mat');
% save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\Janzadeh_Mahtab_Sp.mat');
% save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\arezoo_Sp.mat');
% save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\aghil_Sp.mat');
% save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\shirin_Sp.mat');
save('C:\Users\meisam\Documents\Master_Project\Symmetry Planes\romania_Sp.mat');