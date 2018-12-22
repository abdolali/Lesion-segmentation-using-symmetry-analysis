all_median = mean(X(:));
all_std = std(X(:));
ii=1; jj=1;
for i=17:495
      if (i ~= 66 && i ~= 69 && i ~= 442 && i ~= 445 && i ~= 454 && i ~= 452 && i ~= 447 && i ~= 483 && i ~= 493 && i ~= 457 && i ~= 469 &&  i ~= 482 &&  i ~= 488 &&  i ~= 494 && i ~= 479 )
            for j = 1:x_di
                  Y(ii , j) = X(i ,j); 
            end
            ii=ii+1;
      end
end
ii=ii-1;
xx=1:x_di;       xx=repmat(xx,[1,ii]);
zz=zeros(size(xx));
yy=zeros(size(xx));
for i=1:ii
    zz(x_di*(i-1)+1:x_di*(i-1)+x_di)=i;
    yy(x_di*(i-1)+1:x_di*(i-1)+x_di)=Y(i ,:);
end

sf = fit([xx',zz'],yy','poly11');
plot(sf,[xx',zz'],yy')