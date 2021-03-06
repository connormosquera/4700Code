clearvars
clearvars -GLOBAL
close all

W = 100;
L = 3/2*W;
Vo = 1;
Terms = 60;

V = zeros(L+1,W);

for k=1:Terms
    for j=1:W
        for i=1:L+1
            V(i,j) = V(i,j) + 4*Vo/pi*1/(2*k-1)*cosh((2*k-1)*pi*(i-L/2-1)/W)/cosh((2*k-1)*pi*(L/2)/W)*sin((2*k-1)*pi*j/W);
        end
    end
    surf(V);
    pause(0.05);
end
