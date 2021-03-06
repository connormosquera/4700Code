clearvars
clearvars -GLOBAL
close all

L = 100;
W = 100;
Vo = 1;

G = sparse(L*W,L*W);
B = zeros(1,L*W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        
        if i==L 
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = Vo;
        elseif i==1 
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==1
            G(n,:) = 0;
            G(n,n) = -1;
            G(n,n+1) = 1;     % dV/dy = 0      
        elseif j==W
            G(n,:) = 0;
            G(n,n) = -1; 
            G(n,n-1) = 1;     % dV/dy = 0    
        else
            G(n,:) = 0;
            G(n,n) = -4;
            G(n,n-1) = 1;
            G(n,n+1) = 1;
            G(n,n-W) = 1;
            G(n,n+W) = 1;
        end
        
    end
end

spy(G);

F=G\B';

Vmap = zeros(L,W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        Vmap(i,j) = F(n);
    end
end


surf(Vmap);
pause(0.001);
