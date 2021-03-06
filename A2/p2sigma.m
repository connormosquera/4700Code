clearvars
clearvars -GLOBAL
close all

SigVec = [50 20 10 2 1 0.1 0.01 0.001 0.0001];
steps=length(SigVec);
Curr = zeros(1,steps);

for h=1:steps

    L = 80;
    W = 3/2*L;
    Lb = round(L/3);
    Wb = round(W/3);
    Vo = 1;
    maxI = 200;
    delta = 1;
    
    Sigma = ones(W,L)*SigVec(h);
    
    for i=1:Wb
        for j=round(L/2-Lb/2):round(L/2+Lb/2)
            Sigma(i,j) = 0.01;
        end
    end
    for i=round(W-Wb):W
        for j=round(L/2-Lb/2):round(L/2+Lb/2)
            Sigma(i,j) = 0.01;
        end
    end
    
    G = sparse(L*W,L*W);
    B = zeros(1,L*W);
    
    for i=1:W
        for j=1:L
            n = j + (i-1)*L;
            
            if i==1 && j==L/2
                G(n,:) = 0;
                G(n,n) = 1;
                B(n) = Vo;
            elseif i==W && j==L/2
                G(n,:) = 0;
                G(n,n) = 1;
            elseif j==1 && i==1
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1));
                G(n,n+1) = 1/Sigma(i,j+1);
                G(n,n+L) = 1/Sigma(i+1,j);
            elseif j==L && i==1
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j-1));
                G(n,n-1) = 1/Sigma(i,j-1);
                G(n,n+L) = 1/Sigma(i+1,j);
            elseif j==1 && i==W
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i-1,j) + 1/Sigma(i,j+1));
                G(n,n+1) = 1/Sigma(i,j+1);
                G(n,n-L) = 1/Sigma(i-1,j);
            elseif j==L && i==W
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i-1,j) + 1/Sigma(i,j-1));
                G(n,n-1) = 1/Sigma(i,j-1);
                G(n,n-L) = 1/Sigma(i-1,j);
            elseif j==1
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1) + ...
                    1/Sigma(i-1,j));
                G(n,n+1) = 1/Sigma(i,j+1);
                G(n,n+L) = 1/Sigma(i+1,j);
                G(n,n-L) = 1/Sigma(i-1,j);
            elseif j==L
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j-1) + ...
                    1/Sigma(i-1,j));
                G(n,n-1) = 1/Sigma(i,j-1);
                G(n,n+L) = 1/Sigma(i+1,j);
                G(n,n-L) = 1/Sigma(i-1,j);
            elseif i==1
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1) + ...
                    1/Sigma(i,j-1));
                G(n,n+1) = 1/Sigma(i,j+1);
                G(n,n+L) = 1/Sigma(i+1,j);
                G(n,n-1) = 1/Sigma(i,j-1);
            elseif i==W
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i,j+1) + 1/Sigma(i,j-1) + ...
                    1/Sigma(i-1,j));
                G(n,n-1) = 1/Sigma(i,j-1);
                G(n,n+1) = 1/Sigma(i,j+1);
                G(n,n-L) = 1/Sigma(i-1,j);
            else
                G(n,:) = 0;
                G(n,n) = -(1/Sigma(i,j+1) + 1/Sigma(i,j-1) + ...
                    1/Sigma(i-1,j) + 1/Sigma(i+1,j));
                G(n,n-1) = 1/Sigma(i,j-1);
                G(n,n+1) = 1/Sigma(i,j+1);
                G(n,n-L) = 1/Sigma(i-1,j);
                G(n,n+L) = 1/Sigma(i+1,j);
            end
            
        end
    end
    
    F=G\B';
    
    Vmap = zeros(W,L);
    
    for i=1:W
        for j=1:L
            n = j + (i-1)*L;
            Vmap(i,j) = F(n);
        end
    end
    
    [Ex,Ey]=gradient(-Vmap);
    
    Jx=Sigma.*Ex;
    Jy=Sigma.*Ey;
    
    Cin = sqrt(Jx(1,L/2)^2+Jy(1,L/2)^2);
    Cout = sqrt(Jx(W,L/2)^2+Jy(W,L/2)^2);
    Curr(h) = (Cin+Cout)*0.5;
    h
end

figure(1)
hold on
plot(SigVec, Curr);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title('Conductivity vs. Current Density')
xlabel('Conductivity')
ylabel('Current Density')
hold off
