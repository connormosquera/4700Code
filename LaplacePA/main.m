% Laplace form: (T(m+1,n) - 2*T(m,n) + T(m-1,n))/dx^2 + 
%                       (T(m,n+1) - 2*T(m,n) + T(m,n-1))/dy^2 = 0

clear all

size=50;
nx=size;
ny=size;
V=zeros(nx,ny);
G=zeros(nx,ny);
maxI = 500;
d=1;

V(:,1)=ones(nx,1);
Vtemp=V;

figure(1)
hold on

for h=1:maxI
    for i=2:nx-1
        for j=2:ny-1
            Vtemp(i,j) = 1/4*(V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1));
        end
        Vtemp(1,i) = 1/3*(V(i-1,j) + V(i+1,j) + V(2,i));
        Vtemp(ny,i) = 1/3*(V(i-1,j) + V(i+1,j) + V(ny-1,i));
    end
    V=Vtemp;
    surf(V);
    pause(0.001);
end