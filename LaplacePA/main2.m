clear all

size=50;
nx=size;
ny=size;
V=zeros(nx,ny);
G=zeros(nx,ny);
maxI = 500;
d=1;
x=linspace(1,nx,nx);
y=linspace(1,ny,ny);

V(:,1)=ones(nx,1);
V(:,nx)=ones(nx,1);
Vtemp=V;

figure(1)

for h=1:maxI
    for i=2:nx-1
        for j=2:ny-1
            Vtemp(i,j) = 1/4*(V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1));
        end
    end
    
    V=Vtemp;
    [Ex,Ey]=gradient(V);
        
    surf(V);
    %quiver(x,y,Ex,Ey);
    pause(0.001);
end