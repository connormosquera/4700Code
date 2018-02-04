clearvars
clearvars -GLOBAL
close all

global nElectrons nPlot T L W 
global x y Vx Vy C time Temp dt Vth sigmaMB

C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant

nElectrons = 10000;
nPlot = 20;
T = 300;
L = 200e-9;
W = 100e-9;
MarkerSize = 1;
dt = 1e-15;
TStop = 1e-12;
Vth = sqrt(2*C.kb*T/(C.m_0*0.26));
time = 0;
Temp = T;
taumn = 0.2e-12;
sigmaMB = sqrt(C.kb*T/(C.m_0*0.26));

collisionT = zeros(200,nElectrons);
collisionV = zeros(200,nElectrons);
collisionIndex = ones(1,nElectrons);
collisions = 0;

InitElectrons;

cc = jet(nElectrons);

figure(1)
hold all
plot([0.8,0.8]*1e-7,[0,0.4]*1e-7, 'r-')
plot([0.8,0.8]*1e-7,[0.6,1]*1e-7, 'r-')
plot([1.2,1.2]*1e-7,[0,0.4]*1e-7, 'r-')
plot([1.2,1.2]*1e-7,[0.6,1]*1e-7, 'r-')
plot([0.8,1.2]*1e-7,[0,0]*1e-7, 'r-')
plot([0.8,1.2]*1e-7,[0.4,0.4]*1e-7, 'r-')
plot([0.8,1.2]*1e-7,[0.6,0.6]*1e-7, 'r-')
plot([0.8,1.2]*1e-7,[1,1]*1e-7, 'r-')


for i=0:dt:TStop
    time = i;
    %PlotAll;
    PlotElectrons(cc);
    
    TempCalc();
    
    x = x - dt * Vx; % moving the particles in one time step
    y = y - dt * Vy;
    
    for j=1:nElectrons
        if x(j) > L
            x(j) = x(j) - L;
        elseif x(j) < 0
            x(j) = x(j) + L;
        end
        
         if y(j) > W
             Vy(j) = -Vy(j);
         elseif y(j) < 0
             Vy(j) = -Vy(j);
         end
    end
    
    for j=1:nElectrons
        if (1-exp(-dt/taumn)) > rand()
            collisions = collisions+1;
            collisionT(collisionIndex(j)+1,j) = time;
            collisionV(collisionIndex(j)+1,j) = sqrt(Vx(j)^2+Vy(j)^2);
            collisionIndex(j)=collisionIndex(j)+1;
            
            Theta = rand(1, 1)*2*pi;
            Vx(j) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
            Vy(j) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
        end
    end
    
    BlockBorders(); % uncomment for specular reflection
    %BlockBordersDiffusive(); % uncomment for diffuse reflection
    
    pause(0.001)
    
end


MFP=0;
TBC=0;

for i=1:nElectrons
    for j=1:collisionIndex(i)
        if j ~= 1
            TBC = TBC + collisionT(j,i)-collisionT(j-1,i);
            MFP = MFP + (collisionT(j,i)-collisionT(j-1,i))*(collisionV(j,i)-collisionV(j-1,i));
        end
    end
end

TBC = TBC/collisions
MFP = MFP/collisions



figure(2)
hold on
n=hist3([x',y'],[50 50]);
pcolor(n');
colorbar;
title('Electron Density Map');
hold off

V5050 = zeros(50);

for h=1:nElectrons
    for i=1:50
        for j=1:50
            if x(h)>((i-1)/50*L) && x(h)<(i/50*L) && y(h)>((j-1)/50*W) && y(h)<(j/50*W)
                V5050(i,j)=Vx(h)^2+Vy(h)^2;
            end
        end
    end
end

for i=1:50
    for j=1:50
       if n(i,j)~=0
          V5050(i,j) = V5050(i,j)/n(i,j);  
       else
          V5050(i,j) = 0; 
       end
    end
end


figure(3)
hold on
m=V5050.*0.5*0.26*C.m_0/C.kb;
pcolor(m');
colorbar;
title('Temperature Map');
hold off