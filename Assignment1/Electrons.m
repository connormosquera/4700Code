clearvars
clearvars -GLOBAL
close all

global nElectrons nPlot T L W MarkerSize
global x y Vx Vy C time Temp dt Vth sigmaMB

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s�
C.am = 1.66053892e-27;

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
            Theta = rand(1, 1)*2*pi;
            Vx(j) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
            Vy(j) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
        end
    end
    
    BlockBorders(); % uncomment for specular reflection
    %BlockBordersDiffusive(); % uncomment for diffuse reflection
    
    pause(0.001)
    
end

figure(2)
hold on
n=hist3([x',y'],[50 50]);
pcolor(n');
colorbar;
hold off

% figure(3)
% hold on
% m=hist3([0.5*C.kb*0.26*C.m_0*Vx'.^2,0.5*C.kb*0.26*C.m_0*Vy'.^2],[50 50]);
% pcolor(m');
% colorbar;
% hold off