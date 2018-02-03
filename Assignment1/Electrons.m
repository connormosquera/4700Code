clearvars
clearvars -GLOBAL
close all

global nElectrons nPlot T L W
global x y Vx Vy C time dt Vth Temp

C.q_0 = 1.60217653e-19;             % electron charge
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant

nElectrons = 10000;
nPlot=20; % number of electrons to actually plot
T = 300;
L = 200e-9;
W = 100e-9;
dt = 1e-15;
TStop = 1e-12; %1000 timesteps
Vth = sqrt(2*C.kb*T/(C.m_0*0.26)); %using 2 degrees of freedom
time = 0;
Temp = T; % temperature variable that updates in TempCalc
taumn = 0.2e-12;

InitElectrons; % initialize positions and velocities

figure(1)
hold all
for i=0:dt:TStop
    time = i;
    TempCalc();
    
    %PlotTemp; % uncomment for temperature vs. time plot
    PlotElectrons(); % uncomment for trajectories plot
    
    x = x - dt * Vx; % moving the particles in one time step
    y = y - dt * Vy;
    
    for j=1:nElectrons % boundary condition definitions
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

    pause(0.001)
    
end