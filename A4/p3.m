%% Part 2 (DO NOT FORGET FOURIER COMPONENT)
% a) linear circuit
% 

G1=1/1;
C=0.25;
G2=1/2;
L=0.2;
G3=1/10;
alpha=100;
G4=1/0.1;
G0=1/1000;

Cn = 0.00001;

% V = [V1
%       V2
%       IL
%       V3
%       I3
%       V4
%       V0]

% F = [Vin
%       0
%       0
%       In
%       0
%       0
%       0]

G = [1 0 0 0 0 0 0;
    -G2 G1+G2 -1 0 0 0 0;
    0 1 0 -1 0 0 0;
    0 0 -1 G3 0 0 0;
    0 0 0 0 -alpha 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+G0];

C = [0 0 0 0 0 0 0;
    -C C 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 Cn 0 0 0; % Cn added to C matrix
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];

dt = 1e-3;
T = 1;
time = linspace(dt,T,T/dt);
V = zeros(7,T/dt);
F = zeros(7,T/dt);
Ap = inv(C/dt + G);

In = normrnd(0.001, 0.0003); % In randomly picked from normal dist


%%% 2d i)C %%%
% Vin defined as a gaussian
V = zeros(7,T/dt);
F = zeros(7,T/dt);
F = normpdf(time,0.06,0.03)*max(normpdf(time,0.06,0.03))^(-1);

for i=2:T/dt
    V(:,i) = Ap*(C*V(:,i-1)/dt + F(:,i));
end

% figure (5)
% plot(time,F(1,:)','LineWidth',2)
% hold on
% grid on
% plot(time,V(7,:)','LineWidth',2)
% title('\fontsize{22}Gaussian Pulse')
% xlabel('\fontsize{18}Time (s)')
% ylabel('\fontsize{18}Voltage (V)')
% legend('\fontsize{18}Vin','\fontsize{18}V0')


