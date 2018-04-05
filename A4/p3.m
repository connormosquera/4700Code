%% Part 3 (DO NOT FORGET FOURIER COMPONENT)
% 
 
G1=1/1;
C=0.25;
G2=1/2;
L=0.2;
G3=1/10;
alpha=100;
G4=1/0.1;
G0=1/1000;
 
%Cn = 0.00001;
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
fs = T/dt;
time = 0:(1/fs):(1-1/fs);
n=fs;
f = (0:n-1);

%%%%%%%%% Cn = 0.00001 %%%%%%%%%%
V = zeros(7,T/dt);
F = zeros(7,T/dt);
Ap = inv(C/dt + G);
 
% Vin defined as a gaussian
F(1,:) = normpdf(time,0.06,0.03)*max(normpdf(time,0.06,0.03))^(-1);
F(4,:) = normrnd(0.001, 0.0003,1,T/dt); % In randomly picked from normal dist
 
for i=2:T/dt
    V(:,i) = Ap*(C*V(:,i-1)/dt + F(:,i));
end
 
figure (1)
plot(time,F(1,:)','LineWidth',2)
hold on
grid on
plot(time,V(7,:)','LineWidth',2)
title('\fontsize{22}Gaussian Pulse, Cn = 0.00001')
xlabel('\fontsize{18}Time (s)')
ylabel('\fontsize{18}Voltage (V)')
legend('\fontsize{18}Vin','\fontsize{18}V0')

XF1 = (fft(F(1,:))).^2/n;
YF1 = fftshift(XF1);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftF1 = abs(YF1).^2/n;     % zero-centered power

XV7 = (fft(V(7,:))).^2/n;
YV7 = fftshift(XV7);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftV7 = abs(YV7).^2/n;     % zero-centered power

figure (2)
plot(fshift,powershiftF1,'LineWidth',2)
hold on
grid on
plot(fshift,powershiftV7,'LineWidth',2)
title('\fontsize{22}Frequency Composition of Gaussian Pulse, Cn = 0.00001')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-20 20])

%%%%%%%%% Cn = 0.1 %%%%%%%%%%
Cn = 0.1;
C(4,4) = Cn; % updating matrix

V = zeros(7,T/dt);
F = zeros(7,T/dt);
Ap = inv(C/dt + G);
 
% Vin defined as a gaussian
F(1,:) = normpdf(time,0.06,0.03)*max(normpdf(time,0.06,0.03))^(-1);
F(4,:) = normrnd(0.001, 0.0003,1,T/dt); % In randomly picked from normal dist
 
for i=2:T/dt
    V(:,i) = Ap*(C*V(:,i-1)/dt + F(:,i));
end
 
figure (3)
plot(time,F(1,:)','LineWidth',2)
hold on
grid on
plot(time,V(7,:)','LineWidth',2)
title('\fontsize{22}Gaussian Pulse, Cn = 0.1')
xlabel('\fontsize{18}Time (s)')
ylabel('\fontsize{18}Voltage (V)')
legend('\fontsize{18}Vin','\fontsize{18}V0')

XF1 = (fft(F(1,:))).^2/n;
YF1 = fftshift(XF1);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftF1 = abs(YF1).^2/n;     % zero-centered power

XV7 = (fft(V(7,:))).^2/n;
YV7 = fftshift(XV7);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftV7 = abs(YV7).^2/n;     % zero-centered power

figure (4)
plot(fshift,powershiftF1,'LineWidth',2)
hold on
grid on
plot(fshift,powershiftV7,'LineWidth',2)
title('\fontsize{22}Frequency Composition of Gaussian Pulse, Cn = 0.1')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-20 20])

%%%%%%%%% Cn = 0.0000001 %%%%%%%%%%
Cn = 0.0000001;
C(4,4) = Cn; % updating matrix

V = zeros(7,T/dt);
F = zeros(7,T/dt);
Ap = inv(C/dt + G);
 
% Vin defined as a gaussian
F(1,:) = normpdf(time,0.06,0.03)*max(normpdf(time,0.06,0.03))^(-1);
F(4,:) = normrnd(0.001, 0.0003,1,T/dt); % In randomly picked from normal dist
 
for i=2:T/dt
    V(:,i) = Ap*(C*V(:,i-1)/dt + F(:,i));
end
 
figure (5)
plot(time,F(1,:)','LineWidth',2)
hold on
grid on
plot(time,V(7,:)','LineWidth',2)
title('\fontsize{22}Gaussian Pulse, Cn = 0.0000001')
xlabel('\fontsize{18}Time (s)')
ylabel('\fontsize{18}Voltage (V)')
legend('\fontsize{18}Vin','\fontsize{18}V0')

XF1 = (fft(F(1,:))).^2/n;
YF1 = fftshift(XF1);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftF1 = abs(YF1).^2/n;     % zero-centered power

XV7 = (fft(V(7,:))).^2/n;
YV7 = fftshift(XV7);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftV7 = abs(YV7).^2/n;     % zero-centered power

figure (6)
plot(fshift,powershiftF1,'LineWidth',2)
hold on
grid on
plot(fshift,powershiftV7,'LineWidth',2)
title('\fontsize{22}Frequency Composition of Gaussian Pulse, Cn = 0.0000001')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-20 20])

%%%%%%%%% dt = 25e-3 %%%%%%%%%%
Cn = 0.0001;
C(4,4) = Cn; % updating matrix

dt = 25e-3;
T = 1;
fs = T/dt;
time = 0:(1/fs):(1-1/fs);
n=fs;
f = (0:n-1);

V = zeros(7,T/dt);
F = zeros(7,T/dt);
Ap = inv(C/dt + G);
 
% Vin defined as a gaussian
F(1,:) = normpdf(time,0.06,0.03)*max(normpdf(time,0.06,0.03))^(-1);
F(4,:) = normrnd(0.001, 0.0003,1,T/dt); % In randomly picked from normal dist
 
for i=2:T/dt
    V(:,i) = Ap*(C*V(:,i-1)/dt + F(:,i));
end
 
figure (7)
plot(time,F(1,:)','LineWidth',2)
hold on
grid on
plot(time,V(7,:)','LineWidth',2)
title('\fontsize{22}Gaussian Pulse, dt=25e-3')
xlabel('\fontsize{18}Time (s)')
ylabel('\fontsize{18}Voltage (V)')
legend('\fontsize{18}Vin','\fontsize{18}V0')

XF1 = (fft(F(1,:))).^2/n;
YF1 = fftshift(XF1);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftF1 = abs(YF1).^2/n;     % zero-centered power

XV7 = (fft(V(7,:))).^2/n;
YV7 = fftshift(XV7);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftV7 = abs(YV7).^2/n;     % zero-centered power

figure (8)
plot(fshift,powershiftF1,'LineWidth',2)
hold on
grid on
plot(fshift,powershiftV7,'LineWidth',2)
title('\fontsize{22}Frequency Composition of Gaussian Pulse, dt=20e-3')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-20 20])