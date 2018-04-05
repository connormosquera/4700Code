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
 
% V = [V1
%       V2
%       IL
%       V3
%       I3
%       V4
%       V0]
 
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
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];
 
dt = 1e-3;
T = 1;
fs = T/dt;
time = 0:(1/fs):(1-1/fs);
n=fs;
f = (0:n-1);

V = zeros(7,T/dt);
F = zeros(7,T/dt);
Ap = inv(C/dt + G);
 
%%% 2d i)A %%%
% Vin defined as 0 until 0.03 when source turns on
F(1,0.03/dt:T/dt) = ones([1,T/dt - 0.03/dt+1]);
 
for i=2:T/dt
    V(:,i) = Ap*(C*V(:,i-1)/dt + F(:,i));
end
 
figure (1)
plot(time,F(1,:)','LineWidth',2)
hold on
grid on
plot(time,V(7,:)','LineWidth',2)
title('\fontsize{22}Step Transition Input')
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
title('\fontsize{12}Frequency Composition of Step Transition Input')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-50 50])

%%% 2d i)B %%%
% Vin defined as a sinusoid
f = 1/0.03;
F1 = zeros(7,T/dt);
F2 = zeros(7,T/dt);
V = zeros(7,T/dt);
V1 = zeros(7,T/dt);
V2 = zeros(7,T/dt);
 
for i=dt/dt:T/dt
    F(1,i) = sin(2*pi*f*dt*i);
    F1(1,i) = sin(2*pi*(f*10)*dt*i); % 10 times the frequency
    F2(1,i) = sin(2*pi*(f*0.1)*dt*i); % 1/10 of the frequency
end
 
for i=2:T/dt
    V(:,i) = Ap*(C*V(:,i-1)/dt + F(:,i));
    V1(:,i) = Ap*(C*V1(:,i-1)/dt + F1(:,i));
    V2(:,i) = Ap*(C*V2(:,i-1)/dt + F2(:,i));
end
 
figure (3)
plot(time,F(1,:)','LineWidth',2)
hold on
grid on
plot(time,V(7,:)','LineWidth',2)
title('\fontsize{22}Sinusoid Input, f=1/0.03')
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
title('\fontsize{22}Frequency Composition of Sinusoid Input, f=1/0.03')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-50 50])

figure (5)
plot(time,F1(1,:)','LineWidth',2)
hold on
grid on
plot(time,V1(7,:)','LineWidth',2)
title('\fontsize{22}Sinusoid Input, f=1/0.003')
xlabel('\fontsize{18}Time (s)')
ylabel('\fontsize{18}Voltage (V)')
legend('\fontsize{18}Vin','\fontsize{18}V0')
axis([0 0.2 -1 1]);

XF1 = (fft(F1(1,:))).^2/n;
YF1 = fftshift(XF1);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftF1 = abs(YF1).^2/n;     % zero-centered power

XV7 = (fft(V1(7,:))).^2/n;
YV7 = fftshift(XV7);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftV7 = abs(YV7).^2/n;     % zero-centered power

figure (6)
plot(fshift,powershiftF1,'LineWidth',2)
hold on
grid on
plot(fshift,powershiftV7,'LineWidth',2)
title('\fontsize{22}Frequency Composition of Sinusoid Input, f=1/0.003')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-500 500])

figure (7)
plot(time,F2(1,:)','LineWidth',2)
hold on
grid on
plot(time,V2(7,:)','LineWidth',2)
title('\fontsize{22}Sinusoid Input, f=1/0.3')
xlabel('\fontsize{18}Time (s)')
ylabel('\fontsize{18}Voltage (V)')
legend('\fontsize{18}Vin','\fontsize{18}V0')
 
XF1 = (fft(F2(1,:))).^2/n;
YF1 = fftshift(XF1);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftF1 = abs(YF1).^2/n;     % zero-centered power

XV7 = (fft(V2(7,:))).^2/n;
YV7 = fftshift(XV7);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftV7 = abs(YV7).^2/n;     % zero-centered power

figure (8)
plot(fshift,powershiftF1,'LineWidth',2)
hold on
grid on
plot(fshift,powershiftV7,'LineWidth',2)
title('\fontsize{22}Frequency Composition of Sinusoid Input, f=1/0.3')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-50 50])

%%% 2d i)C %%%
% Vin defined as a gaussian
V = zeros(7,T/dt);
F = zeros(7,T/dt);
F(1,:) = normpdf(time,0.06,0.03)*max(normpdf(time,0.06,0.03))^(-1);
 
for i=2:T/dt
    V(:,i) = Ap*(C*V(:,i-1)/dt + F(:,i));
end
 
figure (9)
plot(time,F(1,:)','LineWidth',2)
hold on
grid on
plot(time,V(7,:)','LineWidth',2)
title('\fontsize{22}Gaussian Pulse')
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

figure (10)
plot(fshift,powershiftF1,'LineWidth',2)
hold on
grid on
plot(fshift,powershiftV7,'LineWidth',2)
title('\fontsize{22}Frequency Composition of \n Gaussian Pulse')
xlabel('\fontsize{18}Frequency (Hz)')
ylabel('\fontsize{18}Amplitude')
legend('\fontsize{18}Vin','\fontsize{18}V0')
xlim([-50 50])
 
