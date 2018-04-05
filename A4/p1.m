G1=1/1;
C=0.25;
G2=1/2;
L=0.2;
G3=1/10;
alpha=100;
G4=1/0.1;
G0=1/1000;

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

% V = [V1
%       V2
%       IL
%       V3
%       I3
%       V4
%       V0]

F = [0; 0; 0; 0; 0; 0; 0];

%%% DC SWEEP %%%

vinvec = zeros(1,21); % vin vector
v0vec = zeros(1,21); % v0 vector
v3vec = zeros(1,21); % v3 vector

for i=1:21 % sweeping vin from -10 to 10
    F = [i-11; 0; 0; 0; 0; 0; 0];
    V = G\F; % DC solution 
    
    vinvec(i) = i-11;
    v0vec(i) =  V(7);
    v3vec(i) = V(4);
end

figure (1)
plot(vinvec,v0vec)
hold on
grid on
plot(vinvec,v3vec)
title('\fontsize{22}DC Vin Sweep')
xlabel('\fontsize{18}Vin (V)')
ylabel('\fontsize{18}Voltage (V)')
legend('\fontsize{18}v0','\fontsize{18}v3')

%%% AC SWEEP %%%

omegavec = logspace(-3,5,30); % 30 log spaced values from 10^-3 to 10^5
v0vec = zeros(1,30);

F = [1; 0; 0; 0; 0; 0; 0];

for i=1:30
    
    V = (G+i*omegavec(i)*C)\F; % AC solution 
    
    v0vec(i) =  V(7);
end

figure (2)
loglog(omegavec,v0vec)
hold on
grid on
title('\fontsize{22}Frequency vs. V0')
xlabel('\fontsize{18}Frequency (rads)')
ylabel('\fontsize{18}V0 (V)')

figure (3)
semilogx(omegavec,20*log10(v0vec))
hold on
grid on
title('\fontsize{22}Frequency vs. Gain (V0/V1)')
xlabel('\fontsize{18}Frequency (rads)')
ylabel('\fontsize{18}Gain (dB)')

%%% PERTURBATIONS %%%
ncaps = 200;

dist = normrnd(.25,0.05,ncaps,1)
v0vec = zeros(1,ncaps);

for i=1:ncaps
    C(2,1) = -dist(i);
    C(2,2) = dist(i);
    
    V = (G+i*pi*C)\F;
    
    v0vec(i) =  V(7);
end

figure (4)
histogram(real(20*log10(v0vec)),15)
hold on
grid on
title('\fontsize{22}Capacitor Variations')
xlabel('\fontsize{18}Gain (dB)')

