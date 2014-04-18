% clear all
clearvars
clearvars -GLOBAL
close all
format shorte
set(0,'DefaultFigureWindowStyle','docked')
global C V Mun Mup Gv Dn Dp Bv Em DnM MunM DpM MupM np pp x xm n p
global Rho divFp divFn niSi TwoCarriers t

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb*2*pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.Mun_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

EpiSi = C.eps_0*11.68;
MunSi = 1400*1e-4; % cm2 V-1s-1 * 1/(100 cm/m)^2
DnSi = MunSi*C.kb*300/C.q_0; % D = kt/q Mun
MupSi = 450*1e-4; % cm2 V-1s-1 * 1/(100 cm/m)^2
DpSi = MunSi*C.kb*300/C.q_0; % D = kt/q Mun
tauSi = 1e-8;

niSi = 1e10*1e6; % 1/cm^3 * (100 cm/m)^3 intrinsic concentration

Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 201;
l = 1e-6;
x =linspace(0,l,nx);

dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

% Poisson equation d^2V/dx^2 = -1/EpiSi rho(x)
% Gv V = -dx^2/EpiSi rho(x)
% E = - dV/dx


FormGv(nx,0);
nV = zeros(1,nx);
[L,U] = lu(Gv);

Mun = ones(1,nx)*MunSi;
Dn = ones(1,nx)*DnSi;
MunM(1:nx-1) = (Mun(1:nx-1) + Mun(2:nx))/2;
DnM(1:nx-1) = (Dn(1:nx-1) + Dn(2:nx))/2;
n = zeros(1,nx);

Mup = ones(1,nx)*MupSi;
Dp = ones(1,nx)*DnSi;
MupM(1:nx-1) = (Mup(1:nx-1) + Mup(2:nx))/2;
DpM(1:nx-1) = (Dp(1:nx-1) + Dp(2:nx))/2;
p = zeros(1,nx);

Nd = 1e16 * 1e6; % Const. 1/cm3 (100 cm/m)^3
NetDoping = ones(1,nx).*Nd; % doping


if TwoCarriers == 1
    n0 = (NetDoping + sqrt(NetDoping.^2 + 4* niSi*niSi))/2;
    p0 = niSi^2./n0;
else
    n0  = NetDoping;
    p0 = zeros(1,nx);
end

x0 = l/2;
nw = l/20;
n0 = n0 + 10e16*1e6*exp(-((x-x0)/nw).^2);

if TwoCarriers == 1
    p0 = p0 + 10e16*1e6*exp(-((x-x0)/nw).^2);
end

divFn = zeros(1,nx);
divFp = zeros(1,nx);

dtMax = min(dx^2/2/max(Dn),dx^2/2/max(Dp));

Rho = zeros(1,nx);
if (Coupled)
    Rho = C.q_0*(NetDoping - n0 + p0); % update Rho
    Rho(1) = 0;
end

Rho(1) = 0;
LVbc = 0;
Bv(1) = LVbc;

V = U\(L\(-dx^2/EpiSi*Rho' + Bv'));
Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;
MaxEm = max(abs(Em));
Maxn = max(n0);
Ld = sqrt(EpiSi/(C.q_0*Maxn));
dxMax = Ld/5;

if MaxEm > 0
    dt = min([2*dx/MaxEm dtMax])/4;
else
    dt = dtMax/4
end

t = 0;
n = n0;
np = n0;

p = p0;
pp = p0;

PlotVals(nx,dx,'on',[]);

TStop = 1200000*dt;
PlDelt = 1000*dt;
Plt0 = PlDelt;

while t < TStop

    MaxEm = max(abs(Em));
    
    if MaxEm > 0
        dt = min([2*dx/MaxEm dtMax])/4;
    else
        dt = dtMax/4;
    end
    
    V = U\(L\(-dx^2/EpiSi*Rho' + Bv'));
    Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;
    
    DivFn(nx,dx);
    n = np + dt*divFn;
    Maxn = max(n);
    
    if TwoCarriers == 1
        DivFp(nx,dx);
        p = pp + dt*divFp;
        Maxp = max(p);
        if RC
            Urc = (n.*p - n0.*p0)./(n + p + 2*niSi)/tauSi;
            n = n - dt*Urc;
            p = p - dt*Urc;
        end
    else
        Maxp = 0;
    end
    
    Ldn = sqrt(EpiSi/(C.q_0*Maxn));
    
    if TwoCarriers == 1
        Ldp = sqrt(EpiSi/(C.q_0*Maxp));
    else
        Ldp = 1e6;
    end
    
    dxMax = min(Ldn,Ldp)/5;
    
    if (Coupled)
        Rho = C.q_0*(NetDoping - n + p); % update Rho
        Rho(1) = 0;
    end
    
    
    np = n;
    pp = p;
    
    t = t + dt;
    if t > Plt0
        fprintf('time: %g (%5.2g %%)\n',t,t/TStop*100);
        PlotVals(nx,dx,'off',[]);
        Plt0 = Plt0 + PlDelt;
        pause(0.0001)
    end
    
    
end

figure
PlotVals(nx,dx,'off',[]);