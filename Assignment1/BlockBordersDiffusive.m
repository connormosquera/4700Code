function [ ] = BlockBordersDiffusive()
global nElectrons T L W MarkerSize
global x y Vx Vy dt sigmaMB Vth

for i=1:nElectrons % rethermalized when hit boundary, theta defines scattering angle so it reflects away from boundary
    if Vy(i)<0 && y(i)>0.6e-7 && y(i)<0.61e-7 && x(i)<1.2e-7 && x(i)>0.8e-7
        Theta = rand*pi;
        Vx(i) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
        Vy(i) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
    elseif Vy(i)>0 && y(i)<0.4e-7 && y(i)>0.39e-7 && x(i)<1.2e-7 && x(i)>0.8e-7
        Theta = rand*pi+pi;
        Vx(i) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
        Vy(i) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
    elseif Vx(i)<0 && x(i)>0.8e-7 && x(i)<0.81e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)
        Theta = rand*pi-pi/2;
        Vx(i) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
        Vy(i) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
    elseif Vx(i)>0 && x(i)<1.2e-7 && x(i)>1.19e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)
        Theta = rand*pi+pi/2;
        Vx(i) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
        Vy(i) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
    end
end

end