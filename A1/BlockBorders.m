function [ ] = BlockBorders()
global nElectrons T L W MarkerSize
global x y Vx Vy dt

for i=1:nElectrons % conditions for meeting a boundary, specular reflection by inverting x or y velocity
    if Vy(i)<0 && y(i)>0.6e-7 && y(i)<0.61e-7 && x(i)<1.2e-7 && x(i)>0.8e-7
        Vy(i)=-Vy(i);
    elseif Vy(i)>0 && y(i)<0.4e-7 && y(i)>0.39e-7 && x(i)<1.2e-7 && x(i)>0.8e-7
        Vy(i)=-Vy(i);
    elseif Vx(i)<0 && x(i)>0.8e-7 && x(i)<0.81e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)
        Vx(i)=-Vx(i);
    elseif Vx(i)>0 && x(i)<1.2e-7 && x(i)>1.19e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)
        Vx(i)=-Vx(i);
    end
end

end