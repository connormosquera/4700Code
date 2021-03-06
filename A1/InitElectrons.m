function [x,y,Vx,Vy] = InitElectrons()
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)

%------------- BEGIN CODE --------------
global nElectrons T L W x y Vx Vy C Vth sigmaMB

x = rand(1, nElectrons)*L; % assigning random initial particle positions
y = rand(1, nElectrons)*W;

for i=1:nElectrons
   while(1) 
      if ( x(i)<1.2e-7 && x(i)>0.8e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)) 
          x(i) = rand*L;
          y(i) = rand*W;
      else
        break
      end
   end
end

Theta = rand(1, nElectrons)*2*pi;
Vx = cos(Theta).*(Vth + sigmaMB*randn(1, nElectrons));
Vy = sin(Theta).*(Vth + sigmaMB*randn(1, nElectrons));

%figure(1)
%histogram(sqrt(Vx.^2+Vy.^2), 20);
%hold off

avgV = sum(sqrt(Vx.^2+Vy.^2))/nElectrons



%------------- END OF CODE --------------
%Please send suggestions for improvement of the above template header 
%to Denis Gilbert at this email address: gilbertd@dfo-mpo.gc.ca.
%Your contribution towards improving this template will be acknowledged in
%the "Changes" section of the TEMPLATE_HEADER web page on the Matlab
%Central File Exchange