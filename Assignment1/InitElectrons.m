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

Theta = rand(1, nElectrons)*2*pi;
Vx = cos(Theta).*(Vth + sigmaMB*randn(1, nElectrons));
Vy = sin(Theta).*(Vth + sigmaMB*randn(1, nElectrons));

% figure(1)
% histogram(sqrt(Vx.^2+Vy.^2), 20);
% hold on
% xlabel('Binned velocities (m/s)');
% ylabel('Frequency');
% grid on
% title('Histogram of velocities of 10000 electrons');
% hold off

avgV = sum(sqrt(Vx.^2+Vy.^2))/nElectrons



%------------- END OF CODE --------------
%Please send suggestions for improvement of the above template header 
%to Denis Gilbert at this email address: gilbertd@dfo-mpo.gc.ca.
%Your contribution towards improving this template will be acknowledged in
%the "Changes" section of the TEMPLATE_HEADER web page on the Matlab
%Central File Exchange