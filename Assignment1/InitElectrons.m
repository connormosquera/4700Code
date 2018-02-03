function [x,y,Vx,Vy] = InitElectrons()
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)

%------------- BEGIN CODE --------------
global nElectrons T L W x y Vx Vy C Vth

x = rand(1, nElectrons)*L; % assigning random initial particle positions
y = rand(1, nElectrons)*W;

Theta = rand(1, nElectrons)*2*pi; % random velocity direction
Vx = cos(Theta)*Vth;
Vy = sin(Theta)*Vth;

