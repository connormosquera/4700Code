%%%%% Part 1 %%%%%

Is = 0.01e-12;
Ib = 0.1e-12;
Vb = 1.3;
Gp = 0.1;

Ifunc = @(V) Is*(exp(1.2*V/0.025)-1) + Gp*V - Ib*exp(-1.2*(V+Vb)/0.025);

V = linspace(-1.95, 0.7, 200);
I = Ifunc(V);
In = I.*(rand(1,200)*2/5 + 0.8);

%%%%% Part 2 %%%%%

P4 = polyfit(V,I,4);
P4n = polyfit(V,In,4);
P8 = polyfit(V,I,8);
P8n = polyfit(V,In,8);

IP4 = polyval(P4,V);
IP4n = polyval(P4n,V);
IP8 = polyval(P8,V);
IP8n = polyval(P8n,V);

% figure(1)
% plot(V,I)
% hold on
% plot(V,IP4)
% plot(V,IP8)
% grid on
% title('Linear Plot of V vs. I');
% legend('I Data','4th order fit','8th order fit')
% hold off
% 
% figure(2)
% semilogy(V,abs(I))
% hold on
% semilogy(V,abs(IP4))
% semilogy(V,abs(IP8))
% grid on
% title('Log Plot of V vs. I');
% legend('I Data','4th order fit','8th order fit')
% hold off
% 
% figure(3)
% plot(V,In)
% hold on
% plot(V,IP4n)
% plot(V,IP8n)
% grid on
% title('Linear Plot of V vs. In');
% legend('I Data','4th order fit','8th order fit')
% hold off
% 
% figure(4)
% semilogy(V,abs(In))
% hold on
% semilogy(V,abs(IP4n))
% semilogy(V,abs(IP8n))
% grid on
% title('Log Plot of V vs. In');
% legend('I Data','4th order fit','8th order fit')
% hold off

%%%%% Part 3 %%%%%

fo = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+D))/25e-3)-1)');
ff = fit(V',I',fo);
If = ff(V)';

f1 = fittype('A.*(exp(1.2*x/25e-3)-1) + 0.1.*x - C*(exp(1.2*(-(x+1.3))/25e-3)-1)');
ff1 = fit(V',I',f1);
If1 = ff1(V)';

f2 = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+1.3))/25e-3)-1)');
ff2 = fit(V',I',f2);
If2 = ff2(V)';

f3 = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+D))/25e-3)-1)');
ff3 = fit(V',I',f3);
If3 = ff3(V)';

f1n = fittype('A.*(exp(1.2*x/25e-3)-1) + 0.1.*x - C*(exp(1.2*(-(x+1.3))/25e-3)-1)');
ff1n = fit(V',I',f1n);
If1n = ff1n(V)';

f2n = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+1.3))/25e-3)-1)');
ff2n = fit(V',I',f2n);
If2n = ff2n(V)';

f3n = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+D))/25e-3)-1)');
ff3n = fit(V',I',f3n);
If3n = ff3n(V)';

% figure(5)
% plot(V,I)
% hold on
% plot(V,If1)
% plot(V,If2)
% plot(V,If3)
% grid on
% title('Linear Plot of V vs. I for 3 different "fit" types');
% legend('I Data','Fit for A and C','Fit for A, B, and C','Fit for A, B, C, and D')
% hold off
% 
% figure(6)
% semilogy(V,abs(I))
% hold on
% semilogy(V,abs(If1))
% semilogy(V,abs(If2))
% semilogy(V,abs(If3))
% grid on
% title('Log Plot of V vs. I for 3 different "fit" types');
% legend('I Data','Fit for A and C','Fit for A, B, and C','Fit for A, B, C, and D')
% hold off
% 
% figure(7)
% plot(V,In)
% hold on
% plot(V,If1n)
% plot(V,If2n)
% plot(V,If3n)
% grid on
% title('Linear Plot of V vs. In for 3 different "fit" types');
% legend('I Data','Fit for A and C','Fit for A, B, and C','Fit for A, B, C, and D')
% hold off
% 
% figure(8)
% semilogy(V,abs(In))
% hold on
% semilogy(V,abs(If1n))
% semilogy(V,abs(If2n))
% semilogy(V,abs(If3n))
% grid on
% title('Log Plot of V vs. In for 3 different "fit" types');
% legend('I Data','Fit for A and C','Fit for A, B, and C','Fit for A, B, C, and D')
% hold off

%%%%% Part 4 %%%%%

% inputs = V;
% targets = I;
% hiddenLayerSize = 10;
% net = fitnet(hiddenLayerSize);
% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 15/100;
% [net,tr] = train(net,inputs,targets);
% outputs = net(inputs);
% errors = gsubtract(outputs,targets);
% performance = perform(net,targets,outputs)
% view(net)
% Inn = outputs

inputs = V;
targets = In;
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net,tr] = train(net,inputs,targets);
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs)
view(net)
Inn = outputs

