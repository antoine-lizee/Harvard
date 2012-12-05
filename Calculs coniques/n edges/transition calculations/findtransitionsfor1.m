%cone.m
%
% This program finds some of the simple solution for dividing a cone

figure(10)

res = 100;
R = 1;
beta = [0:2*pi/(res-1):2*pi];
l1 = (2*pi*ones(1,res) - beta)*R/(2)^(1/2);
l2 = 2*R*ones(1,res);

gamma = [0:0.95*pi/(res-1):0.95*pi];
beta2 = 2*pi*ones(1,res) - 2*(gamma-sin(gamma)+(pi*ones(1,res)-gamma-sin(gamma)).*(tan(gamma/2).^2));
l3 = R*(pi*ones(1,res)-gamma).*tan(gamma/2);

plot(beta*180/pi,l1,'g',beta*180/pi,l2,'b',beta2*180/pi,l3,'r')