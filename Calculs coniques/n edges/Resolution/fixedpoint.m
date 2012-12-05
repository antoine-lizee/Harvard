function [ solution values  ] = fixedpoint( n, gamma, H1 )
%FIXEDPOINT Summary of this function goes here
%This function finds the different values that characterize the solution of
%the fixed point problem on the cone depending on the number "n" of edges
%of the solution, and the angle "gamma" of the cone. 
%The normalizing parameter
%"H1" is only needed for the output "values", which compute the numerical
%values needed to draw the solution in a drawing software like geogebra
%(see the script "fixedpointforn").
%
%Copyright Antoine Lizée @ Dumais Laboratory, Harvard University
%July 2011

r=sqrt(2);
r2=r^(n-1);
solution.gamma=gamma;
gamma=gamma/180*pi;

    function [delta]=delta1(theta1)
        delta=1-sqrt((1+r^2)/(1+r^2-2*r*cos(theta1)));
    end

    function [delta]=delta2(theta2)
        delta=1-sqrt((1+r2^2)/(1+r2^2-2*r2*cos(theta2)));
    end

    function [theta2]=theta1_2(theta1)
        theta2=2*pi-(n-1)*theta1-gamma;
    end

theta0=(2*pi-gamma)/n;
options=optimset('Display', 'notify');
[theta1]=fzero(@(theta1) delta1(theta1)-delta2(theta1_2(theta1)), theta0, options);

solution.theta1=theta1*180/pi;
theta2=theta1_2(theta1);
solution.thetan=theta2*180/pi;
solution.alpha=(theta2-theta1)*180/pi;
solution.delta=delta1(theta1);

x=-sign(solution.delta)*(1-solution.delta);

values=zeros(1,4);
values(1)=H1/abs(solution.delta);
values(2)=x*values(1);
values(3)=theta1;
values(4)=theta2; %definition of theta between h_n and h_n+1

end


