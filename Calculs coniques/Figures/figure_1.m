%   This script create the final plot of a solution for n edges on a cone
%   of "forbidden angle" gamma. A reduced version of this code is found in
%   the gui script 'solution_on_the_code.

%   Copyright Antoine Lizée July 2011 
%   @Dumais lab, OEB dpt, Harvard University.

gamma=120;
n=3;

[ V A, ~ ]=findfixedpoint(n,gamma,1);
[VF VFd]=draw_daughter_on_cone(V,A,1);

figure(101);
clf
subplot(1,2,1);
hold on
h=draw_2D(V,A);
axis([-1 1.5 -1 1]);
axis fill;
set(h.edges,'Linewidth',2,'color','b');

subplot(1,2,2);
hold on

% Creation of the cone
[X Y]=meshgrid(-1:0.05:1,-1:0.05:1);
Rho=sqrt(X.^2+Y.^2);
a=(2*pi-gamma*pi/180)/(2*pi);
t=sqrt(1/a^2-1);
d=0.05;
Z=-Rho*t;
h=surf(X,Y,Z);

%3D Plot edition
set(h,'EdgeColor','None','facecolor', 'g') ;
%set(h,'BackFaceLighting','unlit') %Default value
lighting phong % Define the way of computing the lightning effects. 'phong' is an interpolation. Try flat to understand... This function set the FaceLighting and EdgeLighting properties of surfaces and patches 
l=light('Position',[-2,0,2]);
%delete(findobj('Type','light')); %Useful line !
material([0.5, 0.9, 0.3]); %Set the reflectance properties of the surface
alpha(0.7); % Set the transparency
view([-27 24]);
axis equal

%Plotting the mother cell :
plot3(VF(:,1),VF(:,2),VF(:,3),'LineWidth',2);

