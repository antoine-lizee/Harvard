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
hold on
draw_2D(V,A);
plot(VF(:,1),VF(:,2));
axis equal

figure(102);
clf
hold on

% Creation of the cone
a=(2*pi-gamma*pi/180)/(2*pi);
h=sqrt(1/a^2-1);
bd=1.3/sqrt(1+h);
[X Y]=meshgrid(-bd:0.05:bd,-bd:0.05:bd);
Rho=sqrt(X.^2+Y.^2);
d=0.05;
Z=-Rho*h-h*d*exp(-Rho/d);
h=surf(X,Y,Z);
Cmap=[
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.4000         0    0.4000
    0.5765    0.0039    0.5765
    0.6173    0.0031    0.6173
    0.6580    0.0024    0.6580
    0.6988    0.0016    0.6988
    0.7396    0.0008    0.7396
    0.7804         0    0.7804];


%3D Plot edition
colormap(Cmap); %define the custom colormap as the colormap ued by the entire figure
shading interp %Define the way the Surface is shaded. It sets the parameters of Edgecolor and Facecolor altogether of on axes object. interp means interpolation, and enables each of the colors to be smooth.
%set(h,'EdgeColor','None') % That is already the case !
%set(h,'BackFaceLighting','unlit') %Default value
lighting phong % Define the way of computing the lightning effects. 'phong' is an interpolation. Try flat to understand... This function set the FaceLighting and EdgeLighting properties of surfaces and patches 
l=light('Position',[-2,0,2]);
%delete(findobj('Type','light')); %Useful line !
material([0.5, 0.9, 0.3]); %Set the reflectance properties of the surface
alpha(0.9); % Set the transparency
view([-27 24]);
axis equal

set(gca,'Ztick',[],'XTick',[],'YTick',[]);

%Plotting the mother cell :
plot3(VF(:,1),VF(:,2),VF(:,3),'LineWidth',2,'Color',[1 0.776 1]);
% plot(VF(:,1),VF(:,2));

%Plotting the daughter cell :
plot3(VFd(:,1),VFd(:,2),VFd(:,3),'r','LineWidth',1.5,'Color',[1 0.8 1]);
% plot(VFd(:,1),VFd(:,2));

