%  Function version of the script for profiling
%
%   Copyright Antoine Lizée July 2011 
%   @Dumais lab, OEB dpt, Harvard University.

function [h1 h2]=test_3D_function_2(gamma, n, h1, h2)

[ V A, ~ ]=findfixedpoint(n,gamma,1);
[VF VFd]=draw_daughter_on_cone_2(V,A,1);

figure(h1);
clf
hold on
draw_2D(V,A);
plot(VF(:,1),VF(:,2));
axis equal


figure(h2);
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
cell_cone=[[linspace(0,0.3,58)';linspace(0.3,0.5,6)'],[linspace(0.5,0.7,58)';linspace(0.7,0.8,6)'], 0.1*ones(64,1)];
a=Z(2,floor(size(Z,2)/2));
bool=Z<a;
Z(bool)=a;
h=surf(X,Y,Z);


%3D Plot edition
colormap(cell_cone); %define the custom colormap as the colormap ued by the entire figure
shading interp %Define the way the Surface is shaded. It sets the parameters of Edgecolor and Facecolor altogether of on axes object. interp means interpolation, and enables each of the colors to be smooth.
%set(h,'EdgeColor','None') % That is already the case !
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
%plot(VF(:,1),VF(:,2));

%Plotting the daughter cell :
plot3(VFd(:,1),VFd(:,2),VFd(:,3),'r','LineWidth',1.5);
%plot(VFd(:,1),VFd(:,2));


end

