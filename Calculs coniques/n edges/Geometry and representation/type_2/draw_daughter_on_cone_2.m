function [ VF VFd] = draw_daughter_on_cone_2( V, A, p, gen )
%DRAW_ON_CONE_2 Summary of this function goes here
%   This is the version of draw_daughter_on_cone for type 2 solutions. The
%   only difference lay in the drawing of the daughter cell (line 35&36)
%
%   Copyright Antoine Lizée July 2011 
%   @Dumais lab, OEB dpt, Harvard University.

if nargin==3
    gen=1;
end

n=length(A); % number of actual edges (real number of edges + 1)
np=p*1000; % number of points in the drawing
npe=ceil(np/n); % number of points per edge
Vf=zeros(np,2);
gamma=abs(V(1))*2;

% From polar to cartesian coordinates :
[Vx Vy]=pol2cart(V(:,1),-V(:,2)); %We put the gamma 'forbiden' angle of the cone on the left.
V=[Vx Vy];
% Computation of the set of points :
for ii=1:n 
    Vf(npe*(ii-1)+1:npe*ii,:)=findP(V(ii:(ii+1),:),A(ii),linspace(0,1,npe));
end
% From cartesian to polar coordinates :
[Vftheta, Vfrho]=cart2pol(Vf(:,1),Vf(:,2)); 

% From the plan to the cone :
r=pi/(pi-gamma/2);
Vftheta=r*Vftheta ;
Vfrho=Vfrho/r;

% In 3D cartesian coordinates
[VFx VFy]=pol2cart(Vftheta,Vfrho);
VF=[VFx VFy -Vfrho*sqrt(r^2-1)];

% For the daughter, assuming the area for the mother is one.
solution=fixedpoint_2(n-1,gamma*180/pi,1);
Vfdtheta=Vftheta-gen*(solution.theta1+solution.thetan)*pi/180*r;
Vfdrho=Vfrho/(sqrt(2)^gen);
[VFdx VFdy]=pol2cart(Vfdtheta,Vfdrho);
VFd=[VFdx VFdy -Vfdrho*sqrt(r^2-1)];


end

