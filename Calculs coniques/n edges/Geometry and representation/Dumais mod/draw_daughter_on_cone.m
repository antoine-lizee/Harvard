function [ VF VFd] = draw_daughter_on_cone( V, A, p, turn )
%DRAW_ON_CONE Summary of this function goes here
%   This Function computes the set of points in cartesian coordinates that 
% describes the solution in 3D, on the cone, from the vertex and angles
% which describe the solution in a geometric way, on the cutted plan.
%
% LEt n be the actual number of edges, which is greater by one than the 
% number of final edges of the solution we are looking for.
% The input args must be a n+1-by-2 matrix V and a n vector A describing
% the n edges with curvature in a way that the gamma angle of the cone is
% on the right. Therefore, the forbidden angles are for theta between 
% [-gamma/2, gamma/2]. These inputs are matched with the outputs of the
% function 'findfixedpoint'.
%
%DRAW DAUGHTER ON CONE is an evolution of the draw_on_cone function, which
%give as an additional output the solution reduced in a sqrt(2) ratio, and
%rotated by the appropriate theta angle (see end of the function). It is a
%powerful demonstration of the success of the solution determination.
%
% N.B.: Could be faster with a all-polar version of the transformation
%
%   Copyright Antoine Lizée July 2011 
%   @Dumais lab, OEB dpt, Harvard University.

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
solution=fixedpoint(n-1,gamma*180/pi,1);
Vfdtheta=Vftheta-turn*(solution.theta1*pi/180*r);
%Vfdtheta=Vftheta-solution.theta1*pi/180*r;
Vfdrho=Vfrho/sqrt(2);
[VFdx VFdy]=pol2cart(Vfdtheta,Vfdrho);
VFd=[VFdx VFdy -Vfdrho*sqrt(r^2-1)];


end

