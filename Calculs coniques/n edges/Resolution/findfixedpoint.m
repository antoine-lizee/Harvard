function [ V beta_final perim ] = findfixedpoint( n, gamma, a )
%FINDFIXEDPOINT
%This function computes the set of parameters needed to caracterize the
%fixed point solution of a given area, and gives also the perimeter of the
%solution.
% 'fixedpoint' gives us the solution but in a set of parameters (hi, Ri,
% thetai) that are not directly usable for matlab in order to trace and
% transform the solution to the cone. 'Findfixedpoint' is the global function
% that uses 'fixedpoint' then 'getparam' in order to deliver the solution in
% the proper set of parameters.
%
%Copyright Antoine Lizée @ Dumais Laboratory, Harvard University
%July 2011


% Management of 'Type 2' cases:
if n-0.2==floor(n)
    [ V beta_final perim ] = findfixedpoint_2( n, gamma, a );
    return
end


% Definition of the two central parameters of the transformation : 
%here rotation one edge by one edge
r1=sqrt(2);
r2=r1^(n-1);

% Compute the first set of parameters, those which caracterize entirely the
% solution and let anyone draw it with (very precise) rular, compass and
% angle mesurement :
solution = fixedpoint(n,gamma,1);

%Infer the second set of parameters, which make us able to draw it easily
%in matlab and compute the area :
%First type of quarter :
[beta1  rho1 theta_b1]=getparam(solution.delta,r1,solution.theta1);
%Second type of quarter (in the orientation inverse)
[betan rhon theta_bn]=getparam(solution.delta,r2,solution.thetan);
%Get all the parameters, (for now normalized with R1=1)
if theta_b1(2)+theta_bn(2)<0
    error('FIND:endofsol',['the solution exceed the limit of the ' num2str(n) ' edges : gamma is too big']);
end
rho=[rho1*r1.^(0:n-2) rhon];
r=[ abs(solution.delta)*r1.^(0:n-1) ; rho]; %hi and rho_i for i=1:n
r=[r(:)', r(1)];
R=r1.^(0:n-1);
Ri=[R;circshift(R,[0,-1])];
Ri=Ri(:)';
beta=[beta1'*ones(1,n-1) [betan(2);betan(1)] ];
beta=beta(:)';
theta=[theta_b1'*ones(1,n-1) [theta_bn(2);theta_bn(1)] ];
theta=theta(:)';

% Compute area :
area=getareageom(r,theta,Ri,beta);

%Normalize to the area and compute the values that characterize the cell.
R_final=1/sqrt(area/a);
theta_f=cumsum([gamma/2*pi/180 theta]);
theta_final=theta_f([1,2:2:2*n,2*n+1])';
rho_final=r([1,2:2:2*n,2*n+1]);
V=[theta_final, rho_final'*R_final]; %Vertices of the cell in polar coordinates
beta_f=beta(3:2:2*n-1)+beta(2:2:2*n-2);
beta_final=[beta(1),beta_f,beta(end)]; %angle between each edge and the cordlength

%Calculation of the perimeter
perim=sum(abs(2*beta_final.*(R_final.*[Ri(1:2:end) 1])));

end

