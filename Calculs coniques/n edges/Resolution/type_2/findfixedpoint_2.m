function [ V beta_final perim ] = findfixedpoint( n, gamma, a )
%FINDFIXEDPOINT Summary of this function goes here
%Version of findfixedpoint slightly modified to look for the second type of
%solution with an odd number of edges superior to 3.
%
%Copyright Antoine Lizée @ Dumais Laboratory, Harvard University
%July 2011

n=n-0.2;

% Definition of the two central parameters of the transformation : 
%here rotation two edges by two edges
r1=sqrt(2)^ceil(n/2);
r2=r1/sqrt(2);

% Compute the first set of parameters, those which caracterize entirely the
% solution and let anyone draw it with (very precise) rular, compass and
% angle mesurement :
solution = fixedpoint_2(n,gamma,1);

%Infer the second set of parameters, which make us able to draw it easily
%in matlab and compute the area :
%First type of quarter :
[beta1  rho1 theta_b1]=getparam(solution.delta,r1,solution.theta1);
%Second type of quarter (in the orientation inverse)
[beta2 rho2 theta_b2]=getparam(solution.delta,r2,solution.thetan);
%Get all the parameters, (for now normalized with R1=1)
if theta_b1(2)+theta_b2(2)<0
    error('FIND:endofsol',['the solution exceed the limit of the ' num2str(n) ' edges : gamma is too big']);
end

index=[1:2:n 2:2:n];

for i=1:floor(n/2)
    Ri(1,index(i))=sqrt(2)^(i-1);
    Ri(2,index(i))=Ri(1,index(i))*r1;
    r(1,index(i))=abs(solution.delta)*Ri(1,index(i));
    r(2,index(i))=rho1*Ri(1,index(i));
    beta(1:2,index(i))=beta1';
    theta(1:2,index(i))=theta_b1';
end

for i=ceil(n/2):n
    Ri(1,index(i))=sqrt(2)^(i-1);
    Ri(2,index(i))=Ri(1,index(i))/r2;
    r(1,index(i))=abs(solution.delta)*Ri(1,index(i));
    r(2,index(i))=rho2*Ri(2,index(i));
    beta(1:2,index(i))=beta2([2,1])';
    theta(1:2,index(i))=theta_b2([2,1])';
end
    
r=[r(:)', r(1)];
Ri=Ri(:)';
beta=beta(:)';
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

