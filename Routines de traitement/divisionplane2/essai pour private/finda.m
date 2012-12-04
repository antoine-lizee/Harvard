function [ s1 s2 a V A]=finda( s1, M, i, j )
%FINDA is the function that compute the geometrical core and give for any
%abscissa s1 of the start edge the full characterization of the dividing wall
%which start from that abscissa to another edge and the area ratio obtained
%by this dividing wall.
%
% Remember that we always turn counter clockwise. This is crucial for all
% geometrics calculous
%
%%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.

%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu



%% Initialisation
r1=M.R_e(i);
r2=M.R_e(j);
V1=M.C(i,:)+p2c(M.alpha(i,1)+2*s1*M.a(i),r1);
T1=p2c(M.alpha(i,1)+2*s1*M.a(i)+M.s(i)*pi/2,1);
R1=p2c(M.alpha(i,1)+2*s1*M.a(i),1);

%% First computations
C1C2=M.C(j,:)-M.C(i,:);
L=C1C2*R1';
l=C1C2*[-R1(2);R1(1)];

%% Test for curvature and run computation for straight edges cases


%% Find s4 such that the arc s3-s4 cross the two edges with perpendicular angles
%Test

% separation of the four cases
if M.a(i)>0 && M.a(j)>0
    or=[-1,-1];
elseif M.a(i)>0 && M.a(j)<0
    or=[-1,1];
elseif M.a(i)<0 && M.a(j)>0
    or=[1,-1];
elseif M.a(i)<0 && M.a(j)<0
    or=[1,1];
else
    error('Finda:edgesorting','The treatment of edge with curvature must not be run if at least one of the edge is straight')
end


     % We first do the "radius computation"
    %We do the heart of the computation, the calculous of 'r', which will
    %give us the point R, the latter being the intersection between the two
    %radius of the two circles 1 and 2. You will have to draw a figure to get even
    %a small part of the (quite simple) reasoning here.
    r=(l^2-r2^2+(L-r1)^2)/(2*(or(1)*(L-r1)+or(2)*r2)); % This arise from simple pythagore-like calculation
    
    VR=M.C(i,:)+R1*(r1+or(1)*r); %These are the coordinates of the point R, given by the calculation of r.
    CR2=VR-M.C(j,:); % This is the vector between the Center of the target edge, C2, and R
    A2=atan2(CR2(2),CR2(1))+(or(2)==-1 && r>r2)*pi +(or(2)==1 && r<-r2)*pi;
    d2=A2-M.alpha(j,1);
    
%     % We want d2 to be in an interval such that s2, and therefore the area
%     % calculation make no "jump" between to solutions
%     C2vi2=M.v2(i,:)-M.C(j,:); %This is the vector between the Center of the target edge, C2, and the "last" (and second) vertex of the first edge.
%     if or(2)<0
%         asup=atan2(C2vi2(2), C2vi2(1))-M.alpha(j,1);
%         if asup
%         asup=asup+(asup<=0)*2*pi+(asup<=-2*pi)*2*pi; % From [-3pi pi] to ]0 2pi]
%     else
%         asup=pi;
%     end
%     ainf=asup-2*pi;    
    %disp([asup; ainf]);
    %disp(d2);
%     d2=d2+(d2<=(ainf-2*pi))*2*pi+(d2<=ainf)*2*pi-(d2>asup)*2*pi-(d2>(asup+2*pi))*2*pi; %d2=d2 modulo2pi between ainf and asup;
    d2=d2+(d2<-3*pi)*2*pi+(d2<-pi)*2*pi-(d2>pi)*2*pi-(d2>3*pi)*2*pi;
    %disp(d2);
    
    
    s2=d2/(2*M.a(j));
    V2=M.C(j,:)+p2c(A2,r2);
    V=[V2;V1];
    Exy=V1-V2;
    angle=atan2(Exy(2),Exy(1))-A2+pi*(or(2)==-1);
    angle=angle+(angle<-pi)*2*pi-(angle>pi)*2*pi;
    
    A=[M.a(i)*(1-s1); M.a(j)*s2; angle]; 
    
    VR2=M.C(j,:)+p2c(A2,1)*(r2+or(2)*r);
    
    % Debug : Uncomment to see points of interest, after having plotted
    % the tissue under treatment.
%     plot([V1(1) VR(1)],[V1(2) VR(2)],'-og');
%     plot([V2(1) VR2(1)],[V2(2) VR2(2)],'-*b');
%     plot([VR2(1)],[VR2(2)],'ob');
%     plot([M.C(j,1) V2(1)],[M.C(j,2) V2(2)],'-xm');

    

deltaR=VR2-VR; %Detect problem of convergence with the first method

if norm(deltaR)>0.0000001 ;  % Then we do the "tangent computation"
    disp(['WARNING : For cell number ' num2str(M.cell,' %.0f') ' of the tissue '' ' M.tissue ''''])
    disp(['Tangent computation is enabled for the computation of the division wall from edge '...
        num2str(i,' %.0f') ' to edge ' num2str(j,' %.0f') ' at abscissa ' num2str(s1) ]);
    disp(['deltaR is equal to' num2str(deltaR)]);
    x0 = [0.5; 1];           % Make a starting guess at the solution
    options=optimset('Display','off');   % Option to display output
    [x,fval,exitflag,output] = fsolve(@myfun,x0,options);  % Call solver
    if exitflag~=1, display(output.message); end
    s2= x(1);
    m=x(2);
    
    V2=M.C(j,:)+p2c(M.alpha(j,1)+2*s2*M.a(j),r2);
    T2=p2c(M.alpha(j,1)+2*s2*M.a(j)+M.s(j)*pi/2,1); %Tangent vector of the edge going in the way of the curvilign abscissa (counter clockwise of the cell)
    
    

    V=[V2;V1];
    Exy=V1-V2;
    
    angle=asin(norm(Exy)/(2*m)); %*sign(det([M.T(i,:)',Exy'])));
    
    A=[M.a(i)*(1-s1); M.a(j)*s2; angle]; 
    
    % Debug : Uncomment to see points of interest, after having plotted
    % the tissue under treatment.
    plot([V1(1) V1(1)+m*T1(1)],[V1(2) V1(2)+m*T1(2)],'-og');
    %VR2=M.C(j,:)+p2c(A2,1)*(r2+or(2)*r);
    plot([V2(1) V2(1)-m*T2(1)],[V2(2) V2(2)-m*T2(2)],'-*b');
    %plot([VR2(1)],[VR2(2)],'ob');
    %plot([M.C(j,1) V2(1)],[M.C(j,2) V2(2)],'-xm');

    disp(fval)    
end  
%     % Debug 2 : Uncomment to see the boundaries of the new cells that are
%     % used for area calculation.
%     P=findP([V(2,:);M.v1(i+1,:)],A(1),0:0.01:1);
%     plot(P(:,1),P(:,2),'-k','Linewidth',2);
%     P=findP([M.v1(j,:);V(1,:)],A(2),0:0.01:1);
%     plot(P(:,1),P(:,2),'-k','Linewidth',2);

%% Compute area
a3=area([V1;M.v1(i+1:j,:);V2],[A(1);M.a(i+1:j-1);A(2:3)]);
a=a3/M.area;

%% Nested functions
function F=myfun(x)
    F=V1+T1*x(2)- ...
        (M.C(j,:)+p2c(M.alpha(j,1)+2*x(1)*M.a(j),r2)-...
        x(2)*p2c(M.alpha(j,1)+2*x(1)*M.a(j)+M.s(j)*pi/2,1));
end


end

function x=p2c(THETA,RHO)
% We need this function to skip the inherent nargout check of the classical
% plot2cart function.
[a b]=pol2cart(THETA,RHO);
x=[a b];
end
