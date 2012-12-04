function [ s1 s2]=finds2( s1, M, i, j )
%FINDS2 is a lightweight version of finda that computes the abscissa s2 on the
%ending edge of the dividing wall starting on s1 on the starting edge.
%
% Remember that we always turn counter clockwise. This is crucial for all
% geometrics calculous
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu


%% Exceptions (neighbour edges)
if (s1==1) && isequal(circshift(M.I_e,[j,0]),circshift(M.I_e,[i+1,0]))
    s2=0;
    return;
elseif (s1==0) && isequal(circshift(M.I_e,[j,0]),circshift(M.I_e,[i-1,0]))
    s2=1;
    return;
end

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

%% Find s4 such that the arc s3-s4 cross the two edges with perpendicular angles
%Test
%% separation of the four cases
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
d2=d2+(d2<=-3*pi)*2*pi+(d2<=-pi)*2*pi-(d2>pi)*2*pi-(d2>3*pi)*2*pi;


s2=d2/(2*M.a(j));

VR2=M.C(j,:)+p2c(A2,1)*(r2+or(2)*r);

% Debug : Uncomment to see points of interest, after having plotted
% the tissue under treatment.
%     plot([V1(1) VR(1)],[V1(2) VR(2)],'-og');
%     plot([V2(1) VR2(1)],[V2(2) VR2(2)],'-*b');
%     plot([VR2(1)],[VR2(2)],'ob');
%     plot([M.C(j,1) V2(1)],[M.C(j,2) V2(2)],'-xm');

    

deltaR=VR2-VR;

if norm(deltaR)>10^(-9) ;  % Then we do the "tangent computation"
    disp(['WARNING : For cell number ' num2str(M.cell,' %.0f') ' of the tissue '' ' M.tissue ''''])
    disp(['Tangent computation is enabled for the computation of the division wall from edge '...
        num2str(i,' %.0f') ' to edge ' num2str(j,' %.0f') ' at abscissa ' num2str(s1) ]);
    disp(['deltaR is equal to' num2str(deltaR)]);
    x0 = [0.5; 1];           % Make a starting guess at the solution
    options=optimset('Display','off');   % Option to display output
    [x,fval,exitflag,output] = fsolve(@myfun,x0,options);  % Call solver
    if exitflag~=1, display(output.message), end
    s2= x(1);
 
end  

%%Nested functions
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
