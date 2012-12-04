function [ sa sb V A af ] = findplaneii( M, i, j, a )
%FINDPLANEII is a version of FINDPLANE which tackle the search of dividing
%plane in the case of the same edge for starting and ending of this wall.
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu

if i~=j
    error('findplaneii must be called only for division on the same edge (i=j)')
end

if nargin<4, a=0.5; end

s1=findboundsii(M,i,j);
if isempty(s1), 
    A=[];
    V=[];
    af=[];
    sa=[];
    sb=[];
    return;
elseif s1==[0, 0],
    A=[];
    V=[];
    af=[];
    sa=1;
    sb=[];
    return;
end
[sa, sb, af, V, A]=findsii(s1(1), s1(2), a, M, i, j, 1000/range(s1));

% Check for intersections with other edges
for k=[1:i-1, i+1:length(M)]
    if abs(M.a(i))>0.00001
        % Calculate the properties of the new edge needed to do the check
        % (see findcarac) :
        new.xy=V(2,:)-V(1,:); % Coordinates of the edge vector
        new.l_e=sqrt(sum(new.xy.^2,2)); % Length of the edge
        new.T=new.xy./(new.l_e*[1,1]); % Tangent vector of the edge
        new.N=[-new.T(:,2) new.T(:,1)]; % normal vector, "left"-rotated (positively)-rotated
        new.r_e=new.l_e./(2*tan(new.a)); % Distance from the edge to the center of the circle
        new.R_e=new.l_e./(2*sin(new.a)); % curvature radius. Can take the "Inf" value
        new.VC=new.T.*(new.l_e/2*[1,1]) + new.N.*(new.r_e*[1,1]); % coordinates of the distance from the edge to the center
        new.C=V(1,:) + M.VC; % Coordinates of the center
        
        
        diff=M.C(i,:)-M
        test=M.C(i,:)

end

