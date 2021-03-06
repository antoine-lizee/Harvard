%% boundary
% Find the intersection of a half line with a cell boundary.
%%

%% Syntax   
% Q = boundary(N,theta);
%
%% Description 
% This program computes the intersection of the line starting at the
% centroid of the cell and forming an angle theta with the horizontal line
% with the boundary of the cell. 
%
%% Inputs
% * N - a celNetwork object
% * theta - the angle formed by the line with the horizontal axis in the
% range [0, 2*pi]
% 
%% Outputs
% * B - the boundary network
%
%% Examples
% >> findQ([1 0; 0 1;0 0],[pi/4  0 0],pi,.5)
%
% >> findQ([1 0; 0 1;0 0],[pi/4 0 0],pi/4)
%
%% TODO
%
%% See also
% * findP
% * cellcentre
%
%% Author    
% Sebastien Besson
% email address: sbesson@oeb.harvard.edu
% November 2007; Last revision: November 8 2007

function Q = findQ(V,A,centre,theta,s)

% Check number of inputs
error(nargchk(4, 5, nargin));

% Check the dimensions of the inputs
isValid = (size(V,2) == 2 & size(A) == size(V,1));
if ~isValid, error('Dimensions of the inputs do not agree with the definition of the function.'); end

% Set the optional argument
if (nargin == 4), s = 1; end

% Computes the centre of the cell
%centre = cellcentre(V,A);

% Length of the cordlength of the edge
T = circshift(V, [-1 0]) - V;
N= [T(:,2) -T(:,1)];

V = [V; V(1,:)];

% Loop over the edges o the cell
for i = 1:size(V)-1
    if A(i)^2 <= 0.00001
        % Look for the intersection between two lines of equations :
        % [x,y] = centre + u*[cos(theta) sint(theta)]
        % (y-V(i,2)) = a*(x-V(i+1,1)) +b
        if (V(i,1) == V(i+1,1))
             % Equation of the type x= xi
             u = (V(i,1) - centre(1))/cos(theta);
             % Test if the solution belongs to the segment [V(i),V(i+1)]
             t = (centre(2) + u*sin(theta) - V(i,2))/(V(i+1,2) - V(i,2));
             if ((t>=0) && (t<=1) && (u>0)), break; end;
        else
            % Equation of the cordlength of the type y-V(i,2) = a*(x-V(i,1))
            a = (V(i+1,2) - V(i,2))/(V(i+1,1) - V(i,1));
            % Test if the lines are parallel
            if (sin(theta) == a*cos(theta)), continue; end;
            u = (a*(centre(1)-V(i,1)) - (centre(2) - V(i,2)))/(sin(theta) - a*cos(theta));

            % Test if the solution belongs to the segment [V(i),V(i+1)]
            t = (centre(1) + u*cos(theta) - V(i,1))/(V(i+1,1) - V(i,1));       
            if ((t>=0) && (t<=1) && (u>0)), break; end;
        end
    else
        % Look for the intersection between  a line and a circle of equations :
        % [x,y] = centre + u*[cos(theta) sint(theta)]
        % (x-C(1))^2+(y-C(2))^2 = R^2
        R = norm(T(i,:))/(2*sin(A(i)));
        C = V(i,:) + T(i,:)/2 -N(i,:)/(2*tan(A(i)));

        % 2nd order equation to solve 
        % u^2 +u*(2*(centre(1)-C(1))*cos(theta)+ 2*(centre(2)-C(2))*sin(theta))
        % + (centre(1)-C(1))^2+ (centre(2)-C(2))^2 -R^2 = 0
        b = 2*(centre(1)-C(1))*cos(theta) + 2*(centre(2)-C(2))*sin(theta);
        c = (centre(1)-C(1))^2 + (centre(2)-C(2))^2 - R^2;
        delta = b^2 - 4*c;

        % Test the existence of solutions
        if (delta<0), continue; end;

        % 1st solution
        u = -b/2 - sqrt(delta)/2;
        Ql = centre + u*[cos(theta) sin(theta)];  
        if (u>0) 
            alpha0 = atan2(V(i,2)-C(2),V(i,1)-C(1));
            alpha = atan2(Ql(2)-C(2),Ql(1)-C(1));
            alpha1 = atan2(V(i+1,2)-C(2),V(i+1,1)-C(1));
            t =  (alpha-alpha0)/(alpha1-alpha0);
            if ((t>=0) && (t<=1)), break; end;
        end

        % 2nd solution
        u = -b/2 + sqrt(delta)/2;
        Ql = centre + u*[cos(theta) sin(theta)];  
        if (u>0) 
            alpha0 = atan2(V(i,2)-C(2),V(i,1)-C(1));
            alpha = atan2(Ql(2)-C(2),Ql(1)-C(1));
            alpha1 = atan2(V(i+1,2)-C(2),V(i+1,1)-C(1));
            t =  (alpha-alpha0)/(alpha1-alpha0);
            if ((t>=0) && (t<=1)), break; end;
        end
    end
    
    % Error message if no intersection point is found
    if (i == size(V)-1), error('No intersection point found'); end;
end

% Find the intersection point
Ql = centre + u*[cos(theta) sin(theta)];

% Find the intermediate point
Q = centre + s*(Ql-centre);

end