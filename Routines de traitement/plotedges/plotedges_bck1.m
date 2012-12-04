function [ plot ] = plotedges( N, varargin )
%% PLOTEDGES Summary of this function goes here
% This function is dedicated to the study of a cell network and its shape.
% It gives the possibility to visualize the way edges are connected
% together, by a plot of each "crossing" between three edges.
%   Detailed explanation goes here
% N is a Cell Network
% varargin is a paramater string, taking the values :
% - 'excl' if you want to plot only one time each edge, thereby plot only one
% quarter of the crossings
% - 'norm' if you want to normalize the right edge
% - 'color' if you want to sort the edges by their size and make it
% - 'edgeswithangles' if you want to plot the real edge length and the real
% angles instead of the points (in case of angled edges)
% appear with a color code.
% - 'low' prevents the algorithm from flipping the crossing patterns : it
% limits the sorting to the first condition, concerning the right edge, and
% do not use all the sorting conditions. Only useful for 'size' and 'angle'
% modes.
% And one of these different ploting strategies :
% - 'angle' if you want to sort by angles (default) : the edge with the smallest angle
% is on the top, the second smallest on the bottom of the graph
% - 'size' if you want to sort by size : the smallest edge is on the right,
% the second smallest on the top of the graph.
% - 'size_angle' if you want to put the smallest vertex on the right, then
% sort by angle the two other edges (the smallest angle will be on the top)


%% Initialize the plot struct :
N_vertices=length(N.v);
plot.v=zeros(N_vertices, 4); % The 3 neighbor vertices followed by the central vertex
plot.x=zeros(N_vertices, 3); % The coordinates of the 3 neighbor vertices
plot.y=plot.x;
plot.a=plot.x; % The angle information of the three edges
plot.r=plot.x; % Polar coordinates of the edges extremities
plot.t=plot.x;
plot.d=plot.x; % Angles between the edges
plot.o=plot.x; % Order by angle
plot.s=plot.x; % Order by size

%% Unpacking the options
mode=1;
norm=0;
excl=0;
color=0;
angle=0;
low=0;
if length(varargin)>=1
    for j=1:length(varargin)
        switch varargin{j}
            case 'excl', excl=1; 
            case 'norm', norm=1; 
            case 'color', color=1; 
            case 'edgeswithangles', angle=1;
            case 'low', low=1;
            case 'angle', mode=1;
            case 'size', mode=2;
            case 'size_angle', mode=3;
            otherwise, error('wrong option name for option input number %d',j);
        end
    end
end


not_plotted=[];

%% Create the arrays with the neighbour vertices
nv=1:N_vertices;
i=0;
while i<length(nv)
    i=i+1;
    [neighbe c]=find(N.e(:,1:2)==nv(i)); % find the position of the vertex index i in the edge array (neighbe is therefore the three neighbor edge indexes)
    if length(neighbe)~=3
        not_plotted=[not_plotted; nv(i), length(neighbe)]; % Store the exception
        nv(i)=[]; % Suppress the faulty point
        i=i-1; % re-indent i
        continue; % End the loop
    end 
    c=1+(c==1); %invert the columns
    neighbv=[N.e(neighbe(1),c(1)) N.e(neighbe(2),c(2)) N.e(neighbe(3),c(3))] ; % compute the vertices index of the neighbor vertices
    plot.x(i,1:3)=N.v(neighbv,1)-[N.v(nv(i),1); N.v(nv(i),1); N.v(nv(i),1)] ;
    plot.y(i,1:3)=N.v(neighbv,2)-[N.v(nv(i),2); N.v(nv(i),2); N.v(nv(i),2)] ;
    plot.v(i,1:4)=[neighbv nv(i)];
    if angle % Edges with Angles option
        plot.a(i)=N.e(neighbe,3)'.*(-1).^(c==1); % Store the angle information with regards to the edge orientation (if c(k)=1, then the kth neighbor vertex was the first of the edge, and we must flip the angle of this edge)
    end
    if excl % Exclusion option
        nv(neighbv)=[]; % Suppress the neighbor vertices indexes from the nv vector.
    end                             
end

%% Reduce the size of the preallocated arrays
N_v=find(plot.x==0 & plot.y==0, 1)-1;
for i=['v', 'x', 'y', 'r', 't', 'a', 'o', 's', 'd']
    plot.(i)=[plot.(i)(1:N_v,:) []];
end
    

%% Convert it to polar coordinates :
[plot.t plot.r]=cart2pol(plot.x, plot.y);

%% Option : Compute the real length and angles :
if angle
    plot.r=plot.r.*plot.a./sin(plot.a);
    plot.t=plot.t+plot.a;
end

%% Compute the angles between the edges :
[plot.t index]=sort(plot.t,2); % Put the edges in the trigonometric order (useful for later sorting)
for i=1:N_v % Sort the size information accordingly (Impossible without for loop ?)
    plot.r(i,:)=plot.r(i,index(i,:));
end
plot.d=plot.t-circshift(plot.t,[0,1]); % Ex : the angle between edge 1 and 2 is the second
plot.d=plot.d+2*pi.*(plot.d<0); % modulo 2pi

%% Give the "orders" : WARNING : it is not the order of the point, but the index of the point corresponding to each position
[~, plot.o]=sort(plot.d,2);
[~, plot.s]=sort(plot.r,2);

%% Sort accordingly of the mode : 
switch mode
    case 1 % 'angle'case
        % That is a faulty version, considering that plot.o was the order
        % of the points :
        % bool=circshift(plot.o==3,[0,1]); % Localise the edgethat will be on the right (i.e., on this case, in front of the biggest angle)
        % minus=(plot.t.*bool)*ones(3,3); % Draw an appropriate array
        minus = zeros(N_v,3);
        for i=1:N_v
            minus(i,:)=plot.t(i,plot.o(i,2))*[1 1 1]; % The biggest angle is between the edges 1 and 3 of plot.o, so we take edge 2
        end
        plot.t(:,1:3)=plot.t(:,1:3)-minus; % Put this edge on the right
        if ~low % if low is not enabled, we must sort the two other edges accordingly
            i_switch=plot.o(bool)==1 ; % the indexes of the crossings we must flip
            plot.t(i_switch)=2*pi*eyes(length(i_switch),3)-plot.t(i_switch); % Flipping of the crossing
        end
    case 2 % "size' case
end
            
            
            
            
            
%% Option : normalisation
if norm
    div=(plot.r.*bool)*ones(3,3);
    plot.r=plot.r./div
end
            
%% Plot
axis on
hold on
polar(plot.t, plot.r,'.r');


%% Exceptions

if ~isempty(not_plotted)
    sprintf('WARNING : %d vertices have not been plotted because they had more or less than 3 edges connected to them', length(not_plotted)),
    % sprintf('Vertex of index %d has %d neighbors\n', not_plotted(:,1), not_plotted(:,2)), 
end

end

