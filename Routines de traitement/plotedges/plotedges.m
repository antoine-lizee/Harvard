function [ graph ] = plotedges( N, varargin )
%% PLOTEDGES 
% This function is dedicated to the study of a cell network and its shape.
% It gives the possibility to visualize the way edges are connected
% together, by a plot of each "crossing" between three edges.
%   Detailed explanation goes here
% N is a Cell Network
% varargin is a paramater string, taking the option values :
%
% Concerning the calculation :
% - 'excl' if you want to plot only one time each edge, thereby plot only one
% quarter of the crossings
% - 'norm' if you want to normalize the right edge
% - 'normmean' if you want to normalize by the mean of the edge sizes
% - 'edgeswithangles' if you want to plot the real edge length and the real
% angles instead of the points (in case of angled edges)
% - 'low' prevents the algorithm from flipping the crossing patterns : it
% limits the sorting to the first condition, concerning the right edge, and
% do not use all the sorting conditions. Only useful for 'size' and 'angle'
% modes.
%
% Plot options :
% - 'rose' plots a polar histogram of the edge repartition.
% - 'hist' plots a classic histogram of the edge angle repartition.
% - 'box' plots a boxplot of the edge angle repartition
% - 'histsize' plots a classic histogram of the edge size repartition.
% - 'sub' plot in the same figure, thanks to the sublplots.
% - 'no' don't plot the edge representation
%
% 3d options :
% - '3d' enable the 3d hillcock visalization computing
% - '3d2' is another version of the latter, with variation of the
% hillcock width depending the edge length
%
% Different sorting strategies :
% - 'angle' if you want to sort by angles (default) : the edge with the smallest angle
% is on the top, the second smallest on the bottom of the graph (DEFAULT)
% - 'size' if you want to sort by size : the smallest edge is on the right,
% the second smallest on the top of the graph.
% - 'size_angle' if you want to put the smallest vertex on the right, then
% sort by angle the two other edges (the smallest angle will be on the top)
%
%See also : hist2, contour3d (script), Stack_... (scripts)
%
%Copyright Antoine Lizée 05/2011 @ Dumais Lab, OEB DPT, Harvard


%% Initialize the plot struct :
N_vertices=length(N.v);
graph.v=zeros(N_vertices, 4); % The 3 neighbor vertices followed by the central vertex
graph.x=zeros(N_vertices, 3); % The coordinates of the 3 neighbor vertices
graph.y=graph.x;
graph.a=graph.x; % The angle information of the three edges
graph.r=graph.x; % Polar coordinates of the edges extremities
graph.t=graph.x;
graph.d=graph.x; % Angles between the edges
graph.o=graph.x; % Order by angle
graph.s=graph.x; % Order by size
% graph.z is the array for 3D Histogram
% graph.stat is a struct with the statistical data
% graph.h regroup the axes_handles plotted

%% Unpacking the options
mode=1;
norm=0;
excl=0;
color=0;
angle=0;
low=0;
ros=0;
sub=0;
his=0;
box=0;
hissize=0;
edges=1;
hist3d=0;
n=0;
stat=0;
if length(varargin)>=1
    for j=1:length(varargin)
        switch varargin{j}
            case 'excl', excl=1; 
            case 'norm', norm=1;
            case 'normmean', norm=2;
            case 'color', color=1; 
            case 'edgeswithangles', angle=1;
            case 'low', low=1;
            case 'angle', mode=1;
            case 'size', mode=2;
            case 'size_angle', mode=3;
            case 'rose', ros=1;
            case 'hist', his=1;
            case 'box', box=1;
            case 'histsize', hissize=1;
            case '3d', hist3d=1;
            case '3d2', hist3d=2;
            case 'no', edges=0;
            case 'sub', sub=1;
            case 'stat', stat=1;
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
    graph.x(i,1:3)=N.v(neighbv,1)-[N.v(nv(i),1); N.v(nv(i),1); N.v(nv(i),1)] ;
    graph.y(i,1:3)=N.v(neighbv,2)-[N.v(nv(i),2); N.v(nv(i),2); N.v(nv(i),2)] ;
    graph.v(i,1:4)=[neighbv nv(i)];
    if angle % Edges with Angles option
        graph.a(i,1:3)=N.e(neighbe,3).*(-1).^(c==1); % Store the angle information with regards to the edge orientation (if c(k)=1, then the kth neighbor vertex was the first of the edge, and we must flip the angle of this edge)
    end
    if excl % Exclusion option
        nv(nv(:)==neighbv(1))=[]; % Suppress the neighbor vertices indexes from the nv vector.
        nv(nv(:)==neighbv(2))=[];
        nv(nv(:)==neighbv(3))=[];
    end                             
end

%% Reduce the size of the preallocated arrays
N_v=find(graph.x==0 & graph.y==0, 1)-1;
for i=['v', 'x', 'y', 'r', 't', 'a', 'o', 's', 'd']
    graph.(i)=[graph.(i)(1:N_v,:) []];
end
    

%% Convert it to polar coordinates :
[graph.t graph.r]=cart2pol(graph.x, graph.y);

%% Option : Compute the real length and angles :
if angle
    boola=find(graph.a);
    graph.r(boola)=graph.r(boola).*graph.a(boola)./sin(graph.a(boola));
    graph.t(boola)=graph.t(boola)-graph.a(boola); %The - comes from the choice of orientation of the program
end

%% Compute the angles between the edges :
[graph.t index]=sort(graph.t,2); % Put the edges in the trigonometric order (useful for later sorting)
for i=1:N_v % Sort the size information accordingly (Impossible without for loop ?)
    graph.r(i,:)=graph.r(i,index(i,:));
end
graph.d=graph.t-circshift(graph.t,[0,1]); % Ex : the angle between edge 1 and 2 is the second
graph.d=graph.d+2*pi.*(graph.d<0); % modulo 2pi

%% Give the "orders" : WARNING : it is not the order of the point, but the index of the point corresponding to each position
[~, graph.o]=sort(graph.d,2);
[~, graph.s]=sort(graph.r,2);

%% Sort accordingly of the mode : 
switch mode
    case 1 % 'angle'case
        bool=zeros(N_v,3);
        for i=1:N_v % Localise the edgethat will be on the right (i.e., on this case, in front of the biggest angle)
            bool(i,graph.o(i,3))=1;
        end
        bool=circshift(bool,[0,1]); % Array of 0 and one to localize the "base" edges (that will go on the right)
        minus=(graph.t.*bool)*ones(3,3); % Draw an appropriate array
        graph.t(:,1:3)=graph.t(:,1:3)-minus; % Put this edge on the right
        graph.t(graph.t<-pi)=graph.t(graph.t<-pi)+2*pi; % Recompute the angles to fit in [-Pi, Pi]
        graph.t(graph.t>pi)=graph.t(graph.t>pi)-2*pi;
        j=0;
        if ~low % if low is not enabled, we must flip the Vertex to put the biggest angle on the top
            for i=1:N_v
                if graph.t(i,graph.o(i,2))==0
                    graph.t(i,:)=-graph.t(i,:);
                    j=j+1;
                end
            end
        end
    case 2 % "size' case
        bool=zeros(N_v,3);  % bool is the array of 0 and one to localize the "base" edges (that will go on the right)
        for i=1:N_v % Localise the edgethat will be on the right (i.e., on this case, the bigger one)
            bool(i,graph.s(i,3))=1;
        end
        minus=(graph.t.*bool)*ones(3,3); % Draw an appropriate array
        graph.t(:,1:3)=graph.t(:,1:3)-minus; % Put this edge on the right
        graph.t(graph.t<-pi)=graph.t(graph.t<-pi)+2*pi; % Recompute the angles to fit in [-Pi, Pi]
        graph.t(graph.t>pi)=graph.t(graph.t>pi)-2*pi;
        j=0;
        if ~low % if low is not enabled, we must flip the Vertex to put the biggest angle on the top
            for i=1:N_v
                if graph.t(i,graph.s(i,2))==0
                    graph.t(i,:)=-graph.t(i,:);
                    j=j+1;
                end
            end
        end
end

display(sprintf('%d vertices have been flipped',j)),
            
            
            
            
%% Option : normalisation
switch norm
    case 1        
        div=(graph.r.*bool)*ones(3,3);        
        graph.r=graph.r./div;            
    case 2
        div=graph.r*ones(3)./3;
        graph.r=graph.r./div;
end

%% Option : color sorting
%given up

%% Plots

% Plot the classic representation
if edges
    if sub, h(n+1)=subplot(2,4,2); end;
    polar(graph.t', graph.r','.r');
    title('Plot of the extremeties of the edges');
    n=n+1;
end


% Plot rose
if ros
    if sub, h(n+1)=subplot(2,4,3); 
    elseif n==1;
    else figure(2),
    end;
    plot_rose=graph.t(:);
    plot_rose(plot_rose==0)=[]; % Suppress the zero values
    rose(plot_rose,200);
    title('Polar Histogram of the edges angles');
    n=n+1;
end

% Plot hist
if his
    if sub, h(n+1)=subplot(2,4,4); else figure(3),end;
    hist_graph=graph.t(:)*180/pi;
    hist_graph(hist_graph==0)=[];
    [Na xouta]=hist2(hist_graph,100,1);
    Na=Na*100/N_v; % Convert in "%" for each edge (bottom and top)
    plot(xouta, Na,'-');
    title('Histogram of the angles of the up and down edges');
    n=n+1;
end
if hissize
    if sub, h(n+1)=subplot(2,4,5); 
    elseif n==0;
    else figure(4),
    end;
    [graph.t index]=sort(graph.t,2); % Sort the edges by their position on the vertex;
    for i=1:N_v 
        graph.r(i,:)=graph.r(i,index(i,:));
    end
    [Ns xouts]=hist2(graph.r,60,5);
    Ns=Ns*100/N_v; % Convert in "%" for each edge (bottom, top and right)
    %Plot of the histogram plus the lines :
    %bar(xouts, Ns, 'BarWidth', 0.5,'LineStyle', 'none');
    %hold on;
    %colors=[1 0.2 0.2;0 1 0;0.2 0.2 1];
    %colormap(colors); % For the bar plot
    %set(gca,'ColorOrder', colors); % For the plot
    %dec=(xouts(2)-xouts(1))/5; % Plot on the bar plot -> we must move the curves
    %xouts=[xouts-dec xouts xouts+dec];
    plot(xouts, Ns, '-');
    legend('bottom edge', 'right edge', 'top edge');
    title('Histogram of the edges size');
    n=n+1;
end


%% Statistical Data computation

% Sort the colums with the position of the edges (if not already done in
% histsize)
if stat
    if ~hissize
        [graph.t index]=sort(graph.t,2); % Sort the edges by their position on the vertex;
        for i=1:N_v 
            graph.r(i,:)=graph.r(i,index(i,:));
        end
    end

    % Compute the statistical values
    stats=zeros(2,3); % Array of the stat data concerning the size of the edges
    stata=zeros(2,3); % Array of the stat data concerning the angle of the edges
    stats(1,:)=mean(graph.r);
    stats(2,:)=std(graph.r);
    stata(1,:)=mean([-graph.t(:,1) 2*pi+graph.t(:,1)-graph.t(:,3) graph.t(:,3)]);
    stata(2,:)=std([graph.t(:,1) 2*pi+graph.t(:,1)-graph.t(:,3) graph.t(:,3)]);
    stata=stata./pi*180;

    % Return the values
    graph.stat.a=stata;
    graph.stat.s=stats;

    % Draw the stats in a plot
    if sub, h(n+1)=subplot(2,4,6); else figure(11),end;
    axis off;
    txa{1}='\bfConcerning the angles :\rm';
    txa{2}='             bottom angle |  left angle  |  top angle  ';
    txa{3}=['means :     ' num2str(stata(1,:), '%15.1f')];
    txa{4}=['std dev. :   ' num2str(stata(2,:), '%15.1f')];

    txs{1}='\bfConcerning the sizes :\rm';
    txs{2}='             bottom edge |  right edge  |  top edge  ';
    txs{3}=['means :     ' num2str(stats(1,:), '%15.2f')];
    txs{4}=['std dev. :   ' num2str(stats(2,:), '%15.2f')];

    text(0.5,0.75,txa, 'EdgeColor', 'k',...
        'BackgroundColor', 'w','HorizontalAlignment', 'Center',...
        'VerticalAlignment', 'Middle', 'Margin', 5);
    text(0.5,0.25,txs, 'EdgeColor', 'k',...
        'BackgroundColor', 'w','HorizontalAlignment', 'Center',...
        'VerticalAlignment', 'Middle', 'Margin', 5);
    n=n+1;

    % Draw boxplots
    if box
        if sub, h(n+1)=subplot(2,4,7); else figure(12),end;
        boxplot([-graph.t(:,1) 2*pi+graph.t(:,1)-graph.t(:,3) graph.t(:,3)]*180/pi,...
            'label', {'bottom angle','left angle','top angle'});
        title('Boxplot of the distribution of the angles');
        n=n+1;
    end
end
    
%% Managing the subplots
if sub
    figure(999); % Initialize the figure for resizing plotting
    switch n
        case {1,2}
            for i=1:n
                h2(i)=subplot(1,n,i);
            end
        case 3
            h2=[subplot(2,2,[1,3]); subplot(2,2,2); subplot(2,2,4)];
        case 4
            for i=1:n
                h2(i)=subplot(2,2,i);
            end
    end
    drawnow;
    for i=1:n
        set(h(i),'Position',get(h2(i),'Position'))
    end

    delete(999);
    figure(gcf);
    graph.h=h;
end  

%% 3D Histogram Plot
if hist3d
    res=500; %Resolution of the plot
    
     % Preparing for computation
    [graph.x graph.y]=pol2cart(graph.t, graph.r);
    suppr=graph.t(:)==0; % Suppress the edges on the right
    for i=['x', 'y', 'r', 't']
        h3d.(i)=graph.(i)(:);
        h3d.(i)(suppr)=[];    
    end
    N_h3d=length(h3d.r);
    
    % Computing the appropriate width of each "hillock":
    if norm~=0 % Compute the longitudinal total width
        l=1;
    elseif norm==0
        l=sqrt(mean(mean(h3d.r)));
    end
    a=50; 
    w=a*l/1000;
    
   
    % Compute the coordinates to be associated with each point :
    y=linspace(max(h3d.y)+4*w,min(h3d.y)-4*w,res);
    x=linspace(min(h3d.x)-4*w,max(h3d.x)+4*w,res);
    [X Y]=meshgrid(x,y);

    % Compute the map of values :
    bar = waitbar(0,['Preparing for computation of ' num2str(length(N_h3d)) ' hillocks']);
    r2=zeros(res,res);
    z=zeros(res, res);
    totz=zeros(res, res);
    c1=sqrt(2*pi)*w*ones(N_h3d,1);
    c2=2*w^2*ones(N_h3d,1);
    cnorm=1000/N_h3d;
    if hist3d==2, % variation of the hillock width with the length of the edge
        c1=sqrt(2*pi)*w*h3d.r;
        c2=2*(w*h3d.r).^2;
    end 
    for k=1:N_h3d
        waitbar(k/N_h3d,bar,['Computing the hillock number ' num2str(k) ' out of ' num2str(N_h3d)]);
        r2(:,:)=(Y-h3d.y(k)).^2+(X-h3d.x(k)).^2;
        z(:,:)=1./c1(k)*exp(-r2(:,:)./c2(k));
        totz=z+totz;
    end
    totz=totz*cnorm;
    deb=sum(sum(z(:,:)));
    close(bar);
    % totz=sum(z,3); % Old version, with z in 3d, faster but occupying to
    % much memory...
    
    %Plot the surface
    if ~isempty(get(gcf,'CurrentAxes')),
        figure(50); 
    end
    surf(x,y,totz,'EdgeColor','none');
    colormap hsv;
    daspect([1 1 50]);
    view(2);
    n=n+1;
        
    % Returning the information
    graph.z.x=x;
    graph.z.y=y;
    graph.z.totz=totz;
    graph.z.N_h3d=N_h3d;
    graph.z.deb=deb;
end

%% Exceptions

if ~isempty(not_plotted)
    display(sprintf('WARNING : %d vertices have not been plotted because they had more or less than 3 edges connected to them \n', length(not_plotted))),
    % sprintf('Vertex of index %d has %d neighbors\n', not_plotted(:,1), not_plotted(:,2)), 
end

end

