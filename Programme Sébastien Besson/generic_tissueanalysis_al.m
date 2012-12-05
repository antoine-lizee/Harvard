%% Routine for treating an extracted tissue - revised version
% Assume you have a fully defined tissue

%clear
clc
close all
warning on

% Name of the folder & file
[file folder]=uigetfile;

mergetissue = 1;
simplify =1;
analyse =0;
compare =0;

%% Remove the youngest edges

if mergetissue % Remove youngest edges 
    load([folder file]); % Load the file
    
    %sauvegarde_bis; % To check the relevance of the new version of the program - Debugging fonctionality
        
    % Case of a treated tissue :
    % tissue=originaltissue;
        
    originaltissue = tissue; % Useless for treated cases
    check(tissue);
    tissue.e(:,3) = 0; % Remove the angle information

    Io_edgestomerge = find(tissue.e(:,4)==2); % Create the vector with the indexes of the edge to suppress (in originaltissue)
    N_edgestomerge = length(Io_edgestomerge); % Determine number of edges to suppress
    Iom_verticestosuppress(:,1:2)= originaltissue.e(Io_edgestomerge,1:2); % Get the list of couple of vertices index for the edge we are supressing
    for i=1:N_edgestomerge
        next_edge = find(tissue.e(:,4)==2,1); % Find the indice of the next recent edge
        tissue = merge(tissue,neighbors(tissue,next_edge)); % Merge the two neighbor cells
    end
    mergedtissue =tissue; % Save the merged tissue
    save([file(1:end-4) '_bis.mat'],'originaltissue','mergedtissue','Iom_verticestosuppress');
end


%% Remove the vertices and create a new edge from the two bounded

if simplify
    %load([file(1:end-4) '-treated.mat']);
    
    figure('Position',[20 20 1000 500])    
    axes('OuterPosition',[0 0 1/2 1]);
    h0 = plot(originaltissue); axis equal;   
    axes('OuterPosition',[1/2 0 1/2 1]);
    h = plot(mergedtissue); axis equal;
    
    tissue = mergedtissue; % Retrieve the merged version of the tissue
    Is_newedges = zeros(size(Iom_verticestosuppress)); % Initialise the edge index list (in simplifiedtissue)
    
    % Sort Iom_vertices in a vector and suppress vertices which are employed on two
    % young edges (exception of the three edge per vertex law) :    
    nv=unique(Iom_verticestosuppress)';
    
    bar = waitbar(0,['Removing' num2str(length(nv)) 'vertices']);
    
    for i=nv(end:-1:1) % i indice of the vertex to be removed, initially in mergedtissue, but also in tissue thanks to the backward looping
        waitbar(1-i/length(nv),bar,['Removing vertex: ' num2str(i)]);
        
        % Find the neighbors of the vertex
        [neighbe c]=find(tissue.e(:,1:2)==i); %find the position of the vertex index i in the edge array (neighbe is therefore the two neighbor edge indexes)
        if length(neighbe)~=2, error('More or less than 2 edges neighbor to the vertex to suppress'); end
        c=1+(c==1); %invert the columns
        neighbv=[tissue.e(neighbe(1),c(1)) tissue.e(neighbe(2),c(2))] ; % compute the vertices index of the neighbor vertices
        
        % Compute an interlligent version of "neighbe" to be used in cell
        % reindexing (see later). The vertices to be suppressed should be
        % junction of the two edges :
        neighbe_forcell=[neighbe(1)*(-1)^(c(1)==2) neighbe(2)*(-1)^(c(2)==1)]
       
                
        %Add the edge between the two vertices
        tissue = addEdge(tissue,[neighbv 0 1]);
        
        % Getting the newedges array wich give the edge indexes in simplified
        % tissue relative to the vertices we've suppressed
        index = Iom_verticestosuppress == i; 
        Is_newedges(index) = length(tissue.e); % Corresponding to the vertex we've suppressed, the index of the edge we've created is the length of the current edge array     
        neighbe=sort(neighbe,'descend'); % Just for the next step to work properly
        for j = 1:length(neighbe)
                Is_newedges(Is_newedges==neighbe(j)) =  length(tissue.e); % The only line which I don't get...
                Is_newedges(Is_newedges>neighbe(j)) =  Is_newedges(Is_newedges>neighbe(j)) -1; % because of the suppression of the two neighbor edges in the next step (see removeVertex)

        end
        
        
        % Create the two new cells, with the brand new edges :
        nc = neighbors(tissue,neighbe(1));
        for j=nc
            if ~isempty(find(tissue.c{j}==neighbe_forcell(1),1)) % Case 1, the neighbe_forcell are in the good order
                ind1 = find(tissue.c{j}==neighbe_forcell(1),1);
                oldcell = circshift(tissue.c{j},[0 -ind1+1]);
                newcell = [length(tissue.e) oldcell(3:end)];
            elseif ~isempty(find(tissue.c{j}==-neighbe_forcell(2),1))
                ind2 = find(tissue.c{j}==-neighbe_forcell(2),1);
                oldcell = circshift(tissue.c{j},[0 -ind2+1]);
                newcell = [-length(tissue.e) oldcell(3:end)];
            end 
            tissue.c{j}=newcell;
            if isempty(newcell), ['bug pour la cellule' num2str(nc) 'de l''edge' num2str(neighbe(1))], end
        end
        % Suppress the two points, thereby the two "wrong" cells
        tissue = removeVertex(tissue,i); 
    end
    close(bar);
    simplifiedtissue = tissue;
    h=replot(simplifiedtissue,h);
    save([file(1:end-4) '-treated.mat'],'originaltissue','mergedtissue','simplifiedtissue','Iom_verticestosuppress','Is_newedges');
end

%% Analyse to compute the theoritical division planes

if analyse
    load([file(1:end-4) '-treated.mat']);
    tissue = simplifiedtissue;  % Retrieve the simplified version of the tissue
    ncells = length(originaltissue.c)-length(tissue.c);  % Determines the number of cells to be treated
    planes = divisionplanes(tissue,length(tissue.c)-ncells+1:length(tissue.c)); % Compute the division planes
    
    idealtissue = divide(tissue,planes(:,1)); % Create an ideal divided tissue
    save([file(1:end-4) '-treated.mat'],'originaltissue','mergedtissue','simplifiedtissue','idealtissue','Iom_verticestosuppress','Is_newedges','planes');
end


%% Compare the computed planes to the real one - topological method

if compare
    load([file(1:end-4) '-treated.mat']);
    figure('Position',[20 20 1000 500])
    axes('OuterPosition',[0 0 1/2 1]);
    plot(originaltissue);
    axes('OuterPosition',[1/2 0 1/2 1]);
    h = plot(simplifiedtissue);
    plot(idealtissue,'--');
    
    tissue = simplifiedtissue;
    
    ncells = length(originaltissue.c)-length(tissue.c); 
    real_edges = sort(Is_newedges,2);
    
    cell_list = 1:ncells;
    divisionmodes = NaN(1,ncells);
    r  = NaN(1,ncells);
    for i=1:ncells
        ideal_edges =sort(abs([tissue.c{planes(i,1).cell}([planes(i,:).i]); tissue.c{planes(i,1).cell}([planes(i,:).j])]),1)'; %Sort the indexes of the two edges of the computed planes and create an array of planesx2 indexes (see the "'")
        [a j] = intersect(ideal_edges,real_edges(i,:),'rows'); % Compare the edges of the ideal division and real division. a give the couple of edges, b the index of the division plane.
        if ~isempty(j)
            divisionmodes(i)=j;
            r(i) = (planes(i,j).walllength-planes(i,1).walllength)/planes(i,1).walllength; % Compute the probability weight
        end
    end
    
    b=[0 0 1];
    r=[1 0 0];
    
    set(h.cells([planes(divisionmodes==1,1).cell]),'FaceColor',b);%,'FaceAlpha',.5);
    set(h.cells([planes(divisionmodes==2,2).cell]),'FaceColor',2/3*b+1/3*r);%,'FaceAlpha',.5);
    set(h.cells([planes(divisionmodes==3,2).cell]),'FaceColor',1/3*b+2/3*r);%,'FaceAlpha',.5);
    set(h.cells([planes(divisionmodes>=4,2).cell]),'FaceColor',r);%,'FaceAlpha',.5);
    
    % Treatment
    A =  abs(cellArea(tissue,1:ncells))';
    D1 = (divisionmodes==1);
    D2 = (divisionmodes==2);
    D3 = (divisionmodes==3);
    I =  arrayfun(@(x) isempty(x.walllength),planes(:,3))';

    l1 = [planes(:,1).walllength];
    l2 = [planes(:,2).walllength];
    l3 = [planes(:,3).walllength];
    delta12 = 2*(l2-l1)./sqrt(A);
    delta13 = 2*(l3-l1(~I))./sqrt(A(~I));
    delta23 = 2*(l3-l2(~I))./sqrt(A(~I));

    delta12= delta12(D1|D2);
    D12 = D1(D1|D2);

    N=50;
    [delta12 index] = sort(delta12);
    D12 =D12(index);
    dl12 = zeros(floor(length(delta12)/N),1);
    ddl12 = zeros(size(dl12));
    r12 = zeros(size(dl12));
    dr12 = zeros(size(dl12));
    for i =1:floor(length(delta12)/N)
        dl12(i) = mean(delta12((i-1)*N+1:i*N));
        ddl12(i) = std(delta12((i-1)*N+1:i*N));
        r12(i) = mean(D12((i-1)*N+1:i*N));
    end

    delta13= delta13(D1(~I)|D3(~I));
    d1 = D1(~I);
    D13 = d1(D1(~I)|D3(~I));

    [delta13 index] = sort(delta13);
    D13 =D13(index);
    dl13 = zeros(floor(length(delta13)/N),1);
    ddl13 = zeros(size(dl13));
    r13 = zeros(size(dl13));
    dr13 = zeros(size(dl13));
    for i =1:floor(length(delta13)/N)
        dl13(i) = mean(delta13((i-1)*N+1:i*N));
        ddl13(i) = std(delta13((i-1)*N+1:i*N));
        r13(i) = mean(D13((i-1)*N+1:i*N));
    end

    delta23= delta23(D2(~I)|D3(~I));
    d2 = D2(~I);
    D23 = d2(D2(~I)|D3(~I));

    N=30;
    [delta23 index] = sort(delta23);
    D23 =D23(index);
    dl23 = zeros(floor(length(delta23)/N),1);
    ddl23 = zeros(size(dl23));
    r23 = zeros(size(dl23));
    dr23 = zeros(size(dl23));
    for i =1:floor(length(delta23)/N)
        dl23(i) = mean(delta23((i-1)*N+1:i*N));
        ddl23(i) = std(delta23((i-1)*N+1:i*N));
        r23(i) = mean(D23((i-1)*N+1:i*N));
    end
    
    figure
    hold on
    msize =  10;
    deltamax = 1;
    errorbar2(dl12,ddl12,r12,dr12,'ok','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',msize);
    errorbar2(dl13,ddl13,r13,dr13,'ok','MarkerEdgeColor','k','MarkerFaceColor','none','MarkerSize',msize);
    errorbar2(dl23,ddl23,r23,dr23,'sk','MarkerEdgeColor','k','MarkerFaceColor','none','MarkerSize',msize);
    deltafit = 0:.01:1;
    rfit = 1./(1+exp(-10*deltafit));
    plot(deltafit,rfit,'-k');
end
    
