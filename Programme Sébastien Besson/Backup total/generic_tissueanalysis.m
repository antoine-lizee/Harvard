% Routine for treating an extracted tissue
% Assume you have a fully defined tissue

clear
clc
close all
warning %off all

% Name of the folder & file
[folder file]=uigetfile;

mergetissue = 1;
simplify =1;
analyse =1;
compare =1;

if mergetissue % Remove youngest edges 
    load([folder file]); % Load the file
    originaltissue = tissue; % Save the original tissue information
    check(tissue);
    tissue.e(:,3) = 0; % Remove the angle information

    N_edgestomerge = sum(tissue.e(:,4)==2); % Determine number of edges to merge
    vertices_list = zeros(N_edgestomerge,2); 
    for i=1:N_edgestomerge
        ne = find(tissue.e(:,4)==2,1);
        vertices_list(i,:) = tissue.e(ne,1:2); % Get the list of vertices 
        tissue = merge(tissue,neighbors(tissue,ne)); 
    end
    mergedtissue =tissue; % Save the merged tissue
    save([file(1:end-4) '-treated.mat'],'originaltissue','mergedtissue','vertices_list');
end

if simplify
    load([file(1:end-4) '-treated.mat']);
    figure('Position',[20 20 1000 500])
    axes('OuterPosition',[0 0 1/2 1]);
    h0 = plot(originaltissue);
    axes('OuterPosition',[1/2 0 1/2 1]);
    h = plot(mergedtissue);
    tissue = mergedtissue; % Retrieve the merged version of the tissue
    edges_list = zeros(size(vertices_list)); % Initialise the edges list
    
    old_connectivity = arrayfun(@(x) sum(sum(originaltissue.e(:,1:2)==x)),1:length(originaltissue.v));
    new_connectivity = arrayfun(@(x) sum(sum(tissue.e(:,1:2)==x)),1:length(tissue.v));
    nv = find(new_connectivity==2 & new_connectivity~=old_connectivity); % Find the vertices to be removed

    bar = waitbar(0,'Removing vertices');
    for i=nv(end:-1:1)
        waitbar(1-i/length(nv),bar,['Removing vertex: ' num2str(i)]);
        i1 = find(tissue.e(:,1)==i);
        i2 = find(tissue.e(:,2)==i);
        if length([i1;i2])~=2, error('More than 2 edges'); end

        allv = unique(tissue.e([i1;i2],1:2));
        v = allv(allv~=i);

        % Order edges to find them in the cell array
        if isempty(i1)
            edges = [i2(1) -i2(2)];
            v1 = tissue.e(i2(1),1);
            v2 = tissue.e(i2(2),1);
        elseif isempty(i2)
            edges = [-i1(1) i1(2)];
            v1 = tissue.e(i1(1),2);
            v2 = tissue.e(i1(2),2);
        else
            edges = [i2 i1];
            v1 = tissue.e(i2,1);
            v2 = tissue.e(i1,2);
        end
        tissue = addEdge(tissue,[v1 v2 0 1]);

        index = find(vertices_list == i);
        edges_list(index) = length(tissue.e);

        % List connected edges
        e_index = find(transpose(tissue.e(:,1:2) == i));
        e_index = round(e_index/2);

        e_index = sort(e_index, 'descend');
        for j = 1:length(e_index)
                edges_list(edges_list==e_index(j)) =  length(tissue.e);
                edges_list(edges_list>e_index(j)) =  edges_list(edges_list>e_index(j)) -1;

        end

        nc = neighbors(tissue,edges(1));
        for j=nc
            if ~isempty(find(tissue.c{j}==edges(1),1))
                ind1 = find(tissue.c{j}==edges(1),1);
                oldcell = circshift(tissue.c{j},[0 -ind1+1]);
                newcell = [length(tissue.e) oldcell(3:end)];
            elseif ~isempty(find(tissue.c{j}==-edges(2),1))
                ind2 = find(tissue.c{j}==-edges(2),1);
                oldcell = circshift(tissue.c{j},[0 -ind2+1]);
                newcell = [-length(tissue.e) oldcell(3:end)];
            end
            tissue.c{j}=newcell;
        end
        tissue = removeVertex(tissue,i);
    end
    close(bar);
    simplifiedtissue = tissue;
    h=replot(simplifiedtissue,h);
    save([file(1:end-4) '-treated.mat'],'originaltissue','mergedtissue','simplifiedtissue','vertices_list','edges_list');
end

if analyse
    load([file(1:end-4) '-treated.mat']);
    tissue = simplifiedtissue;  % Retrieve the simplified version of the tissue
    ncells = length(originaltissue.c)-length(tissue.c);  % Determines the number of cells to be treated
    planes = divisionplanes(tissue,length(tissue.c)-ncells+1:length(tissue.c)); % Compute the division planes
    
    idealtissue = divide(tissue,planes(:,1)); % Create an ideal divided tissue
    save([file(1:end-4) '-treated.mat'],'originaltissue','mergedtissue','simplifiedtissue','idealtissue','vertices_list','edges_list','planes');
end

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
    real_edges = sort(edges_list,2);
    
    cell_list = 1:ncells;
    divisionmodes = NaN(1,ncells);
    r  = NaN(1,ncells);
    for i=1:ncells
        ideal_edges =sort(abs([tissue.c{planes(i,1).cell}([planes(i,:).i]); tissue.c{planes(i,1).cell}([planes(i,:).j])]),1)';
        [a j] = intersect(ideal_edges,real_edges(i,:),'rows');
        if ~isempty(j)
            divisionmodes(i)=j;
            r(i) = (planes(i,j).walllength-planes(i,1).walllength)/planes(i,1).walllength;
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
    
