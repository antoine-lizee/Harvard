 tissue=mergedtissue;
for i=nv(end) % i indice of the vertex to be removed
        
        i1 = find(tissue.e(:,1)==i); % Locate the edge indexes using the vertex to be removed
        i2 = find(tissue.e(:,2)==i); % i1, i2 are vector of edge indexes
        if length([i1;i2])~=2, error('More than 2 edges'); end
        
        % Locate the neigbor vertices of the vertex i %completely
        % redondant with the following part. v is unused.
        allv = unique(tissue.e([i1;i2],1:2));
        v = allv(allv~=i);

        % Order edges to find them in the cell array in the good order (see
        % below)
        if isempty(i1) %then i2 has two rows of one edge index each
            edges = [i2(1) -i2(2)];
            v1 = tissue.e(i2(1),1); %v1, v2 are vertex indexes of the neighbor vertices of vertex i.
            v2 = tissue.e(i2(2),1); %redondant with v
        elseif isempty(i2)
            edges = [-i1(1) i1(2)];
            v1 = tissue.e(i1(1),2);
            v2 = tissue.e(i1(2),2);
        else
            edges = [i2 i1];
            v1 = tissue.e(i2,1);
            v2 = tissue.e(i1,2);
        end
        [v1 v2]
        edges        
        tissue = addEdge(tissue,[v1 v2 0 1]); %Add the edge between the two vertices

        index = find(vertices_list_bis == i);
        edges_list(index) = length(tissue.e);

        % List connected edges
        e_index = find(transpose(tissue.e(:,1:2) == i));
        e_index = round(e_index/2);

        e_index = sort(e_index, 'descend');
        for j = 1:length(e_index)
                edges_list(edges_list==e_index(j)) =  length(tissue.e); 
                edges_list(edges_list>e_index(j)) =  edges_list(edges_list>e_index(j)) -1; % because of the suppression of the two neighbor edges later

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
            tissue.c{j}=newcell; % Add, for each jth neighbor cell, the cell with the new edge 
        end
        tissue = removeVertex(tissue,i); % Suppress the two points, thereby the two "wrong" cells
    end