%{
for i= nv(end:-1:1)
    [neighbe c]=find(mergedtissue.e(:,1:2)==i); %find the position of the vertex index i in the edge array
    if length(neighbe)~=2, error('More or less than 2 edges neighbor to the vertex to suppress'); end
    c=(1+(~logical(c-1)))' %invert the column between 1 and 2
    Iom_neighbv=[mergedtissue.e(neighbe(1),c(1)), mergedtissue.e(neighbe(2),c(2))]; % compute the vertices index
    Is_neighbv=find((simplifiedtissue.v(:,1)==mergedtissue.v(Iom_neighbv(1),1)) & (simplifiedtissue.v(:,1)==mergedtissue.v(Iom_neighbv,2)));
    
    %{
get the Iom_neighbv, find their coordinates, find Is_neighbv, get Is_newedges.
    Or, do it his way.
    %}
end
  %}  
%{
test=-100*ones(length(simplifiedtissue.c),1);
for j=1:length(simplifiedtissue.c)
for i=1:length(simplifiedtissue_bis.c),
if isequal(simplifiedtissue_bis.c(i), simplifiedtissue.c(j)), test(j)=i; i, end
end
end
test
%}

%{
tissue=mergedtissue;

for i=nv(end) % i indice of the vertex to be removed, initially in mergedtissue, but also in tissue thanks to the backward looping
                
        % Find the neighbors of the vertex
        [neighbe c]=find(tissue.e(:,1:2)==i); % locate the indexes in a boolean array); %find the position of the vertex index i in the edge array (neighbe is therefore the two neighbor edge indexes)
        if length(neighbe)~=2, error('More or less than 2 edges neighbor to the vertex to suppress'); end
        c=1+(c==1); %invert the columns
        neighbv=[tissue.e(neighbe(1),c(1)) tissue.e(neighbe(2),c(2))]  % compute the vertices index of the neighbor vertices
        
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
%}

% Draw an array for cell consistency debugging

A=zeros(length(simplifiedtissue.c),18);
for i=1:length(simplifiedtissue.c)
    for j=1:length(simplifiedtissue.c{i})
        if simplifiedtissue.c{i}(j)>0
            A(i,[(2*j-1),(2*j)])=simplifiedtissue.e(simplifiedtissue.c{i}(j),1:2);
        else
            A(i,[(2*j-1),(2*j)])=simplifiedtissue.e((-simplifiedtissue.c{i}(j)),[2,1]);
        end
    end
end

