function [ statdata ] = stats_from_treated( treated )
%%STATS_FROM_TREATED Summary of this function goes here
%Enable one to extract from "-treated.mat" files the data needed for
%statistical studies.
%%   Detailed explanation goes here
%The output arg is a 1xN struct for N divided cells with five fields : the
%first "planesnumber" is the number of computed minimal
%planes, the second field "planes", is an array with the length of the
%computed planes, the third one, "divisionmode" is an integer which
%indicates which one has eventually been chosen by the cell. 
%We use only the planes.mat, originaltissue.mat, and eges_list variable
%
%Copyright Antoine Lizée @ Dumais Lab, 2011

%% Errors detection & workspace finding then loading
if nargin>1, 
    error('Wrong number of input arguments : zero or one argument is required')
end

if nargin==0, file_path=uigetfile;
else file_path=treated;
end

save('workspace_temp')
clear
load(file_path)

%% Assigning the different fields inherited from 'planes'
lp=size(planes);
statdata(lp)=struct('planesnumber',[], 'planes',[], 'divisionmode', [], 'SI_planes', []); %preallocation of statdata

for i=1:lp
    %redirection of the fields already enclosed in the "planes" struct
    statdata(i).planesnumber=length([planes(i,1:end).walllength]);
    statdata(i).planes=[planes(i,1:end).walllength];
    statdata(i).SI_planes=planes(i);
end
    
%% Calculation of the division modes

real_edges = sort(edges_list,2);

for i=1:lp
        ideal_edges =sort(abs([simplifiedtissue.c{planes(i,1).cell}([planes(i,:).i]); simplifiedtissue.c{planes(i,1).cell}([planes(i,:).j])]),1)'; %Sort the indexes of the two edges of the computed planes and create an array of planesx2 indexes (see the "'")
        [a j] = intersect(ideal_edges,real_edges(i,:),'rows'); % Compare the edges of the ideal division and real division. a give the couple of edges, b the index of the division plane.
        if ~isempty(j)
            statdata(i).divisionmode=j;
            r(i) = (planes(i,j).walllength-planes(i,1).walllength)/planes(i,1).walllength; % Compute the probability weight
        end
end

%% Maping the vertex index from 'simplifiedtissue' (therefore from 'planes') to the index of 'original- or 'mergedtissue 
%% Creating the array

for i=1:lp
    statdata_array(i)=[statdata(i).reallength, statdata(i).planesnumber, statdata(i).divisionmode, [statdata(i).planes] ];
end


%%Exporting and save of the different data

save


%% Reload the latter workspace
clear
load('workspace_temp')
delete('workspace_temp.mat')

end

