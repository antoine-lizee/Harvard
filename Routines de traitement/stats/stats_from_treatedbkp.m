function [ statdata ] = stats_from_treated( treated )
%%UNTITLED3 Summary of this function goes here
%Enable one to extract from "-treated.mat" files the data needed for
%statistical studies : computed planes length and choosen plane by the
%cell.
%%   Detailed explanation goes here
%The output arg is a 1xN struct for N divided cells with four fields : the
%first "planesnumber" is an array with the number of computed minimal
%planes, the second field "planes", is an array with the length of the
%computed planes, the third one, "divisionplane" is an integer which
%indicates which one has eventually been chosen by the cell. 
%We use only the planes.mat and originaltissue.mat
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

%% Assigning the different fields
lp=size(planes);
statdata(lp)=struct('planesnumber',[], 'planes',[], 'divisionplane', [], 'SI_planes', []); %preallocation of statdata

for i=1:lp
    statdata(i).planesnumber=length([planes(i,1:end).walllength]);
    statdata(i).planes=[planes(i,1:end).walllength];
    statdata(i).SI_planes=planes(i);
end

%% Calculation of the division planes indexes

%eventuellement, recalculer les variables utilisée avec des bouts de
%generic tissue analysis, puis faire le bout de truc ci-dessous. Ne pas
%oublier de sauver plus que cela à la toute fin !


real_edges = sort(edges_list,2);

for i=1:ncells
        ideal_edges =sort(abs([tissue.c{planes(i,1).cell}([planes(i,:).i]); tissue.c{planes(i,1).cell}([planes(i,:).j])]),1)'; %Sort the indexes of the two edges of the computed planes and create an array of planesx2 indexes (see the "'")
        [a j] = intersect(ideal_edges,real_edges(i,:),'rows'); % Compare the edges of the ideal division and real division. a give the couple of edges, b the index of the division plane.
        if ~isempty(j)
            statdata(i).divisionplane=j;
            r(i) = (planes(i,j).walllength-planes(i,1).walllength)/planes(i,1).walllength; % Compute the probability weight
        end
end

%% Reload the latter workspace
clear
load('workspace_temp')
delete('workspace_temp.mat')

end

