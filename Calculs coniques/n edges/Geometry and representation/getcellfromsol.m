function [ cell ] = getcellfromsol( V, A )
%GETCELLFROMSOL Summary of this function goes here
%   Create the CellNetwork from the points and angles output of the
%   'findfixedpoint' function.
%
%   Copyright Antoine Lizée July 2011 
%   @Dumais lab, OEB dpt, Harvard University.

[V1 V2]=pol2cart(V(:,1),-V(:,2)); %Put the gamma angle of the cone on the left side during the change in coodinates
V=[[0,0];[V1 V2];[0,0]];
E=[(1:(length(V)-1))',(2:length(V))', [0; A'; 0]];
C={1:(length(V))-1};
cell=cellNetwork(V,E,C);

end

