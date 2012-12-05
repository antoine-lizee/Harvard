function [ h cell ] = draw_2D( V, A )
%DRAW_ Summary of this function goes here
%   This function enalbes the user to plot on the 2D-plan the solution 
%issued by 'findfixedpoint' for example.
%   Detailed explanation goes here
%   We do some debugging tasks here by the fact that we use outer and old
%   functions that rely on cellNetworks to plot ("plot_local" is just the plot 
%   overloaded function for cellnetwork, renamed) and compute cell area.
%
%   Copyright Antoine Lizée July 2011 
%   @Dumais lab, OEB dpt, Harvard University.

cell=getcellfromsol(V,A);
h=plot_local(cell,'-or');
axis equal

da=cellArea(cell,1)-1;
if da>0.0000000000001;
    disp(da);
end

end

