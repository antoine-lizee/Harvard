% This script compute the values needed to draw the solution in geometrical
% softwares like geogebra. The solution is normalized to the first length
% from the edge to the center (h1).
%
%Copyright Antoine Lizée @ Dumais Laboratory, Harvard University
%July 2011

n=3;

gamma0=0:0.1:266;
rg=zeros(length(gamma0),4);

for i=1:length(gamma0)
    [~, rg(i,:)]=fixedpoint(n,gamma0(i),1);
end

eval(['results_for_' num2str(n) '_edges=rg;']);