function [ s1 ] = findboundsii( M,i,j )
%FINDBOUNDSII is the version of find bounds for same edged case.
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu
if i~=j
    error('findboundsii must be called only for division on the same edge (i=j)')
end

s1=[0,0.5];

[ ~, ~, a1 ~, ~]=findaii( s1(1), M, i, j );
if ~(a1>0.5),  s1=[]; end % no area for constraint exception
    
end

