function [ sa sb af V A] = findsii( s11,s12,a,M,i,j,p )
%FINDS Summary of this function goes here
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu
if i~=j
    error('findsii must be called only for division on the same edge (i=j)')
end

[s1 s2 a3 V A]=findaii((s11+s12)/2,M,i,j);
%disp(['debug:    ' num2str(s1) num2str(a3)]);
if abs(s1-s12)<=1/p
    sa=s1;
    sb=s2;
    af=a3;
    return
end

if a3<=a
    [ sa sb af V A]=findsii( s11,s1,a,M,i,j,p );
elseif a3>a
    [ sa sb af V A]=findsii( s1,s12,a,M,i,j,p );
end
    
end

