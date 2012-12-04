function [ sa sb af V A] = finds( s11,s12,a,M,i,j,p )
%FINDS carry out the optimization neede to provide the dividing wall that
%will cut the mother cell in half
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu

[s1 s2 a3 V A]=finda((s11+s12)/2,M,i,j);
%disp(['debug:    ' num2str(s1) num2str(a3)]);
if abs(s1-s12)<=1/p
    sa=s1;
    sb=s2;
    af=a3;
    return
end

if a3<=a
    [ sa sb af V A]=finds( s11,s1,a,M,i,j,p );
elseif a3>a
    [ sa sb af V A]=finds( s1,s12,a,M,i,j,p );
end
    
end

