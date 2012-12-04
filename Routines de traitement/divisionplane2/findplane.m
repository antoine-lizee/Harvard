function [ sa sb V A af length] = findplane( M, i, j, a )
%FINDPLANE is the function that seek the divding wall that fulfill Errera's
%rules between two given edges of one cell, if this wall exists.
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu

if nargin<4, a=0.5; end

%% full straight edges computation



%% half straight check



%% general computation 
s1=findbounds(M,i,j);

if isempty(s1), 
    A=[];
    V=[];
    af=[];
    sa=[];
    sb=[];
    length=[];
    return;
elseif s1==[0, 0],
    A=[];
    V=[];
    af=[];
    sa=1;
    sb=[];
    length=[];
    return;
end

[sa, sb, af, V, A]=finds(s1(1), s1(2), a, M, i, j, 1000/range(s1));

%Computation of the length
length = norm(V(1,:)-V(2,:))*A(3)/(sin(A(3)));

end

