function [ s1 ] = findbounds( M,i,j )
%FINDBOUNDS find the boundaries of the otimization performed by finds.
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu

s=zeros(1,2);
s2=s;
s1=s;
[ ~, s2(1) ]=finds2(0,M,i,j);
[ ~, s1(1) ]=finds2(1,M,j,i);
[ ~, s2(2) ]=finds2(1,M,i,j);
[ ~, s1(2) ]=finds2(0,M,j,i);
if all((s1>=0).*(s1<=1)), 
    %display('debug : fulls2');
else
    if  s2(2)>s2(1), s1=circshift(s1, [0,1]); end % Inversion of the abscissas case
    s1=[0,1]+(s1-[0,1]).*(s2>1 | s2<0);
end
if s1(1)<0, s1=[0,0]; return; end % no intersection exception

[ ~, ~, a1 ~, ~]=finda( s1(1), M, i, j );
[ ~, ~, a2 ~, ~]=finda( s1(2), M, i, j );
if a1<a2, display('debug : shift'); s1=circshift(s1,[1,0]); end
if ~(a1>0.5 && a2<0.5), s1=[]; end % no area for constraint exception
    
end

