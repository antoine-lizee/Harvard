function [ output_args ] = divisionplane2( N,n,a )
%DIVISIONPLANE2 is not finished !!
%
%See documentation to understand further the role of the function in the
%entire program divisionplanes. This documentation is available from the
%author or the corresponding professor Jacques Dumais.
%
%Copyright Antoine Lizée @ Dumais Laboratory, OEB Dpt of Harvard University
%June 2011 email :antoine.lizee@polytechnique.edu

% Check the number of inputs
error(nargchk(1, 3, nargin));

% If no optional argument, loops over the cells
if (nargin < 3), ratio = .5; end;
if (nargin < 2), n = 1:length(N.c); end;

% If empty second argument
 if isempty(n), planes = []; return; end;

planes = struct('cell',0,'i',0,'j',0,'s1',0.0,'s2',0.0,'P',[],'Q',[],...
    'angle',0.0,'walllength',0.0,'probability',0.0);



for c = 1:length(n)
    nc = n(c);

    M=findcarac(N,nc);
    
    for i = 1:size(M,1)-1
        for j = i:size(M,1)-1
            
            finds( s1,s2,a,p,Me )
            
        end
    end
     
    
end



end

