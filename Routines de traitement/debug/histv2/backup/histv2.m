function [no,xo] = histv2(varargin)
%HIST V2 "Window" Histogram
%This function takes the syntax from the built in function "hist", but 
%enables to have a "sliding window" count instead of a simple histogram.
%Indeed, you can specify, after the first integer which denotes the window
%width, a second integer which let you fix a number of points you want to
%be plotted.

%For instance, histv2(Y,a,b) will give you b bins, but of a width equal to
%the bins which would have been considered with the command hist(Y,a).

%This function does not support all the features of the built in hist
%function. In particular :
% - you cannot specify a vector for the bins
% - the easy plot abilities are not implemented (but you should prefer the
% classic "plot" function instead of the "bar" function used by the built
% in "hist" for easy plot.

%Nota Bene : The number of point calculated is less then the integer and
%input argument "b", because the extreme right and left point cannot be
%calculated properly, since the bin of width inherited from "a" would exceed
%the range of the supplied vector. The exact amount of point is :
%n=ceil(b*(1-1/a))

%See also : hist, hist2

%   Copyright 2011 Antoine Lizée 
%   $Revision: 1.0 $Date: 2011/05/18$


% Parse possible Axes input %Here, useless concerning the axes, but we use args
% and nargs variables for conveniance :
error(nargchk(1,inf,nargin,'struct'));
[cax,args,nargs] = axescheck(varargin{:});

y = args{1};

if nargs == 1
    x = 10;
else
    x = args{2};
end

if nargs == 3
    z = args{3};
else
    z = x;
end

if isvector(y), y = y(:); end
if ~ishistnumeric(x) || ~ishistnumeric(y)
    error('MATLAB:histv2:InvalidInput', 'Input arguments must be numeric.')
end
if length(x)*length(z)~=1
    error('MATLAB:histv2:InvalidInput', 'Input arguments for bin definition in histv2 must be integers')
end

w=2*lcm(x,z);
wz=w/z;
wx=w/x;
[NN XX]=hist(y,w); 
if isvector(y) && ~isempty(y) % Flip row vectors in column
    NN=NN';
    XX=XX';
end
pi=ceil(wx/wz); %Get the first point with full window
indi=pi*wz/2-wx/2;%Get the indice of the first "underbin" going in the first bin
indf=pi*wz/2+wx/2;
nn=NN(indi+1:wz:end-indf+1,:);
xx=XX(indi+1:wz:end-indf+1);
for i=2:wx
    nn=nn+NN(indi+i:wz:end-indf+i,:);
    xx=xx+XX(indi+i:wz:end-indf+i);
end
xx=xx/wx;

if nargout > 0
    no=nn;
    xo=xx;
else
  error('histv2 does not support easy plotting abilities, you have to specify an output argument');
end

function a = ishistnumeric(b)
% for backward compatibility, logical is allowed in hist.m
a = isnumeric(b) || islogical(b);

