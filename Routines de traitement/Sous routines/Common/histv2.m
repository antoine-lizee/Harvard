function [no,xo] = histv2(varargin)
%HIST V2 "Window" Histogram
%This function takes the syntax from the built in function "hist", but 
%enables to have a "sliding window" count instead of a simple histogram.
%Indeed, you can specify, after the first integer which denotes the window
%width, a second integer which let you fix a number of points you want to
%be plotted.

%For instance, histv3(Y,a,b) will give you b bins, but of a width equal to
%the bins which would have been considered with the command hist(Y,a).

%This function does not support all the features of the built in hist
%function. In particular :
% - you cannot specify a vector for the bins
% - the easy plot abilities are not implemented (but you should prefer the
% classic "plot" function instead of the "bar" function which is used by 
%the built-in "hist" function for easy plot).

%Remark : Because of the renormalization undergone for points whose
%window goes beyond the range of the data sample, the sum of the points (in "no"), even
%if divided by the ratio b/a, will not exactly be equal to the initial
%amount of point initially present in the Y vector.

%See also : hist, hist2, Test_histv2

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

%Transform the NN matrix in order to get a full window at the extremities
n0=zeros(wx/2-wz/2,size(NN,2));
NN=[n0; NN; n0];

nn=NN(1:wz:end-wx+1,:);
xx=XX(1:wz:end-wz+1);
for i=2:wx
    nn=nn+NN(i:wz:end-wx+i,:);
end
for i=2:wz
    xx=xx+XX(i:wz:end-wz+i);
end
xx=xx/wz;

% Normalize the extreme values
hpi=ceil(wx/wz); %Localize the first half-little-bin (inherited from "b") which go past the first half-wide-bin (from"a")
pi=floor(hpi/2)+1; %Localize the first little bin whose middle go beyond the middle of the first wide bin.
for i=1:pi-1
    nn(i,:)=nn(i,:)*wx/(wx/2+(2*i-1)*wz/2);
    nn(end-i+1,:)=nn(end-i+1,:)*wx/(wx/2+(2*i-1)*wz/2);
end

if nargout > 0
    no=nn;
    xo=xx;
else
  error('histv2 does not support easy plotting abilities, you have to specify an output argument');
end

function a = ishistnumeric(b)
% for backward compatibility, logical is allowed in hist.m
a = isnumeric(b) || islogical(b);

