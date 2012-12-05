function hh = errorbar2(x,dx,y,dy,symbol,varargin)
%ERRORBAR2 Error bar plot for both axes based on the errorbar function.
% ERRORBAR2(X,DX,Y,DY) plots the graph of vector X versus vector Y with
%  error bars in both directions specified by the vectors DX and DY.
% ERRORBAR2(X,DX,Y,DY) plots X with error bars [X-DX X+DX] and 
%  plots Y with error bars [Y-DY Y+DY].
% ERRORBAR2(...,'LineSpec') uses the color and linestyle specified by
%  the string 'LineSpec'.  See PLOT for possibilities.
% H = ERRORBAR2(...) returns a vector of line handles.
%
%  Marcos Duarte  mduarte@usp.br Sept1998 

if min(size(x))==1,
  npt = length(x);
  x = x(:);
  dx = dx(:);
  y = y(:);
  dy = dy(:);
else
  [npt,n] = size(x);
end
if nargin == 4
   symbol = '-';
%elseif nargin ~= 5
%   error('Number of arguments msut be 4 or 5.')
end
    
if isstr(x) || isstr(y) || isstr(dx) || isstr(dy)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) || ~isequal(size(x),size(dx))...
        || ~isequal(size(x),size(dy)),
  size(x)
  size(dx)
  size(y)
  size(dy)
    error('The sizes of X, Y, DX and DY must be the same.');
  
end

% Plot graph and bars
hold_state = ishold;
cax = newplot;
next = lower(get(cax,'NextPlot'));

dx = abs(dx(:));
dy = abs(dy(:));

tee = dx/10;
ell = dy/10;

xleft = x - dx;
xright = x + dx;
xltee = x - tee;
xrtee = x + tee;

ytop = y + dy;
ybot = y - dy;
ybell = y - ell;
ytell = y + ell;

vxb = zeros(npt*9,1);
vxb(1:9:end) = x;
vxb(2:9:end) = x;
vxb(3:9:end) = NaN;
vxb(4:9:end) = xltee;
vxb(5:9:end) = xrtee;
vxb(6:9:end) = NaN;
vxb(7:9:end) = xltee;
vxb(8:9:end) = xrtee;
vxb(9:9:end) = NaN;

vyb = zeros(npt*9,1);
vyb(1:9:end) = ytop;
vyb(2:9:end) = ybot;
vyb(3:9:end) = NaN;
vyb(4:9:end) = ytop;
vyb(5:9:end) = ytop;
vyb(6:9:end) = NaN;
vyb(7:9:end) = ybot;
vyb(8:9:end) = ybot;
vyb(9:9:end) = NaN;

hxb = zeros(npt*9,1);
hxb(1:9:end) = xleft;
hxb(2:9:end) = xright;
hxb(3:9:end) = NaN;
hxb(4:9:end) = xleft;
hxb(5:9:end) = xleft;
hxb(6:9:end) = NaN;
hxb(7:9:end) = xright;
hxb(8:9:end) = xright;
hxb(9:9:end) = NaN;

hyb = zeros(npt*9,1);
hyb(1:9:end) = y;
hyb(2:9:end) = y;
hyb(3:9:end) = NaN;
hyb(4:9:end) = ybell;
hyb(5:9:end) = ytell;
hyb(6:9:end) = NaN;
hyb(7:9:end) = ybell;
hyb(8:9:end) = ytell;
hyb(9:9:end) = NaN;

[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

h = plot(vxb,vyb,esymbol,varargin{:}); hold on
h = [h;plot(hxb,hyb,esymbol,varargin{:})]; hold on
h = [h;plot(x,y,symbol,varargin{:})]; 

if ~hold_state, hold off; end

if nargout>0, hh = h; end
