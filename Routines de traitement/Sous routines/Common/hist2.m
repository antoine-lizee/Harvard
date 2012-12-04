function [no,xo] = hist2(varargin)
%%HIST 2 "Softened" Histogram

%SUMMARY:
%This function takes all the abilities of the built-in "hist" function,
%but it computes the numbers of points for bins that are different 
%(somehow wider)than the original bins. It enables the user to compute an
%histogram that get rid of the fast variations. See MEANING for further
%explanation.

%SYNTAX :
%- N=hist2(Y,M) is the same as hist2(Y,M,1)
%- N=hist2(Y,M,z), Y is a vector of points or an array. M is the number of 
%initial bins, and z is the level of "softening". z is also the number of 
%neighbor bins taken into account to calculate the value for each initial 
%bin (see PRECISIONS).
%- N=hist2(Y,X,...), X is a vector specifying bin centers (see "hist").
%Although this function is supported, you are not advised to use it.
%- [N X]=hist2(...) is a supported syntax (see "hist")
%- hist2(...) is a supported syntax (see "hist")

%PRECISIONS
%The value n computed for each bin is the result of a weighting of the 
%neighbor bins depending on the third arg "z" as follows :
%
% z                               plotted bin     3rd right neighbor
%|||                                 |||               |||
% 0     0     0     0     0     0     1     0     0     0     0     0
% 1(def)0     0     0     0     1     1     1     0     0     0     0
% 2     0     0     0     1     2     3     2     1     0     0     0
% 3     0     0     1     3     6     7     6     3     1     0     0
% 4     0     1     4    10    16    19    16    10     4     1     0
% 5     1     5    15    30    45    51    45    30    15     5     1
% etc...
% It is then normalized by the sum of the coefficients used for each bin,
%in order to unchange the total number of points. The computation of the
%extremities of the sample are slightly different.

%THESE coefficients, for z>1, are a very close approximation of a
%normalized gaussian,centered at 0 and of standard deviation equal to :
% sigma=sqrt(2/3*z).

% MEANING : 
%For the first step (z=1), hist2 is a M-point histogram with 
%(range(Y)/M*3)-wide bins (three times larger than the initial ones).
% (See "histv2")
%For any further step, the weigthing carried out by hist2 can have the two 
%following precise meanings :
% - from a probability point of view, it computes the distribution of the
%law resulting from the sum of the initial law corresponding to the classic
%histogram and a normal law, centered in zero and of std sigma. This
%represents a quantification of the error due to discretisation of the
%probability distribution.
% - from an analytic point of view, it computes the convolution of the
% distribution function of the sample with a gaussian function, centered 
%in zero and of standard deviation sigma. As a result, we got  the 
%distribution function in which we have gradually and smoothly erased all
%variations which have a half-period smaller than the typical half period :
% T/2=pi*STD=pi*sigma.
%We put here these values versus some z values :

%     z         | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 15| 22 | 30 | 39 | 61 | 95
%--------------------------------------------------------------------------------
%typical half-  |3.6|4.4|5.1|5.7|6.3|6.8|7.3|7.7|8.1|9.9|12.0|14.0|16.0|20.0|25.0
%period of rejection (in bins)  

% SEE ALSO : hist, histv2, Test_hist2v2

%   Copyright 2011 Antoine Lizée 
%   $Revision: 1.0 $Date: 2011/05/10$


%%Code taken mainly from the "hist" built in function
% Parse possible Axes input
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
    z = 1;
end

if isvector(y), y = y(:); end
if ~ishistnumeric(x) || ~ishistnumeric(y)
    error('MATLAB:hist:InvalidInput', 'Input arguments must be numeric.')
end

% Cache the vector used to specify how bins are created
N = x;
    
if isempty(y),
    if length(x) == 1,
       x = 1:double(x);
    end
    nn = zeros(size(x)); % No elements to count
    %  Set miny, maxy for call to bar below.
    miny = [];
    maxy = [];
    edges = [-Inf Inf];
else
    %  Ignore NaN when computing miny and maxy.
    ind = ~isnan(y);
    miny = min(y(ind));
    maxy = max(y(ind));
    %  miny, maxy are empty only if all entries in y are NaNs.  In this case,
    %  max and min would return NaN, thus we set miny and maxy accordingly.
    if (isempty(miny))
      miny = NaN;
      maxy = NaN;
    end
    
    if length(x) == 1
    	  if miny == maxy,
    		  miny = miny - floor(x/2) - 0.5; 
    		  maxy = maxy + ceil(x/2) - 0.5;
     	  end
        binwidth = (maxy - miny) ./ x;
        xx = miny + binwidth*(0:x);
        xx(length(xx)) = maxy;
        x = xx(1:length(xx)-1) + binwidth/2;
    else
        xx = x(:)';
        binwidth = [diff(xx) 0];
        xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
        xx(1) = min(xx(1),miny);
        xx(end) = max(xx(end),maxy);
    end
    % Shift bins so the interval is ( ] instead of [ ).
    xx = full(real(xx)); y = full(real(y)); % For compatibility
    bins = xx + eps(xx);
    edges = [-Inf bins];
    nn = histc(y,edges,1);
    edges(2:end) = xx;    % remove shift
    
    % Combine first bin with 2nd bin and last bin with next to last bin
    nn(2,:) = nn(2,:)+nn(1,:);
    nn(end-1,:) = nn(end-1,:)+nn(end,:);
    nn = nn(2:end-1,:);

%% Here is the main modification :
    for i=1:z
        nn2=nn(2,:);
        nne=nn(end-1,:);
        nn(2:end-1,:)=(nn(1:end-2,:)+nn(2:end-1,:)+nn(3:end,:))/3;
        nn(1,:)=(2*nn(1,:)+nn2)/3;
        nn(end,:)=(2*nn(end,:)+nne)/3;
    end
    edges(2) = [];
    edges(end) = Inf;
end

%%

if nargout > 0
    if isvector(y) && ~isempty(y) % Return row vectors if possible.
        no = nn';
        xo = x;
    else
        no = nn;
        xo = x';
    end
else
  if ~isempty(cax)
     histPatch = bar(cax,x,nn,[miny maxy],'hist');
     if strcmpi(get(ancestor(cax,'figure'),'Visible'),'off')
         return;
     end     
  else
     histPatch =  bar(x,nn,[miny maxy],'hist');
  end
  
  % Add linked plot and brushing behavior objects to the patch
  varName = inputname(1+~isempty(cax));
  for k=1:length(histPatch)
      linkBehavior = hggetbehavior(histPatch(k),'Linked');
      linkBehavior.DataSourceFcn = {@localDataChange linkBehavior N};
      if ~isempty(varName)
          if ~isvector(y)
              linkBehavior.DataSource = getcolumn(varName,k,'expression');
          else
              linkBehavior.DataSource = varName;
          end
      end
      linkBehavior.BrushFcn = {@localBrushFunc linkBehavior x};
      linkBehavior.LinkBrushFcn = {@localSetLinkedBrushing};
      brushBehavior = hggetbehavior(histPatch(k),'Brush');
      brushBehavior.DrawFcn = {@localDrawFunc brushBehavior histPatch(k)};
      if isempty(get(histPatch(k),'DisplayName'))
          set(histPatch(k),'DisplayName',linkBehavior.DataSource)
      end
      datacursorBehavior = hggetbehavior(histPatch(k),'datacursor');
      set(datacursorBehavior,'UpdateDataCursorFcn',...
          {@localHistUpdateDataCursor,histPatch(k)});
      set(datacursorBehavior,'MoveDataCursorFcn',...
          {@localHistMoveDataCursor,histPatch(k)});
      set(datacursorBehavior,'UpdateFcn',...
          {@localHistDatatipCallback,x,nn(:,k),edges});
  end
  
  % If applying to a linked plot the linked plot graphics cache must
  % be updated manually since there are not yet eventmanager listeners
  % to do this autoamtically.
  f = handle(ancestor(histPatch(1),'figure'));
  if ~isempty(f.findprop('linkPlot')) && f.linkPlot
      datamanager.updateLinkedGraphics(f);
  end
end

function a = ishistnumeric(b)
% for backward compatibility, logical is allowed in hist.m
a = isnumeric(b) || islogical(b);

function v = localGetVert(h,x,binCenters)

v = get(h,'vertices');
n = hist(x,binCenters);
v1 = [zeros(length(n),3) n(:) n(:)]';
v(:,2) = [v1(2:end)';0;0];

% ----------------------------
% Linked plot and data brushing callbacks

function localDataChange(h,data,bh,binCenters)

% Link behavior object DataSourceFcn callback. Responds to workspace events
% by redrawing linked graphics.

% Check that data does not have more than 1 column in which cases linking
% cannot work because each patch cannot be updated independently.
if min(size(data))>=2  || length(findobj(ancestor(h,'axes'),...
        '-and','-not',{'Behavior',struct},...
       '-function',@(h) ~isempty(hggetbehavior(h,'linked','-peek'))))>1   
   error('MATLAB:hist:InvalidLink', 'Cannot link histograms for non-vector data.')
end

bh.UserData = localMakeMap(data,binCenters);

if isscalar(binCenters)
   %  To compute the new vertices for potentially different bin positions,
   %  draw a histogram in an invisible figure and grab the veritices from
   %  the patch(es) created.
   f = figure('Visible','off','HandleVisibility','off');
   ax = axes('visible','off','parent',f);
   hist(ax,data,binCenters);
   hNew = findobj(ax,'type','patch');
   % Use the patch(es) in the invisible figure to reset the vertices in the
   % patches in the current figure.
   if isempty(hNew)
       return
   end
   set(h,'Vertices',get(hNew(1),'Vertices'));
   delete(f);
else
   set(h,'Vertices',localGetVert(h,data,binCenters));
end


function brushStruct = localBrushFunc(I,bobj,binCenters)

% Linked behavior object BrushFcn

% Converts variable brushing arrays (arrays of uint8 the same size as a 
% linked variable which define which subset of that variable is brushed and
% in which color) into generalized brushing data (the generalized form of 
% the BrushData property used by the brush behavior object). For histograms, 
% generalized brushing data is a struct with a field I representing the 
% height of each brushed bin and an index ColorIndex into the figure  
% BrushStyleMap representing the brushing color.

IdentRef = eye(length(binCenters),length(binCenters));
ind = find(I(:)>0);
if ~isempty(ind)
    BrushStyleMapInd = I(ind(1));
else
    BrushStyleMapInd = 0;
end
% Quick return if DataFcn has not yet been executed.
if isempty(I) || isempty(bobj.UserData)
    brushStruct = struct('I',[],'ColorIndex',0);
    return
end
brushBins = bobj.UserData.IndRowsToBins(I>0);
brushBins = brushBins(brushBins>0);
Iout = sum(IdentRef(:,brushBins),2);
brushStruct = struct('I',Iout,'ColorIndex',BrushStyleMapInd);

function bobj = localDrawFunc(brushStruct,b,gobj)

% Brush behavior object DrawFcn callback.

% The DrawFcn creates and draws brushing graphics from generalized brushing
% data (the generalized form of the BrushData property often created by
% the a link behavior object BrushFcn).

if ~feature('HGUsingMATLABClasses')
    fig = ancestor(gobj,'figure');
    bobj = b.UserData;
    if isempty(bobj)
        bobj = copyobj(gobj,get(gobj,'Parent'));
        b.UserData = bobj;
        set(bobj,'HandleVis','off');
    end
    if isempty(brushStruct) || isempty(brushStruct.I)
        set(bobj,'vis','off')
        return
    end
    % Find the brush color
    brushStyleMap = get(fig,'BrushStyleMap');
    if ~isempty(brushStruct.ColorIndex) && brushStruct.ColorIndex>0
        set(bobj,'FaceColor',...
            brushStyleMap(rem(brushStruct.ColorIndex-1,size(brushStyleMap,1))+1,:));
    end
    I = brushStruct.I; % Exclude brushStyleMap index
    v1 = [zeros(length(I),3) I(:) I(:)]';
    v = get(gobj,'vertices');
    v(:,2) = [v1(2:end)';0;0];
    set(bobj,'Vertices',v,'vis','on');
else
    if ~isempty(brushStruct)
        histlocalDrawFuncHGUsingMATLABClasses(...
            brushStruct.I(:)',brushStruct.ColorIndex,gobj);
    end
end

function I = localSetLinkedBrushing(L,region,gObj)

% Linked behavior object LinkBrushFcn callback. Converts region geometry
% into a variable brushing array (arrays of uint8 the same size as a 
% linked variable which define which subset of that variable is brush and 
% in which color). For MCOS graphics, region geometry is specified in figure
% pixel coordinates, otherwise it is specified in data coordinates.

% Respond to brushing gestures in linked plots. Returns brushing array
% for affected variables.
if ~feature('HGUsingMATLABClasses')
    ydata = get(gObj,'ydata');
    xdata = get(gObj,'xdata');
    xdata = xdata(2,:);
    ydata = ydata(2,:);
    if length(region)==4 % ROI brushing
         I1 = (ydata<=region(2)+region(4) & ydata>=region(2)) & ...
            (xdata<=region(1)+region(3) & xdata>=region(1));            
    elseif length(region)==2 % Vertex only, find closest
        if ~any(region(1)<xdata(2:end))
            ind = length(xdata);
        else
            [~,ind] = max(region(1)<xdata(2:end));
        end
        I1 = false(size(ydata));
        I1(ind) = true;
    elseif isempty(region) % Rectangular region
        I1 = false(size(ydata));
    end
else
    
    if ~isempty(region)
        % Get the figure coordinates of the left and right vertices of the
        % top of the non-zero histogram bars 
        vData = brushing.HistBrushing.histBarCameraCoords(gObj);
        vLeft = vData(:,1:2:end);
        vRight = vData(:,2:2:end);
        pixelvLeftLocations = brushing.select.transformCameraToFigCoord(gObj,vLeft);
        pixelvRightLocations = brushing.select.transformCameraToFigCoord(gObj,vRight);
        
        if length(region)==4 % ROI brushing
            % Find which vertices are in the bruishing ROI
            Ileft = brushing.select.inpolygon(region,pixelvLeftLocations);
            Iright = brushing.select.inpolygon(region,pixelvRightLocations);

            % Find which bars have been brushed
            ydata = get(gObj,'ydata');
            nonzerobars = find(any(ydata,1));
            I1 = nonzerobars(unique([Ileft(:);Iright(:)]));
        elseif length(region)==2
            % Find any non-zero bar where the left edge is to the left of the
            % clicked point and the right edge is to its right.
            Ileft = find(pixelvLeftLocations(1,:)<=region(1),1,'last');        
            if ~isempty(Ileft) && pixelvRightLocations(1,Ileft)>=region(1)
                ydata = get(gObj,'ydata');
                nonzerobars = find(any(ydata,1));
                I1 = nonzerobars(Ileft);
            else
                I1 = [];
            end
        end
    else
        I1 = false(size(get(gObj,'ydata')));
    end
end        
I2 = L.UserData.IndBinsToRows(:,I1);
I2(I2==0) = [];
I = I2;

function map = localMakeMap(x,binCenters)

% Build lookup tables for brushing histograms and displaying brushing
% annotations.
[x_,I] = sort(x);
[~,Iinv] = sort(I);
n = cumsum(hist(x_,binCenters));
IndRowsToBins = zeros(size(x));
IndRowsToBins(1:n(1)) = 1;
IBinsToRows = false(length(x),length(n));
IBinsToRows(1:n(1),1) = true; 
for k=2:length(n)
   IndRowsToBins(n(k-1)+1:n(k)) = k;
   IBinsToRows(n(k-1)+1:n(k),k) = true; 
end
IndRowsToBins = IndRowsToBins(Iinv);
IBinsToRows = IBinsToRows(Iinv,:);
IndBinsToRows = zeros(length(x),length(n));
for k=1:length(n)
    ind = find(IBinsToRows(:,k));
    IndBinsToRows(1:length(ind),k) = ind;
end
map = struct('VarData',x,'BinCenters',binCenters,'IndBinsToRows',IndBinsToRows,...
    'IndRowsToBins',IndRowsToBins);

% IndBinsToRows is a numDataPoints x numBins array where the k-th column
% contains the rows for the k-th bin.
% IndRowsToBins is a column numDataPoints x 1 identifying the bin number of
% each row.


% ----------------------------
% Data cursor callbacks
function datatipTxt = localHistDatatipCallback(~,evt,ctrs,counts,edges)
% Specify text displayed in the datatip next to the data cursor.

N_DIGITS = 3;
ind = get(evt,'DataIndex');

% Generate text to display.
datatipTxt = {
    ['Bin Count: ',num2str(counts(ind),N_DIGITS)],...
    '',...
    ['Bin Center: ',num2str(ctrs(ind),N_DIGITS)],...
    ['Bin Edges: [',num2str(edges(ind),N_DIGITS),', ',...
                    num2str(edges(ind+1),N_DIGITS),']'],...
    };


function localHistUpdateDataCursor(hDataCursor,target,hpatch)
% Specify data cursor position based on mouse click.

[~,~,patchInd,~,barface] = vertexpicker(hpatch,target,'-force');

% Specify index into bar series
% patchInd->barInd: 2,3,4,5->1, 7,8,9,10->2, 12,13,14,15->3, etc
barInd = floor((patchInd-1)/5)+1; % 4 patch vertices for each bar
verts = get(hpatch,'Vertices');
%numverts->numbars: 6->1, 11->2, 16->3, etc
numbars = floor(size(verts,1)/5);
if barInd>numbars,
    barInd = numbars;
elseif barInd==0
    barInd = 1;
end

% Specify bar face if we are mouse dragging off the bar plot
if isempty(barface)
    faces = get(hpatch,'Faces');
    [m,~] = find(faces==patchInd,1);
    if ~isempty(m)
        barface = verts(faces(m,:),:)';
    else % If we still can't find a face, take the first indexed face.
        barface = verts(faces(barInd,:),:)';
    end
end

% Specify cursor position
if ~isempty(barface)
   x_max = max(barface(1,:));
   x_min = min(barface(1,:));
   y_max = max(barface(2,:));
   y_min = min(barface(2,:));
   loc = [y_min y_max];
   [~, locInd] = max(abs(loc));
   set(hDataCursor,'Position',[(x_max+x_min)/2, loc(locInd), 0]);
   set(hDataCursor,'DataIndex',barInd);
   pos = hDataCursor.Position;
   set(hDataCursor,'TargetPoint',[pos(1) pos(2)]);
end


function localHistMoveDataCursor(hDataCursor,dir,hPatch)
% Specifies data cursor position when user selects arrows keys 
% (up,down,left,right).

currind = hDataCursor.DataIndex;

verts = get(hPatch,'Vertices');
if strcmp(dir,'up') || strcmp(dir,'right')
    %numverts->numbars: 6->1, 11->2, 16->3, etc
    numbars = floor(size(verts,1)/5);
    if currind < numbars
        currind = currind + 1;
    end
else
    if currind > 1
        currind = currind - 1;
    end
end
hDataCursor.DataIndex = currind;

vertind = (currind-1)*5+1;
binlo = verts(vertind+2,1);
binhi = verts(vertind+3,1);
bincenter = (binhi+binlo)/2;
barheight = verts(vertind+2,2);
currLine = [bincenter barheight];
currLine(~isfinite(currLine)) = 0;
hDataCursor.Position = currLine;

pos = hDataCursor.Position;
set(hDataCursor,'TargetPoint',[pos(1) pos(2)]);


